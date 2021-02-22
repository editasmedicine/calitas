package com.editasmedicine.aligner

import java.util.concurrent._

import com.editasmedicine.aligner.SearchReference.NoDictFastaFile
import com.editasmedicine.aligner.SequentialGuideAligner._
import com.editasmedicine.commons.LazyLogging
import com.editasmedicine.commons.clp.{ClpGroups, EditasTool}
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.{Cigar, CigarElem}
import com.fulcrumgenomics.commons.collection.BetterBufferedIterator
import com.fulcrumgenomics.fasta.ReferenceSequenceIterator
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, Metric, ProgressLogger, Sequences}
import com.fulcrumgenomics.vcf.api.{ArrayAttr, Variant, VcfSource}
import htsjdk.samtools.reference.{FastaSequenceFile, ReferenceSequence, ReferenceSequenceFileFactory}
import htsjdk.samtools.util.{Interval, SequenceUtil, StringUtil}
import htsjdk.samtools.{CigarOperator, SAMSequenceDictionary}
import com.fulcrumgenomics.commons.util.{StringUtil => Strings}
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor

import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import scala.math.min

/**
  * Companion object for SearchReference that provides various functions for windowing
  * through a reference sequence in parallel.
  */
object SearchReference extends LazyLogging {
  /** Represents a window on a reference chromosome along with those bases.
    * Start & End are 1-based closed-ended coordinates. */
  case class RefWindow(chrom: String, start: Int, end: Int, bases: Array[Byte])

  /** Trivial class so we can drop the sequence dictionary as fast as possible. */
  private[aligner] class NoDictFastaFile(ref: PathToFasta) extends FastaSequenceFile(ref, false) {
    override def findAndLoadSequenceDictionary(ref: PathToFasta): SAMSequenceDictionary = null
  }

  /** Generates an iterator over windows on the reference. */
  private[aligner] def windowIterator(fasta: PathToFasta, windowSize: Int, stepSize: Int, chrom: Option[String]): Iterator[RefWindow] = {
    val empty = Array[Byte](0)
    val chromIter = chrom match {
      case Some(chr) => Iterator(ReferenceSequenceFileFactory.getReferenceSequenceFile(fasta).getSequence(chr))
      case None      => ReferenceSequenceIterator(new NoDictFastaFile(fasta))
    }

    chromIter.flatMap { chr =>
      val name  = chr.getName
      val len   = chr.length()
      val bases = chr.getBases

      // All code here is zero-base open ended
      Range(0, len-1, step=stepSize).iterator.map { start =>
        val end = min(len, start + windowSize)

        // Try and avoid large blocks of Ns
        var adjustedStart = start
        var adjustedEnd   = end
        while (adjustedStart < adjustedEnd && bases(adjustedStart) == 'N') adjustedStart += 1
        while (adjustedStart < adjustedEnd && bases(adjustedEnd-1) == 'N') adjustedEnd   -= 1

        val thisWindowLength = adjustedEnd - adjustedStart
        val windowBases      = if (thisWindowLength <= 0) empty else {
          val tmp = new Array[Byte](thisWindowLength)
          System.arraycopy(bases, adjustedStart, tmp, 0, thisWindowLength)
          tmp
        }
        StringUtil.toUpperCase(windowBases)
        RefWindow(name, adjustedStart + 1, adjustedEnd, windowBases)
      }
    }
  }


  /** Generates an executor for submitting alignments jobs to. */
  private def executor(threads: Int): ExecutorService = {
    val queue   = new ArrayBlockingQueue[Runnable](threads * 500)
    val factory = new ThreadFactory() {
      private var index = 1
      override def newThread(r: Runnable): Thread = index.synchronized {
        val t = new Thread(r, s"aligner-pool-thread-${index}")
        t.setDaemon(true)
        index += 1
        t
      }
    }

    new ThreadPoolExecutor(threads, threads, 1, TimeUnit.MINUTES, queue, factory) {
      /** Override to wait until queue has space. */
      override def execute(command: Runnable): Unit = {
        while (queue.remainingCapacity() == 0) Thread.sleep(25)
        super.execute(command)
      }
    }
  }
}


/**
  * Object that contains functions for searching a reference with variants specified by a VCF.
  */
object SearchReferencesWithVariants extends LazyLogging {
  /** Represents a variant allele that has been inserted into the reference sequence.
    * Position is 1-based.
    */
  case class VariantAllele(id: String, pos: Int, ref: String, alt: String, af: Float) {
    def displayString: String = {
      // Note the conversion to zero-based pos for display
      f"${if (id.nonEmpty) id else "."}:${pos-1}:$ref>$alt:${af}%.3f"
    }
  }

  /**
    * Represents a window on the reference that has been modified by the insertion of one or more
    * variant alleles.
    *
    * Start is 1-based.
    */
  case class VariantWindow(chrom: String,
                           start: Int,
                           variants: Seq[VariantAllele],
                           cigar: Cigar,
                           bases: Array[Byte]
                          ) {

    def length: Int = bases.length

    /** Returns the 0-based offset into the reference sequence of the given base within the window.
      *
      * @param offset the 0-based offset of the base within the variant window
      * @param preceding if true inserted bases are reported as the preceding reference base, else the succeeding
      * @return the 0-based offset within the reference genome
      */
    def refOffsetAtBaseOffset(offset: Int, preceding: Boolean): Int = {
      if (offset == bases.length) {
        start - 1 + cigar.lengthOnTarget
      }
        else {
        var refOffset   = start - 1
        var baseOffset  = 0
        val elems       = cigar.iterator
        var currentElem = elems.next()

        // Fast forward to find the cigar element that the offset is within
        while (offset >= baseOffset + currentElem.lengthOnQuery) {
          refOffset  += currentElem.lengthOnTarget
          baseOffset += currentElem.lengthOnQuery
          currentElem = elems.next()
        }

        currentElem.operator match {
          case CigarOperator.I => if (preceding) refOffset - 1 else refOffset
          case CigarOperator.M => refOffset + (offset - baseOffset)
          case op              => unreachable(s"Query bases can't be present at operator $op.")
        }
      }
    }
  }


  /**
    * Class that encapsulates a set of variants and a chosen allele per variant. The alleles are specified
    * as the index of the allele within the respective variant.  I.e. 0 == ref allele, 1 = first ALT, etc.
    *
    * Variants are directly from VCF and pos/end are 1-based closed ended.
    */
  case class VariantSet(variants: IndexedSeq[Variant], alleles: IndexedSeq[Int]) {
    require(variants.size == alleles.size)
    require(alleles.forall(_ > 0)) // no ref alleles

    /** The number of variants in the set. */
    def size: Int = variants.size

    /** The set of indices to index into variants and alleles. */
    def indices: Seq[Int] = Range(0, size)

    def start: Int = variants.head.pos

    def end: Int = variants.last.end

    /** Returns true if there are no conflicts between the set of alleles specified - i.e. no two variants
      * attempt to modify the same reference base. */
    def isValid: Boolean = {
      if (variants.size == 1) true else {
        val spans = indices.filter(i => alleles(i) != 0).map { idx =>
          val variant   = variants(idx)
          val refLength = variant.alleles.ref.length
          new Interval(variant.chrom, variant.pos, variant.pos + refLength - 1)
        }

        val hasOverlaps = Range(0, spans.length - 1).exists(i => spans(i).overlaps(spans(i+1)))
        !hasOverlaps
      }
    }

    /** Returns a VariantAllele object for the ith entry in the set. */
    def variantAllele(i: Int) = {
      val v  = variants(i)
      val a  = alleles(i)
      val af = v.get[ArrayAttr[Float]]("AF").map(_(a-1)).getOrElse(0.0F)
      VariantAllele(v.id.getOrElse(""), v.pos, v.alleles.ref.bases, v.alleles(a).value, af)
    }
  }

  /**
    * Generates an iterator over a set of windows in the reference that contain one or more variants each.
    *
    * @param fasta the fasta file for the reference genome
    * @param vcf the VCF of variants to insert into the reference
    * @param chrom an optional chromosome to restrict analysis to
    * @param padding how much sequence to include upstream and downstream of variants and also the
    *                maximum distance between variants
    * @param maxVariants the maximum number of variants to consider together.  I.e. if there are more than maxVariants
    *                    all within padding distance of one another, treat the variants individually and do not attempt
    *                    to construct all allele permutations.
    * @return an iterator over a set of sequence windows including variants
    */
  private[aligner] def variantWindowIterator(fasta: PathToFasta,
                                             vcf: PathToVcf,
                                             chrom: Option[String],
                                             padding: Int,
                                             maxVariants: Int): Iterator[VariantWindow] = {
    val chromIter: Iterator[ReferenceSequence] = (chrom match {
      case Some(chr) => Iterator(ReferenceSequenceFileFactory.getReferenceSequenceFile(fasta).getSequence(chr))
      case None      => ReferenceSequenceIterator(new NoDictFastaFile(fasta))
    }).map { r => SequenceUtil.upperCase(r.getBases); r }

    val vcfSource = VcfSource(vcf)
    val vcfIter = chrom match {
      case Some(chr) => vcfSource.query(chr, 1, Int.MaxValue).bufferBetter
      case None      => vcfSource.iterator.bufferBetter
    }

    new Iterator[VariantWindow] {
      private var cache = new mutable.Queue[VariantWindow]()
      private var ref: ReferenceSequence = chromIter.next()

      override def hasNext: Boolean = {
        while (cache.isEmpty && vcfIter.hasNext) fillCache()
        cache.nonEmpty
      }

      override def next(): VariantWindow = {
        if (!hasNext) throw new NoSuchElementException("next() called on empty iterator.")
        cache.dequeue()
      }

      def fillCache(): Unit = {
        val vs         = nextChunk(vcfIter, padding)
        val chunks     = reChunk(vs, padding)
        val alleleSets = chunks.flatMap(c => alleleCombos(c, maxVariants))
        while (vs.head.chrom != ref.getName) ref = chromIter.next()  // Advance the reference if needed

        alleleSets.foreach { set => this.cache += buildVariantWindow(set, ref, padding) }
      }
    }
  }


  /** Builds up a VariantWindow object from an ordered set of variants and the appropriate chromosome.
    * Will generate an initial sequence that spans the variant set +/- padding, and then go through
    * and substitute in alternate alleles into the sequence.
    */
  private[aligner] def buildVariantWindow(set: VariantSet, ref: ReferenceSequence, padding: Int): VariantWindow = {
    val windowStart = math.max(1, set.start - padding)
    val windowEnd   = math.min(ref.length(), set.end + padding)
    var bases       = ref.getBases.slice(windowStart - 1, windowEnd)
    val alleles     = set.indices.map(set.variantAllele)

    // Modify the sequence from the end working backwards to keep the math simple
    alleles.reverseIterator.foreach { allele =>
      val startIndex = allele.pos - windowStart

      if (allele.ref.length == allele.alt.length) {
        forloop(from = 0, until = allele.ref.length) { i => bases(startIndex + i) = allele.alt(i).toByte }
      }
      else {
        bases = bases.patch(startIndex, allele.alt.getBytes, allele.ref.length)
      }
    }

    // Build up the CIGAR that maps from the reference into the base string
    val cigar = {
      var refPos     = windowStart
      var baseOffset = 0
      val elems      = IndexedSeq.newBuilder[CigarElem]

      alleles.foreach { allele =>
        val precedingMatch = allele.pos - refPos

        if (precedingMatch > 0) {
          elems      += CigarElem(CigarOperator.M, precedingMatch)
          refPos     += precedingMatch
          baseOffset += precedingMatch
        }

        if (allele.ref.length == allele.alt.length) { // substitution
          elems += CigarElem(CigarOperator.M, allele.ref.length)
        }
        else if (allele.ref.length == 1 && allele.alt.length > 1) { // simple insertion
          elems += CigarElem(CigarOperator.M, 1)
          elems += CigarElem(CigarOperator.I, allele.alt.length - 1)
        }
        else if (allele.ref.length > 1 && allele.alt.length == 1) { // simple deletion
          elems += CigarElem(CigarOperator.M, 1)
          elems += CigarElem(CigarOperator.D, allele.ref.length - 1)
        }
        else { // complicated substitution with length difference
          elems += CigarElem(CigarOperator.D, allele.ref.length)
          elems += CigarElem(CigarOperator.I, allele.alt.length)
        }

        refPos  += allele.ref.length
        baseOffset += allele.alt.length
      }

      // Add the last matching block after the last alt allele
      elems += CigarElem(CigarOperator.M, bases.length - baseOffset)
      Cigar(elems.result()).coalesce
    }

    require(cigar.lengthOnQuery == bases.length, s"Cigar: $cigar, LoQ: ${cigar.lengthOnQuery}, len(bases): ${bases.length}")
    VariantWindow(ref.getName, windowStart, alleles, cigar, bases)
  }

  /** Takes the next chunk of variants that are close together. */
  private def nextChunk(vs: BetterBufferedIterator[Variant], maxDistance: Int): IndexedSeq[Variant] = {
    var last    = vs.next()
    val builder = IndexedSeq.newBuilder[Variant]
    builder += last

    while (vs.hasNext && vs.head.chrom == last.chrom && vs.head.pos <= last.end + maxDistance) {
      last = vs.next()
      builder += last
    }

    builder.result()
  }

  /** Takes a contiguous (and hopefully relatively small) chunk of variants and re-chunks it into chunks
    * where the last variant in each chunk is no more than maxDistance from the end of the first variant
    * in the chunk.
    */
  private def reChunk(vs: IndexedSeq[Variant], maxDistance: Int): IndexedSeq[IndexedSeq[Variant]] = {
    vs.tails.filter(_.nonEmpty).map { sub =>
      sub.takeWhile(_.pos - sub.head.end <= maxDistance)
    }.toIndexedSeq
  }

  /** Generates a sequence of VariantSets where each set represents a chosen set of alleles
    * for the input set of variants. */
  private[aligner] def alleleCombos(vs: IndexedSeq[Variant], maxVariants: Int): IndexedSeq[VariantSet] = {
    if (vs.size > maxVariants) {
      logger.warning(s"Not checking combos for ${vs.size} variants at ${vs.head.chrom}:${vs.head.pos}-${vs.last.end}")
      val v = vs.head
      v.alleles.alts.indices.map(a => VariantSet(IndexedSeq(v), IndexedSeq(a+1)))
    }
    else {
      // Generate the combos and then remove any reference alleles
      alleleCombos(vs.map(_.alleles.size))
        .iterator
        .map { alleles =>
          val (nonRefVs, nonRefAs) = vs.zip(alleles).filter(_._2 != 0).unzip
          VariantSet(nonRefVs, nonRefAs)
        }
        .filter(_.variants.nonEmpty)
        .filter(_.isValid)
        .toIndexedSeq
    }
  }

  /**
    * Generates an array of arrays, where each sub-array is a valid combination of alleles.  E.g for a pair
    * of bi-allelic variants it will produce [[0, 0], [0, 1], [1, 0], [1, 1]].
    *
    * @param alleleCounts the number of alleles for each of a set of variants
    */
  private[aligner] def alleleCombos(alleleCounts: Seq[Int]): Array[Array[Int]] = {
    val results = Array.fill(alleleCounts.product)(new Array[Int](alleleCounts.length))
    var denom   = 1

    forloop (from=0, until=alleleCounts.length) { i =>
      val n = alleleCounts(i)
      denom *= n
      val groupSize = results.length / denom

      var j      = 0
      var allele = 0
      while (j < results.length) {
        val end = j + groupSize
        while (j < end) {
          results(j)(i) = allele
          j += 1
        }
        allele = (allele + 1) % n
      }
    }

    results
  }
}

@clp(group=ClpGroups.Alignment, description=
  """
    |Searches a reference sequence for alignments of a guide+PAM.
    |
    |The search is performed in a sequential manner, first finding all candidate alignments of the guide without
    |the PAM, and then extending those alignments to include an optional PAM.  The search may be performed without
    |a PAM (i.e. PAMless), with a single PAM or with multiple PAMs.  When multiple PAMs are provided the first PAM
    |must be provided as part of the guide, with subsequent PAM being provided via `--auxillary-pams`.  When extending
    |alignments each PAM is considered, and the extension with the best PAM alignment is retained.  Best is defined by
    |maximizing score; when scores are equal PAMs earlier in the list will be preferred.  All PAM sequences must be
    |specified in lower case, while protospacer sequence must be upper case.  E.g.:
    |
    |    --guide ATCGATCGATAGACTGCATnrg --auxiliary-pams nnrg kgg
    |
    |The scoring system is, by default, setup to guarantee that all alignments within the thresholds defined by
    |a) `--max-guide-diffs`, b) `--max-pam-mismatches` and c) `--max-gaps-between-guide-and-pam`, will be discovered
    |for reasonable values of those parameters and for query sequences of 20-40bp (i.e. common protospacer + PAM
    |lengths).  Great care must be taken when adjusting scoring parameters in order not to void this guarantee. Scoring
    |is controlled via four major parameters:
    |
    |1. `--guide-mismatch-net-cost`: the _net_ cost of converting a match to a mismatch within the alignment.
    |2. `--pam-mismatch-net-cost`: as above, but for mismatches within the PAM
    |3. `--genome-gap-net-cost`: the cost of any 1bp gap in the genome (a.k.a. a 1bp bulge in the guide)
    |4. `--guide-gap-net-cost`: the cost of any 1bp gap in the guide (a.k.a. a 1bp bulge in the genome)
    |
    |For longer gaps the cost is simply gap length multiplied by the appropriate net cost.  The 'net' in each
    |name refers to the fact that unlike traditional aligner scores these scores factor in the possible loss of a
    |match in addition to the introduction of a mismatch or gap.  The default values of these cost parameters are set
    |such that mismatches are preferred to gaps in the guide which are in turn preferred to gaps in the genome. In
    |addition mismatches in the PAM, by default, are slightly more than twice as expensive as mismatches in the guide.
    |If changes to these parameters are desired it is important that the following be true:
    |
    |```
    |min_cost = min(guide_mismatch_net_cost, guide_gap_net_cost, genome_gap_net_cost)
    |max_cost = max(guide_mismatch_net_cost, guide_gap_net_cost, genome_gap_net_cost)
    |(max-guide-diffs + 1) * min_cost > max-guide-diffs * max_cost
    |```
    |
    |If the constraint above is violated then the aligner may prefer alignments with too many differences, and then
    |filter those out and as a result not report alignments that are within the bounds specified!
    |
    |Care should be taken when setting the limits on mismatches and gaps.  Different results will be obtained by
    |running with, e.g. `--max-guide-diffs=5 --max-pam-mismatches=1` and the post-filtering vs. running directly with
    |`--max-guide-diffs=3 --max-pam-mismatches=1`.  For example, with the former settings the aligner may prefer to emit
    |an alignment with 4 mismatches in the guide and 0 mismatches in the PAM in preference to a competing alignment with
    |3 mismatches in the guide and 1 mismatch in the PAM.  The latter settings will, on the other hand, emit the
    |3+1 mismatch alignment.  These concerns can be alleviated by setting `--max-overlap` to some value much larger than
    |the guide length (e.g. 100), causing all overlapping alignments to be emitted.
  """)
class SearchReference
(@arg(flag='i', doc="Guide with PAM, PAM must be lower case.") val guide: String,
 @arg(flag='I', doc="ID of the guide.") val guideId: String,
 @arg(flag='x', doc="Additional PAM sequences. Must be lower case.", minElements=0) val auxiliaryPams: Seq[String] = Seq.empty,
 @arg(flag='r', doc="Reference genome fasta.") val ref: PathToFasta,
 @arg(flag='v', doc="Optional VCF of variants to merge into the genome.") val variants: Option[PathToVcf] = None,
 @arg(flag='V', doc="Exclude clusters of this more than this many variants.") val maxVariants: Int = Defaults.MaxVariantsInCluster,
 @arg(flag='o', doc="Output file to write.") val output: FilePath = Io.StdOut,
 @arg(flag='t', doc="Threads to use for alignments.") val threads: Int = 8,
 @arg(flag='w', doc="Window size to align to.") val windowSize: Int = 1000,
 @arg(flag='d', doc="Maximum number of differences (mms+gaps) between guide and genome.") val maxGuideDiffs: Int = Defaults.MaxGuideDiffs,
 @arg(flag='p', doc="Maximum mismatches in the PAM.") val maxPamMismatches: Int = Defaults.MaxPamMismatches,
 @arg(flag='g', doc="Maximum gap bases between guide and PAM") val maxGapsBetweenGuideAndPam: Int = Defaults.MaxGapsBetweenGuideAndPam,
 @arg(flag='D', doc="Maximum total diffs in alignments.") val maxTotalDiffs: Option[Int] = None,
 @arg(flag='O', doc="Maximum overlap allowed between alignments on the same strand.") val maxOverlap: Int = Defaults.MaxOverlap,
 @arg(flag='m', doc="Net cost of going from a match to a mismatch in the guide.") val guideMismatchNetCost: Int = Defaults.MismatchNetCost,
 @arg(flag='M', doc="Net cost of going from a match to a mismatch in the PAM.") val pamMismatchNetCost: Int = Defaults.PamMismatchNetCost,
 @arg(flag='b', doc="Net cost of a 1bp gap in the genome.") val genomeGapNetCost: Int = Defaults.GenomeGapNetCost,
 @arg(flag='B', doc="Net cost of a 1bp gap in the guide.") val guideGapNetCost: Int = Defaults.GuideGapNetCost,
 @arg(flag='c', doc="Examine only the named chromosome.") val chrom: Option[String] = None,
) extends EditasTool {
  import SearchReference._

  Io.assertReadable(ref)
  Io.assertCanWriteFile(output)

  private val dict = SAMSequenceDictionaryExtractor.extractDictionary(ref)

  private val aligner  = new SequentialGuideAligner(
    mismatchNetCost    = guideMismatchNetCost,
    pamMismatchNetCost = pamMismatchNetCost,
    genomeGapNetCost  = genomeGapNetCost,
    guideGapNetCost = guideGapNetCost
  )

  private val maxTotalDiffsActual = this.maxTotalDiffs.getOrElse(this.maxGuideDiffs + this.maxGapsBetweenGuideAndPam + this.maxPamMismatches)

  /** Generates a string of the major parameters that is output with every hit. */
  private def coreParameters: String = Map(
    Strings.camelToGnu("maxVariants")               -> maxVariants,
    Strings.camelToGnu("windowSize")                -> windowSize,
    Strings.camelToGnu("maxGuideDiffs")             -> maxGuideDiffs,
    Strings.camelToGnu("maxPamMismatches")          -> maxPamMismatches,
    Strings.camelToGnu("maxGapsBetweenGuideAndPam") -> maxGapsBetweenGuideAndPam,
    Strings.camelToGnu("maxTotalDiffs")             -> maxTotalDiffsActual,
    Strings.camelToGnu("maxOverlap")                -> maxOverlap,
    Strings.camelToGnu("guideMismatchNetCost")      -> guideMismatchNetCost,
    Strings.camelToGnu("pamMismatchNetCost")        -> pamMismatchNetCost,
    Strings.camelToGnu("genomeGapNetCost")          -> genomeGapNetCost,
    Strings.camelToGnu("guideGapNetCost")           -> guideGapNetCost,
  ).map { case (k,v) => s"$k=$v"}.toSeq.sorted.mkString(";")

  // Construct the Guide object and fail if anything isn't valid
  val query = Guide(guide, auxiliaryPams)

  override def execute(): Unit = {
    val hits       = new ArrayBuffer[ReferenceHit](50000)
    val pool       = executor(threads)
    val hitBuilder = ReferenceHit.Builder(
      guideId   = guideId,
      guide     = query,
      ref       = ref,
      vcf       = variants,
      arguments = coreParameters,
      alignerId = "CALITAS:SearchReference")

    ////////////////////////////////////////////////////////////////////////////
    // Align to the unaltered reference
    ////////////////////////////////////////////////////////////////////////////
    {
      val guideLength   = guide.length
      val windowOverlap = guide.length + maxGuideDiffs + maxGapsBetweenGuideAndPam - 1
      val stepSize      = windowSize - windowOverlap
      val progress      = ProgressLogger(logger=logger, noun="windows", verb="Processed", unit=25000)
      val iter          = windowIterator(ref, windowSize=this.windowSize, stepSize=stepSize, chrom=chrom)

      logger.info("Aligning to reference genome without variants.")
      iter
        .filter(_.bases.length >= guideLength)
        .foreach { window =>
          pool.execute(() => {
            try {
              val results = aligner.align(
                guide                     = query,
                target                    = window.bases,
                targetName                = window.chrom,
                targetOffset              = window.start - 1,
                maxGuideDiffs             = maxGuideDiffs,
                maxPamDiffs               = maxPamMismatches,
                maxGapsBetweenGuideAndPam = maxGapsBetweenGuideAndPam,
                maxTotalDiffs             = maxTotalDiffsActual,
                maxOverlap                = maxOverlap
              )

              results.iterator.map(a => hitBuilder.build(a)).foreach { hit => hits.synchronized(hits += hit) }
              progress.record(window.chrom, window.start)
            }
            catch {
              case ex: Throwable => logger.error(s"Encountered an exception: $ex")
            }
          })
        }

      logger.info("Reference windows processed.")
    }


    ////////////////////////////////////////////////////////////////////////////
    // Align to the reference with variants
    ////////////////////////////////////////////////////////////////////////////
    this.variants.foreach { vcf =>
      val iter = SearchReferencesWithVariants.variantWindowIterator(
        fasta       = ref,
        vcf         = vcf,
        chrom       = this.chrom,
        padding     = query.length - 1 + maxGuideDiffs + maxGapsBetweenGuideAndPam,
        maxVariants = maxVariants
      )

      val progress = ProgressLogger(logger=logger, noun="variant windows", verb="Processed", unit=100000)

      iter.foreach { window =>
        pool.execute(() => {
          try {
            val relativeResults = aligner.align(
              guide                     = query,
              target                    = window.bases,
              targetName                = window.chrom,
              targetOffset              = 0,
              maxGuideDiffs             = maxGuideDiffs,
              maxPamDiffs               = maxPamMismatches,
              maxGapsBetweenGuideAndPam = maxGapsBetweenGuideAndPam,
              maxTotalDiffs             = maxTotalDiffsActual,
              maxOverlap                = maxOverlap
            )

            // If possible insert the flanking bases here because its possible that the alignment starts or
            // ends in an insertion and we should use the remainder of the insertion for the flank
            val flankedResults = relativeResults.map { a =>
              val left10  = if (a.guideStartOffset < 10) None else Some(new String(window.bases.slice(a.guideStartOffset-10, a.guideStartOffset)))
              val right10 = if (window.length - a.guideEndOffset < 10) None else Some(new String(window.bases.slice(a.guideEndOffset, a.guideEndOffset+10)))
              val left8   = if (a.startOffset < 8) None else Some(new String(window.bases.slice(a.startOffset-8, a.startOffset)))
              val right8  = if (window.length - a.endOffset < 8) None else Some(new String(window.bases.slice(a.endOffset, a.endOffset+8)))

              if (a.isPositiveStrand) a.copy(leftOfGuide10bp=left10, rightOfGuide10bp=right10, leftOfFullAln8bp=left8, rightOfFullAln8bp=right8)
              else {
                a.copy(
                  leftOfGuide10bp   = right10.map(Sequences.revcomp),
                  rightOfGuide10bp  = left10.map(Sequences.revcomp),
                  leftOfFullAln8bp  = right8.map(Sequences.revcomp),
                  rightOfFullAln8bp = left8.map(Sequences.revcomp)
                )
              }
            }

            val absoluteResults = flankedResults.map(a => a.copy(
              startOffset      = window.refOffsetAtBaseOffset(a.startOffset, preceding = true),
              endOffset        = window.refOffsetAtBaseOffset(a.endOffset,   preceding = false),
              guideStartOffset = window.refOffsetAtBaseOffset(a.guideStartOffset, preceding = true),
              guideEndOffset   = window.refOffsetAtBaseOffset(a.guideEndOffset,   preceding = false),
            ))

            absoluteResults.iterator.map(a => hitBuilder.build(a, window.variants)).foreach { hit => hits.synchronized(hits += hit) }
            progress.record(window.chrom, window.start)
          }
          catch {
            case ex: Throwable => logger.error(s"Encountered an exception: $ex")
          }
        })
      }
    }

    logger.info("Variant windows processed.")

    pool.shutdown()
    pool.awaitTermination(10, TimeUnit.MINUTES)

    ////////////////////////////////////////////////////////////////////////////
    // Sort, filter and output the alignments
    ////////////////////////////////////////////////////////////////////////////
    logger.info("Sorting and Outputting.")
    val keepers = removeOverlaps(hits, maxOverlap)

    val fwdFraction = keepers.count(_.strand == "+") / keepers.size.toDouble
    if (fwdFraction > 0.52 || fwdFraction < 0.48) logger.warning(f"Strand imbalance: $fwdFraction%2f of alignments are on the F strand.")

    val writer = Metric.writer[ReferenceHit](output)
    writer ++= ReferenceHit.sort(keepers, Left(this.dict))
    writer.close()
  }


  /** Removes overlapping alignments that overlap by at least `maxOverlap` and have lower scores. */
  def removeOverlaps(hits: Seq[ReferenceHit], maxOverlap: Int): Seq[ReferenceHit] = {
    val keepers  = new mutable.ArrayBuffer[ReferenceHit](hits.size)

    hits.groupBy(h => s"{${h.chromosome}:${h.strand}:${h.variant_description.getOrElse("")}").foreach { case (_, hs) =>
      val sorted = ReferenceHit.sort(hs, Left(this.dict))
      val iter   = sorted.iterator.bufferBetter

      while (iter.hasNext) {
        val hit = iter.next()

        // Discard any alignments overlapping too much whose scores are lower
        while (iter.hasNext && iter.head.overlap(hit) >= maxOverlap && iter.head.score <= hit.score) {
          iter.next()
        }

        // Only output if there isn't an overlapping alignment next; because of the loop above we can
        // infer that if an overlapping alignment remains, it must be better
        if (!iter.hasNext || iter.head.overlap(hit) < maxOverlap) keepers += hit
      }
    }

    keepers.toSeq
  }
}
