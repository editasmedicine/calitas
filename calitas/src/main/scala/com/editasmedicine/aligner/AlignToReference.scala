package com.editasmedicine.aligner

import com.editasmedicine.aligner.SequentialGuideAligner.{Defaults, Guide}
import com.editasmedicine.commons.clp.{ClpGroups, EditasTool}
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.util.DelimitedDataParser
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, Metric}
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import com.fulcrumgenomics.commons.util.{StringUtil => Strings}

@clp(group=ClpGroups.Alignment, description=
  """
    |Performs glocal alignment of query sequence to a window on the reference. Input should be
    |a tab-delimited file with the following columns (with headers):
    |  - id: the ID of the query sequence; optional - if not present the sequence itself is used
    |  - query: the query sequence to be aligned
    |  - chrom: the chromosome to align to
    |  - position: the center of the window to align the query to
    |
    |The query sequence may be a mixture of upper and lower case sequence.  The lower case sequence
    |(generally used for the PAM sequence) is scored differently than upper-case sequence, with
    |higher match and mismatch values, to make mismatches in lower-case regions less common.
    |
    |If all three of `--max-guide-diffs`, `--max-pam-mismatches` and `--max-overlap` are specified then
    |all alignments that meet the given criteria for a query will be emitted.  If none of the three parameters
    |are provided, the single best alignment of each query will be reported.  If some but not all three parameters
    |are specified then an error will be generated.
    |
    |The `--window-size` parameter can be used to control the size of the window (centered on the position given for
    |each query) that the queries will be aligned to.  E.g. if a window size of 60bp is given, the target sequence is
    |position +/- 30bp.  If `--window-size` is not given, it will be defaulted to 2 * length(query) for each query.
  """)
class AlignToReference
(@arg(flag='i', doc="Input file of sequence queries and approximate positions.") val input: FilePath,
 @arg(flag='r', doc="Reference genome fasta, must be indexed with faidx.") val ref: PathToFasta,
 @arg(flag='o', doc="Output file to write.") val output: FilePath = Io.StdOut,
 @arg(flag='w', doc="Window size to align to.") val windowSize: Option[Int] = None,
 @arg(flag='d', doc="Maximum number of differences (mms+gaps) between guide and genome.") val maxGuideDiffs: Option[Int] = None,
 @arg(flag='p', doc="Maximum mismatches in the PAM.") val maxPamMismatches: Option[Int] = None,
 @arg(flag='g', doc="Maximum gap bases between guide and PAM") val maxGapsBetweenGuideAndPam: Int = Defaults.MaxGapsBetweenGuideAndPam,
 @arg(flag='D', doc="Maximum total diffs in alignments.") val maxTotalDiffs: Option[Int] = None,
 @arg(flag='O', doc="Maximum overlap allowed between alignments on the same strand.") val maxOverlap: Option[Int] = None,
 // Scoring parameters
 @arg(flag='m', doc="Net cost of going from a match to a mismatch in the guide.") val guideMismatchNetCost: Int = Defaults.MismatchNetCost,
 @arg(flag='M', doc="Net cost of going from a match to a mismatch in the PAM.") val pamMismatchNetCost: Int = Defaults.PamMismatchNetCost,
 @arg(flag='b', doc="Net cost of a 1bp gap in the genome.") val genomeGapNetCost: Int = Defaults.GenomeGapNetCost,
 @arg(flag='B', doc="Net cost of a 1bp gap in the guide.") val guideGapNetCost: Int = Defaults.GuideGapNetCost,

 @arg(flag='t', doc="Threads to use for alignments.") val threads: Int = 8
) extends EditasTool {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  case class Task(id: String, query: String, chrom: String, pos: Int)

  // Setup the aligner
  private val refFile  = ReferenceSequenceFileFactory.getReferenceSequenceFile(ref)
  validate(refFile.isIndexed, s"Reference genome must have a fasta index: $ref")
  validate(refFile.getSequenceDictionary != null && !refFile.getSequenceDictionary.isEmpty,
    s"Reference genome must have a sequence dictionary: $ref")

  private val aligner  = new SequentialGuideAligner(
    refFile            = Some(refFile),
    mismatchNetCost    = guideMismatchNetCost,
    pamMismatchNetCost = pamMismatchNetCost,
    genomeGapNetCost  = genomeGapNetCost,
    guideGapNetCost = guideGapNetCost
  )

  // Dummy guide that is used to initialize the ReferenceHit.Builder but is never really used
  private val dummyGuide = Guide("AAAnnn")

  // The set of parameters that are output in the reference hits
  /** Generates a string of the major parameters that is output with every hit. */
  private val coreParameters: String = Map(
    Strings.camelToGnu("maxGuideDiffs")             -> maxGuideDiffs,
    Strings.camelToGnu("maxPamMismatches")          -> maxPamMismatches,
    Strings.camelToGnu("maxGapsBetweenGuideAndPam") -> maxGapsBetweenGuideAndPam,
    Strings.camelToGnu("maxOverlap")                -> maxOverlap,
    Strings.camelToGnu("guideMismatchNetCost")      -> guideMismatchNetCost,
    Strings.camelToGnu("pamMismatchNetCost")        -> pamMismatchNetCost,
    Strings.camelToGnu("genomeGapNetCost")          -> genomeGapNetCost,
    Strings.camelToGnu("guideGapNetCost")           -> guideGapNetCost,
  ).map { case (k,v) => s"$k=$v"}.toSeq.sorted.mkString(";")

  (maxGuideDiffs, maxPamMismatches, maxOverlap) match {
    case (Some(_), Some(_), Some(_)) => logger.info("Will output all alignments matching given parameters.")
    case (None,    None,    None)    => logger.info("Will output the single best alignment for each query.")
    case _                           => invalid("Must specify all or none of: --max-guide-diffs, --max-pam-mismatches, --max-overlap")
  }

  /** Main method that is executed when the tool is invoked. */
  override def execute(): Unit = {
    logger.info("Reading input file.")
    val parser = DelimitedDataParser(input, '\t')
    val tasks = parser.map { row =>
      val query = row.string("query")
      val id    = row.get[String]("id", allowMissingColumn=true).getOrElse(query)
      Task(id=id, query=query, chrom=row.string("chrom"), pos=row[Int]("position"))
    }

    logger.info("Beginning alignments.")
    val out     = Metric.writer[ReferenceHit](output)
    val builder = ReferenceHit.Builder(guideId="n/a", guide=dummyGuide, ref=ref, vcf=None,
      alignerId="CALITAS:AlignToReference", arguments=coreParameters)

    while (tasks.hasNext) {
      val batch = tasks.take(10000).toIndexedSeq
      val results = batch.parWith(threads).flatMap { task =>
        val guide  = Guide(task.query)

        val alns = (maxGuideDiffs, maxPamMismatches, maxOverlap) match {
          case (Some(mgd), Some(mpm), Some(mo)) =>
            aligner.alignToRef(
              guide                     = guide,
              chrom                     = task.chrom,
              pos                       = task.pos,
              windowSize                = windowSize,
              maxGuideDiffs             = mgd,
              maxGapsBetweenGuideAndPam = maxGapsBetweenGuideAndPam,
              maxPamDiffs               = mpm,
              maxTotalDiffs             = maxTotalDiffs.getOrElse(mgd + maxGapsBetweenGuideAndPam + mpm),
              maxOverlap                = mo
        )
          case _ =>
            Seq(aligner.alignToRefBest(
              guide = guide,
              chrom                     = task.chrom,
              pos                       = task.pos,
              windowSize                = windowSize,
              maxGapsBetweenGuideAndPam = maxGapsBetweenGuideAndPam,
            ))
        }

        val b = builder.copy(guideId=task.id, guide=guide)
        alns.map(aln => b.build(aln))
      }.seq

      ReferenceHit.sort(results, Left(refFile.getSequenceDictionary)).foreach { hit =>
        out.write(hit)
      }
    }

    out.close()
  }
}
