package com.editasmedicine.aligner

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.{Alignment, Cigar, CigarElem, Mode, Aligner => FgAligner}
import com.fulcrumgenomics.util.Sequences
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.reference.ReferenceSequenceFile
import htsjdk.samtools.util.SequenceUtil

import scala.collection.mutable
import scala.math.{max, min, abs}

object SequentialGuideAligner {
  /** Character used to represent a gap in the alignment string. */
  val GapChar: Char = '~'

  object Defaults {
    val MismatchNetCost: Int    = -120
    val GuideGapNetCost: Int = -121
    val GenomeGapNetCost: Int  = -122
    val PamMismatchNetCost: Int = -260

    val MaxGuideDiffs             = 5
    val MaxPamMismatches          = 1
    val MaxGapsBetweenGuideAndPam = 3
    val MaxOverlap                = 10
    val MaxVariantsInCluster      = 16
  }

  /** A case class that represents a guide plus optional PAM and provides easy access to all parts of
    * the sequence in as-is and reverse complemented. */
  case class Guide private (guide: String, pams3Prime: Seq[String], pams5Prime: Seq[String]) {
    require(pams3Prime.isEmpty || pams5Prime.isEmpty, "Guide cannot have both 3' and 5' PAMs.")

    val pamIsFivePrime: Boolean  = pams5Prime.nonEmpty
    val pamIsThreePrime: Boolean = pams3Prime.nonEmpty
    private val pams             = if (pamIsFivePrime) pams5Prime.toArray else pams3Prime.toArray

    private[SequentialGuideAligner] val guideFwBytes: Array[Byte] = guide.getBytes
    private[SequentialGuideAligner] val guideRcBytes: Array[Byte] = Sequences.revcomp(guide).getBytes
    private[SequentialGuideAligner] val pamFwBytes  : Array[Array[Byte]] = pams.map(_.getBytes())
    private[SequentialGuideAligner] val pamRcBytes  : Array[Array[Byte]] = pams.map(Sequences.revcomp(_).getBytes)

    /** The length of the protospacer for this guide. */
    val protospacerLength: Int = this.guideFwBytes.length

    /** The maximum length of the PAM for this guide. */
    val pamLength: Int = if (pams.nonEmpty) pams.map(_.length).max else 0

    /** The maximum length of the whole guide+longest PAM. */
    val length: Int = protospacerLength + pamLength
  }

  object Guide {
    /**
      * Constructs a new guide while forcing all sequences to be the appropriate case.
      *
      * @param guide the guide without any PAM sequences
      * @param pams3Prime zero or more PAMs at the three prime end of the guide
      * @param pams5Prime zero or more PAMs at the five prime end of the guide
      */
    def apply(guide: String, pams3Prime: Seq[String], pams5Prime: Seq[String]): Guide = {
      new Guide(
        guide      = guide.toUpperCase,
        pams3Prime = pams3Prime.map(_.toLowerCase),
        pams5Prime = pams5Prime.map(_.toLowerCase)
      )
    }


    /**
      * Constructs a guide from a single sequence that has the protospacer sequence in upper case and optionally
      * a single PAM sequence in lower case at either the 5' or 3' end of the sequence.
      */
    def apply(sequence: String): Guide = apply(sequence, None)

    /**
      * Constructs a guide from a single sequence that has the protospacer sequence in upper case and optionally
      * a single PAM sequence in lower case at either the 5' or 3' end of the sequence.
      */
    def apply(sequence: String, auxPams: Iterable[String]): Guide = {
      val (guide: String, pam: Option[String], pamIsFivePrime: Boolean, pamIsThreePrime: Boolean) = {
        val parts = splitByCase(sequence.trim)
        require(parts.length <= 2, s"Invalid Guide sequence $sequence.")
        require(parts.length == 2 || parts.head.charAt(0).isUpper, "Guide sequence cannot be all lower case.")
        require(auxPams.isEmpty || parts.length == 2, "Cannot provide auxiliary PAMs without providing a PAM in the guide sequence.")
        require(auxPams.forall(p => p == p.toLowerCase), s"All PAMs must be lower case. PAMs given: ${auxPams.mkString(", ")}")

        if (parts.length == 1) {
          (parts.head, None, false, false)
        }
        else if (parts.head.charAt(0).isUpper) {
          (parts(0), Some(parts(1)), false, true)
        }
        else {
          (parts(1), Some(parts(0)) , true, false)
        }
      }

      if (auxPams.nonEmpty && pam.isEmpty) throw new IllegalStateException("Cannot give aux PAMs without a primary PAM.")
      val pams = pam.toSeq ++ auxPams

      apply(guide      = guide,
            pams3Prime = if (pamIsThreePrime) pams else Nil,
            pams5Prime = if (pamIsFivePrime)  pams else Nil
      )
    }

    /** Splits a String into a series of contiguous chunks of the same case (i.e. upper or lower case). */
    private def splitByCase(bases: String): Seq[String] = {
      val iter = bases.iterator.bufferBetter
      val builder = IndexedSeq.newBuilder[String]

      while (iter.hasNext) {
        val first = iter.head.isLower
        val chunk = iter.takeWhile(_.isLower == first).toArray
        builder += new String(chunk)
      }

      builder.result()
    }
  }


  /** Alignment scorer that charges GAPs differently on the query vs. target to make the net cost
    * of mismatches and gap-bases equivalent.
    */
  final class GuideAlignmentScorer(val matchScore: Int,
                                   val mismatchScore: Int,
                                   val pamMatchScore: Int,
                                   val pamMismatchScore: Int,
                                   val queryGapScore: Int,
                                   val targetGapScore: Int
                                  ) extends FgAligner.AlignmentScorer {
    private val N: Byte = 'N'.toByte
    private val n: Byte = 'n'.toByte

    /** Scores a base pairing as a match if the two bases are compatible, case-insensitive. */
    override def scorePairing(query: Byte, target: Byte): Int = {
      val isPam = Character.isLowerCase(query)
      val m     = if (isPam) pamMatchScore    else matchScore
      val mm    = if (isPam) pamMismatchScore else mismatchScore

      if (target == N || target == n) mm
      else if (Sequences.compatible(query, target)) m
      else mm
    }

    /** Scores a gap solely based on whether it appears on the query or target. */
    override def scoreGap(query: Array[Byte], target: Array[Byte], qOffset: Int, tOffset: Int, gapIsInQuery: Boolean, extend: Boolean): Int = {
      if (gapIsInQuery) queryGapScore
      else targetGapScore
    }
  }
}

/**
  * A class for performing a sequential alignment of a guides plus optional PAMs to other sequences.
  * Alignments are performed by finding all acceptable alignments of the guide to the target sequence,
  * followed by an extension phase in which the PAM is aligned adjacent to the guide with up to some
  * number of gap bases in between.
  *
  * Great care is taken, both in the scoring system and the way in which the sequences are aligned, to
  * ensure that if there is a valid alignment within the given criteria that one will be reported. Note
  * that this is not the same as guaranteeing all possible alignments are returned - e.g. if two alignments
  * start and end in the same position but have gaps in different places, only one will be returned.
  *
  * @param refFile an optional reference sequence file if alignment queries vs. a genome are desired
  */
class SequentialGuideAligner(val refFile: Option[ReferenceSequenceFile] = None,
                             val mismatchNetCost: Int    = SequentialGuideAligner.Defaults.MismatchNetCost,
                             val genomeGapNetCost: Int  = SequentialGuideAligner.Defaults.GenomeGapNetCost,
                             val guideGapNetCost: Int = SequentialGuideAligner.Defaults.GuideGapNetCost,
                             val pamMismatchNetCost: Int = SequentialGuideAligner.Defaults.PamMismatchNetCost
                            ) {
  import SequentialGuideAligner._

  refFile.foreach(r => require(r.isIndexed, s"Cannot work with a non-indexed reference: $refFile"))

  /* Convert the net costs into traditional aligner score parameters.  Allows that each may be positive or
     negative, and uses `abs` and `-` to ensure correct sign (positive for match, negative for everything else).

     Bulges are treated as following:
       position     : 1234567890
       guide/query  : ACGTTA-CGT
       alignment    : ||||~|~|||
       genome/target: ACGT-AACGT

     The first gap, at position 5, inserts a gap in the genome/target due to a bulge (extra T base) in the query/guide.
     The second gap, at position 7, inserts a gap in the guide/query due to a bulge (extra A base) in the genome/target.
   */
  val scorer : GuideAlignmentScorer = {
    val matchScore       = abs(mismatchNetCost) / 2
    val mismatchScore    = - (abs(mismatchNetCost) - matchScore)
    val queryGapScore    = - abs(guideGapNetCost)
    val targetGapScore   = - abs(genomeGapNetCost) + matchScore
    val pamMatchScore    = abs(pamMismatchNetCost) / 2
    val pamMismatchScore = - (abs(pamMismatchNetCost) - pamMatchScore)

    new GuideAlignmentScorer(
      matchScore       = matchScore,
      mismatchScore    = mismatchScore,
      pamMatchScore    = pamMatchScore,
      pamMismatchScore = pamMismatchScore,
      queryGapScore    = queryGapScore,
      targetGapScore   = targetGapScore
    )
  }

  private val fgAligner = new FgAligner(scorer, useEqualsAndX=true, mode=Mode.Glocal)

  /** What's the worst net cost of introducing a difference in the guide region. */
  private val worstGuideDiffScore: Int = Seq(mismatchNetCost, genomeGapNetCost, guideGapNetCost).map(s => -abs(s)).min

  /**
    * TBD
    *
    * @param guide the guide to be aligned
    * @param target the target sequence to be aligned to
    * @param targetName the optional name of the target sequence being aligned to
    * @param targetOffset the optional offset within the target sequence at which the provided target bases occur
    * @param maxGuideDiffs the maximum number of mismatches and inserted/deleted bases in the guide region
    * @param maxPamDiffs the maximum number of mismatches and inserted/deleted bases in the PAM region
    * @param maxOverlap the maximum overlap to be tolerated between alignments that are returned
    *
    * @return a [[GuideAlignment]].  Note that the guide, chrom and strand fields have no useful values.
    */
  def align(guide: Guide,
            target: Array[Byte],
            targetName: String = "n/a",
            targetOffset: Int = 0,
            maxGuideDiffs: Int,
            maxGapsBetweenGuideAndPam: Int,
            maxPamDiffs: Int,
            maxTotalDiffs: Int,
            maxOverlap: Int = 0): Seq[GuideAlignment] = {

    // The minimum score of any alignment we'd like back from the underlying aligner
    val minGuideScore = {
      val matches = scorer.matchScore * guide.guideFwBytes.length
      val diffs   = worstGuideDiffScore * maxGuideDiffs
      matches + diffs
    }

    // Note that during `extendAndFilterRight()` we could use `maxTotalDiffs` as the limit for total differences
    // in the alignment, however this will cause some differences in the results (e.g. it may return a lower scoring
    // alignment that has more total diffs).  As such we use the sum of the three diff totals and apply the
    // `maxTotalDiffs` later as a post-filtering cutoff.
    val maxDiffsDuringFiltering = maxGuideDiffs + maxGapsBetweenGuideAndPam + maxPamDiffs

    // Make a copy of the target sequence that is reverse complemented
    val rcTarget = target.clone()
    Sequences.revcomp(rcTarget)

    // We always want to perform the alignments such that the spot where the PAM is to be aligned
    // subsequently is at the right-hand end of the query sequence. This is due to the fact that we
    // can be guaranteed to find a matching alignment, if one exists, for any given alignment _end_,
    // but the same guarantee is not made for the _start_.  And we want to ensure we examine all possible
    // PAM sites that have a valid alignment terminating at/adjacent to them.
    val (fwdGuideAlns, revGuideAlns) = if (guide.pamIsFivePrime) {
      val fs          = fgAligner.align(query=guide.guideRcBytes, target=rcTarget, minScore=minGuideScore)
      val filteredFs  = extendAndFilterRight(fs, guide.pamRcBytes, rcTarget, maxGuideDiffs, maxPamDiffs, maxGapsBetweenGuideAndPam, maxDiffsDuringFiltering)
      val correctedFs = filteredFs.map { a =>
        val ga = toGuideAlignment(guide, a, targetName, 0, '+')
        ga.copy(
          guide            = rc(ga.guide),
          cigar            = ga.cigar.reverse,
          paddedGuide      = rc(ga.paddedGuide),
          paddedAlignment  = ga.paddedAlignment.reverse,
          paddedTarget     = rc(ga.paddedTarget),
          startOffset      = targetOffset + target.length - ga.endOffset,
          endOffset        = targetOffset + target.length - ga.startOffset,
          guideStartOffset = targetOffset + target.length - ga.guideEndOffset,
          guideEndOffset   = targetOffset + target.length - ga.guideStartOffset
        )
      }

      val rs          = fgAligner.align(query=guide.guideRcBytes, target=target, minScore=minGuideScore)
      val filteredRs  = extendAndFilterRight(rs, guide.pamRcBytes, target, maxGuideDiffs, maxPamDiffs, maxGapsBetweenGuideAndPam, maxDiffsDuringFiltering)
      val correctedRs = filteredRs.map { a =>
        val ga = toGuideAlignment(guide, a, targetName, targetOffset, '+')
        ga.copy(
          guide           = rc(ga.guide),
          cigar           = ga.cigar.reverse,
          strand          = '-',
          paddedGuide     = rc(ga.paddedGuide),
          paddedAlignment = ga.paddedAlignment.reverse,
          paddedTarget    = rc(ga.paddedTarget)
        )
      }

      (correctedFs, correctedRs)
    }
    else {
      val fs          = fgAligner.align(query=guide.guideFwBytes, target=target, minScore=minGuideScore)
      val filteredFs  = extendAndFilterRight(fs, guide.pamFwBytes, target, maxGuideDiffs, maxPamDiffs, maxGapsBetweenGuideAndPam, maxDiffsDuringFiltering)
      val correctedFs = filteredFs.map(a => toGuideAlignment(guide, a, targetName, targetOffset, targetStrand='+'))

      val rs          = fgAligner.align(query=guide.guideFwBytes, target=rcTarget, minScore=minGuideScore)
      val filteredRs  = extendAndFilterRight(rs, guide.pamFwBytes, rcTarget, maxGuideDiffs, maxPamDiffs, maxGapsBetweenGuideAndPam, maxDiffsDuringFiltering)
      val correctedRs = filteredRs.map { a =>
        val ga = toGuideAlignment(guide, a, targetName, 0, targetStrand='+')
        ga.copy(
          strand           = '-',
          startOffset      = targetOffset + target.length - ga.endOffset,
          guideStartOffset = targetOffset + target.length - ga.guideEndOffset,
          endOffset        = targetOffset + target.length - ga.startOffset,
          guideEndOffset   = targetOffset + target.length - ga.guideStartOffset,
        )
      }

      (correctedFs, correctedRs)
    }

    val retval = new mutable.ArrayBuffer[GuideAlignment]()
    for (alns <- Seq(fwdGuideAlns, revGuideAlns); aln <- alns.sorted) {
      if (aln.edits <= maxTotalDiffs && !retval.exists(k => k.strand == aln.strand && k.overlap(aln) > maxOverlap)) {
        retval += aln
      }
    }

    retval.toSeq
  }


  /**
    * Aligns two sequences against one another in Glocal mode.
    *
    * @param guide the guide to be aligned
    * @param target the target sequence to be aligned to
    * @return a [[GuideAlignment]].  Note that the guide, chrom and strand fields have no useful values.
    */
  def alignBest(guide: Guide,
                target: Array[Byte],
                maxGapsBetweenGuideAndPam: Int = Defaults.MaxGapsBetweenGuideAndPam): GuideAlignment = {
    val alns = align(
      guide                     = guide,
      target                    = target,
      maxGuideDiffs             = guide.protospacerLength,
      maxGapsBetweenGuideAndPam = maxGapsBetweenGuideAndPam,
      maxPamDiffs               = guide.pamLength,
      maxTotalDiffs             = guide.protospacerLength + maxGapsBetweenGuideAndPam + guide.pamLength
    )
    alns.maxBy(_.score)
  }


  /**
    * Aligns a guide sequence to a region around the given position.  Attempts alignment on both the F and R
    * strands and returns the best scoring alignment of the pair.
    *
    * @param guide the guide sequence; must be only A/C/G/T and N.
    * @param chrom the chromosome on which the putative hit is found
    * @param pos the approximate location at which to align the guide
    * @param windowSize the number of bases centered around `pos` to look at when aligning. If not provided the
    *                   default is to use `guide.length * 4`.
    * @return the best alignment of the guide to either the F or R strand
    */
  def alignToRef(guide: Guide,
                 chrom: String,
                 pos: Int,
                 windowSize: Option[Int] = None,
                 maxGuideDiffs: Int,
                 maxGapsBetweenGuideAndPam: Int,
                 maxPamDiffs: Int,
                 maxTotalDiffs: Int,
                 maxOverlap: Int = 0): Seq[GuideAlignment] = {
    require(refFile.isDefined, "Cannot perform alignments to ref without a ref fasta!")
    val ref = refFile.get.getSequenceDictionary.getSequence(chrom)
    require(ref != null, s"Unknown chromosome: $chrom")

    val padding = windowSize.map(_ / 2).getOrElse(guide.length * 2)
    val (regionStart, regionEnd) = (max(pos-padding, 1), min(pos+padding, ref.getSequenceLength))
    val target = refFile.get.getSubsequenceAt(chrom, regionStart, regionEnd).getBases

    align(
      guide                     = guide,
      target                    = target,
      targetName                = chrom,
      targetOffset              = regionStart-1,
      maxGuideDiffs             = maxGuideDiffs,
      maxGapsBetweenGuideAndPam = maxGapsBetweenGuideAndPam ,
      maxPamDiffs               = maxPamDiffs,
      maxTotalDiffs             = maxTotalDiffs,
      maxOverlap                = maxOverlap
    ).sorted
  }


  /**
    * Aligns a guide sequence to a region around the given position and returns the best alignment.
    * Attempts alignment on both the F and R strands and returns the best scoring alignment of the pair.
    *
    * @param guide the guide sequence; must be only A/C/G/T and N.
    * @param chrom the chromosome on which the putative hit is found
    * @param pos the approximate location at which to align the guide
    * @param windowSize the number of bases centered around `pos` to look at when aligning. If not provided the
    *                   default is to use `guide.length * 4`.
    * @param maxGapsBetweenGuideAndPam the maximum gap size allowed between the guide and PAM
    * @return the best alignment of the guide to either the F or R strand
    */
  def alignToRefBest(guide: Guide,
                     chrom: String,
                     pos: Int,
                     windowSize: Option[Int] = None,
                     maxGapsBetweenGuideAndPam: Int = Defaults.MaxGapsBetweenGuideAndPam): GuideAlignment = {
    alignToRef(
      guide                     = guide,
      chrom                     = chrom,
      pos                       = pos,
      windowSize                = windowSize,
      maxGuideDiffs             = guide.protospacerLength,
      maxGapsBetweenGuideAndPam = maxGapsBetweenGuideAndPam,
      maxPamDiffs               = guide.pamLength,
      maxTotalDiffs             = guide.protospacerLength + maxGapsBetweenGuideAndPam + guide.pamLength,
      maxOverlap=0
    ).head
  }


  /**
    * Extends alignments to the right with the PAM.  Filters out any alignments that, after extension,
    * have too many mismatches in the PAM.  Will return at most one output alignment per input alignment,
    * prefering alignments with fewer gaps between the guide and PAM.
    *
    * @param alns the set of alignments to extend and filter
    * @param pams the PAM to use to extend the alignment
    * @param target the target sequence that the alignments use
    * @param maxPamMismatches the maximum number of mismatches allowed in the PAM alignment
    * @param maxGapBeforeExtending the maximum gap bases allowed before aligning the PAM
    * @return the filtered set of alignments who meet the given thresholds after extension
    */
  private def extendAndFilterRight(alns: Seq[Alignment],
                                   pams: Array[Array[Byte]],
                                   target: Array[Byte],
                                   maxGuideDiffs: Int,
                                   maxPamMismatches: Int,
                                   maxGapBeforeExtending: Int,
                                   maxTotalDiffs: Int): Seq[Alignment] = {
    // Calculate the number of diffs per guide
    val alnsWithDiffs = alns.iterator.map { aln =>
      (aln, aln.cigar.iterator.filter(_.operator != CigarOperator.EQ).sumBy(_.length))
    }

    // If there are no PAMs just short circuit
    if (pams.isEmpty || pams.length == 1 && pams(0).length == 0) {
      alnsWithDiffs.filter(_._2 <= maxGuideDiffs).map(_._1).toIndexedSeq
    }
    else {
      alnsWithDiffs.filter(_._2 <= maxGuideDiffs).flatMap { case (aln, guideDiffs) =>
        // If the alignment already ends with a gap, we have to reduce the max gap we're allowed to add here
        val terminalGap = if (aln.cigar.last.operator.isIndel) aln.cigar.last.length else 0
        val maxExtraGap = min(maxGapBeforeExtending - terminalGap, maxTotalDiffs - guideDiffs)

        pams.flatMap { pam =>
          val pamLen = pam.length
          val extended = Range.inclusive(0, maxExtraGap).flatMap { offset =>
            val tOffset = aln.targetEnd + offset // targetEnd is 1-based so gives us the 0-based offset of the base after the alignment
            val pamMismatchLimit = min(maxPamMismatches, maxTotalDiffs - guideDiffs - offset)

            if (tOffset + pamLen > target.length || pamMismatchLimit < 0) None else {
              val ops = new Array[Char](pamLen)
              var score = 0

              forloop(from = 0, until = pamLen) { i =>
                val addend = scorer.scorePairing(pam(i), target(tOffset + i))
                score += addend
                ops(i) = if (addend > 0) '=' else 'X'
              }

              if (ops.count(_ == 'X') > pamMismatchLimit) None else {
                val newCigar = Cigar(
                  aln.cigar.elems ++
                    (if (offset > 0) Seq(CigarElem(CigarOperator.D, offset)) else Nil) ++
                    ops.map(o => CigarElem(CigarOperator.characterToEnum(o), 1))
                ).coalesce

                Some(aln.copy(
                  query = concatenate(aln.query, pam),
                  queryStart = 1,
                  cigar = newCigar,
                  score = aln.score + score + (offset * scorer.queryGapScore)
                ))
              }
            }
          }

          if (extended.isEmpty) None else Some(extended.maxBy(_.score))
        }
      }.toIndexedSeq
    }
  }

  /** Concatenates two arrays together. */
  private def concatenate(a1: Array[Byte], a2: Array[Byte]): Array[Byte] = {
    val out = new Array[Byte](a1.length + a2.length)
    System.arraycopy(a1, 0, out, 0, a1.length)
    System.arraycopy(a2, 0, out, a1.length, a2.length)
    out
  }

  /** Converts an Alignment into a GuideAlignment.
    * Alignment uses 1-based closed-ended coordiantes and GuideAlignment uses 0-based open-ended.
    *  */
  private def toGuideAlignment(guide: Guide,
                               alignment: Alignment,
                               targetName: String = "n/a",
                               targetOffset: Int = 0,
                               targetStrand: Char = '.',
                              ): GuideAlignment = {
      val Seq(paddedGuide, alignString, paddedTarget) = alignment.paddedString(gapChar=SequentialGuideAligner.GapChar)
      GuideAlignment(
        guide           = new String(alignment.query),
        chrom           = targetName,
        startOffset     = targetOffset + alignment.targetStart - 1,
        endOffset       = targetOffset + alignment.targetEnd,
        strand          = targetStrand,
        score           = alignment.score,
        cigar           = alignment.cigar,
        paddedGuide     = paddedGuide,
        paddedAlignment = alignString,
        paddedTarget    = paddedTarget
      )
  }

  /** Reverse complements a string that has padding in it. */
  private def rc(s: String): String = {
    val bs = s.getBytes
    SequenceUtil.reverse(bs, 0, bs.length)
    forloop (from=0, until=bs.length) { i =>
      val b = bs(i)
      if (b != '-') bs(i) = Sequences.complement(b)
    }

    new String(bs)
  }
}
