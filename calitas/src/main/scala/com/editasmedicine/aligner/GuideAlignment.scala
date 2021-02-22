package com.editasmedicine.aligner

import java.io.PrintStream

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.{Aligner, Cigar}

object GuideAlignment {
  /** Constructs a GuideAlignment while calculating the guide-only coordinates. */
  def apply(guide: String,
            chrom: String,
            startOffset: Int,
            endOffset: Int,
            strand: Char,
            score: Int,
            cigar: Cigar,
            paddedGuide: String,
            paddedAlignment: String,
            paddedTarget: String): GuideAlignment = {

    val (guideOnlyStart, guideOnlyEnd) = {
      val paddedStart = paddedGuide.indexWhere(_.isUpper)
      val paddedEnd   = paddedGuide.lastIndexWhere(_.isUpper)
      val leftDelta   = Range(0, paddedStart).count(i => paddedTarget(i).isLetter)
      val rightDelta  = Range(paddedEnd+1, paddedTarget.length).count(i => paddedTarget(i).isLetter)

      strand match {
        case '+' => (startOffset + leftDelta, endOffset - rightDelta)
        case '-' => (startOffset + rightDelta, endOffset - leftDelta)
      }
    }

    require(guideOnlyStart >= startOffset)
    require(guideOnlyEnd   <= endOffset)

    apply(
      guide            = guide,
      chrom            = chrom,
      startOffset      = startOffset,
      endOffset        = endOffset,
      strand           = strand,
      guideStartOffset = guideOnlyStart,
      guideEndOffset   = guideOnlyEnd,
      score            = score,
      cigar            = cigar,
      paddedGuide      = paddedGuide,
      paddedAlignment  = paddedAlignment,
      paddedTarget     = paddedTarget
    )
  }
}

/** Represents the alignment of a guide to a section of a reference sequence / chromosome. The alignment
  * is always returned in the same orientation as the guide.  I.e. if the guide matched the negative
  * strand of the genome, the target/genome sequence will be reverse complemented, and the guide returned
  * as is.
  *
  * @param guide the guide sequence (including the PAM if given)
  * @param chrom the chromosome to which the guide was aligned
  * @param startOffset the 0-based offset of the first base on the chromosome to which the guide+PAM is aligned (always < end)
  * @param endOffset the 0-based offset of the chromosomal base after the last aligned base of the guide+PAM
  * @param guideStartOffset the 0-based offset of the first base on the chromosome to which the guide (w/o PAM) is aligned
  * @param guideEndOffset the 0-based offset of the chromosomal base after the last aligned base of the guide (w/o PAM)
  * @param strand the strand of the chromosome to which the guide is aligned
  * @param score the score of the alignment
  * @param cigar the cigar for the alignment
  * @param paddedGuide the padded guide sequence
  * @param paddedAlignment a string the same length as the paddedGuide and paddedTarget which represents
  *                        whether each base is a match/mismatch or part of a gap
  * @param paddedTarget the padded target sequence
  */
case class GuideAlignment(guide: String,
                          chrom: String,
                          startOffset: Int,
                          endOffset: Int,
                          guideStartOffset: Int,
                          guideEndOffset: Int,
                          strand: Char,
                          score: Int,
                          cigar: Cigar,
                          paddedGuide: String,
                          paddedAlignment: String,
                          paddedTarget: String,
                          leftOfGuide10bp: Option[String] = None,
                          rightOfGuide10bp: Option[String] = None,
                          leftOfFullAln8bp: Option[String] = None,
                          rightOfFullAln8bp: Option[String] = None
                         ) extends Ordered[GuideAlignment] {
  require(paddedGuide.length  == paddedAlignment.length, "Padded guide and alignment string are different lengths.")
  require(paddedTarget.length == paddedAlignment.length, "Padded target and alignment string are different lengths.")
  require(strand == '+' || strand == '-' || strand == '.', "Strand must be one of [+-.].")

  def isPositiveStrand: Boolean = strand == '+' || strand == '.'
  def isNegativeStrand: Boolean = !isPositiveStrand

  /* For debugging only - print the alignment to stdout. */
  def print(out: PrintStream = System.out): Unit = Seq(paddedGuide, paddedAlignment, paddedTarget).foreach(out.println)

  def mismatches: Int = paddedAlignment.count(_ == '.')
  def gapBases:   Int = paddedAlignment.count(_ == SequentialGuideAligner.GapChar)
  def edits:      Int = paddedAlignment.count(ch => ch == '.' || ch == SequentialGuideAligner.GapChar)

  def guideMismatches: Int   = count(lower=false, bothSides=false, mms=true,  gaps=false)
  def guideGapBases  : Int   = count(lower=false, bothSides=false, mms=false, gaps=true )
  def guideMmsPlusGaps : Int = count(lower=false, bothSides=false, mms=true,  gaps=true )
  def pamMismatches:   Int   = count(lower=true,  bothSides=true,  mms=true,  gaps=false)
  def pamGapBases:     Int   = count(lower=true,  bothSides=true,  mms=false, gaps=true )
  def pamMmsPlusGaps:  Int   = count(lower=true,  bothSides=true,  mms=true,  gaps=true )

  /** Retrieves the target sequence excluding the PAM portion of the alignment, and without padding. */
  def unpaddedTargetWithoutPam: String = {
    val paddedStart = paddedGuide.indexWhere(_.isUpper)
    val paddedEnd   = paddedGuide.lastIndexWhere(_.isUpper)
    paddedTarget.substring(paddedStart, paddedEnd+1).filter(_.isLetter)
  }


  /** Computes the number of bases overlapping between two guide alignments. */
  def overlap(that: GuideAlignment): Int = if (this.chrom != that.chrom) 0 else {
    val o = math.min(this.endOffset, that.endOffset) - math.max(this.startOffset, that.startOffset)
    if (o > 0) o else 0
  }

  /** Sorts alignments by highest score and the fewest gaps if scores are the same. */
  override def compare(that: GuideAlignment): Int = {
    var retval = that.score - this.score                    // Max by score first
    if (retval == 0) retval = this.gapBases - that.gapBases // Then tie break on gaps
    retval
  }

  /** Counts up the number of gaps and or mismatches in regions of lower or upper case guide sequence.
    *
    * @param lower whether to count where the guide is lower case or upper case
    * @param bothSides whether to require a case match on one or both sides of a gap
    * @param mms true to count mismatches, false to ignore mismatches
    * @param gaps true to count gap bases, false to ignore gaps
    * @return the count of relevant bases
    */
  private def count(lower: Boolean, bothSides: Boolean, mms: Boolean, gaps: Boolean): Int = {
    var n = 0
    val len = paddedAlignment.length

    forloop (from=0, until=len) { i =>
      if (mms && paddedAlignment(i) == '.' && paddedGuide.charAt(i).isLower == lower) n += 1
      else if (gaps && paddedAlignment(i) == SequentialGuideAligner.GapChar) {
        val guideBase = paddedGuide.charAt(i)
        val countMe = (guideBase != '-' && guideBase.isLower == lower) || {
          val prev = previousNonDash(i, paddedGuide)
          val next = nextNonDash(i, paddedGuide)

          if (bothSides) {
            (prev == '-' || prev.isLower == lower) && (next == '-' || next.isLower == lower)
          } else {
            (prev.isLetter && prev.isLower == lower) || (next.isLetter && next.isLower == lower)
          }
        }

        if (countMe) n += 1
      }
    }

    n
  }

  /** Returns the first non-hyphen character from an index prior to `from` in the string.
    * If there is no non-hypen prior, will return a hyphen!!
    */
  def previousNonDash(from: Int, s: String): Char = {
    var i = from
    while (i > 0 && s.charAt(i) == '-') i -= 1
    s.charAt(i)
  }

  /** Returns the first non-hyphen character from an index prior to `from` in the string.
    * If there is no non-hypen prior, will return a hyphen!!
    */
  def nextNonDash(from: Int, s: String): Char = {
    var i = from
    val lastIndex = s.length - 1
    while (i < lastIndex && s.charAt(i) == '-') i += 1
    s.charAt(i)
  }
}
