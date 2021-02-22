package com.editasmedicine.aligner

import com.editasmedicine.commons.testing.UnitSpec
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.Cigar

class GuideAlignmentTest extends UnitSpec {
  /** Filters an alignment string to remove pads. */
  def s(s: String): String = s.filter(_.isLetter)

  "GuideAlignment" should "calculate the number of gaps and mismatches correctly for a perfect alignment" in {
    val paddedQuery  = "GCTGACTGCATGACTATAnrg"
    val paddedAlign  = "|||||||||||||||||||||"
    val paddedTarget = "GCTGACTGCATGACTATAnrg"
    val aln = GuideAlignment(s(paddedQuery), "chr1", 1, 21, '+', 100, Cigar("25M"), paddedQuery, paddedAlign, paddedTarget)

    aln.guideMismatches  shouldBe 0
    aln.guideGapBases    shouldBe 0
    aln.guideMmsPlusGaps shouldBe 0
    aln.pamMismatches    shouldBe 0
    aln.pamGapBases      shouldBe 0
    aln.pamMmsPlusGaps   shouldBe 0
    aln.mismatches       shouldBe 0
    aln.gapBases         shouldBe 0
    aln.edits            shouldBe 0
    aln.guideStartOffset shouldBe 1
    aln.guideEndOffset   shouldBe 18
  }

  it should "correctly count single bp gaps and mismatches in the guide region" in {
    val paddedQuery  = "GCTGACT-GCATGACTATAnrg"
    val paddedAlign  = "||.||||~|||.||~|||||||"
    val paddedTarget = "GCAGACTCGCACGA-TATAnrg"
    val aln = GuideAlignment(s(paddedQuery), "chr1", 1, 21, '+', 100, Cigar("7M1D6M1I7M"), paddedQuery, paddedAlign, paddedTarget)

    aln.guideMismatches  shouldBe 2
    aln.guideGapBases    shouldBe 2
    aln.guideMmsPlusGaps shouldBe 4
    aln.pamMismatches    shouldBe 0
    aln.pamGapBases      shouldBe 0
    aln.pamMmsPlusGaps   shouldBe 0
    aln.mismatches       shouldBe 2
    aln.gapBases         shouldBe 2
    aln.edits            shouldBe 4
    aln.guideStartOffset shouldBe 1
    aln.guideEndOffset   shouldBe 18
  }

  it should "correctly count single bp gaps and mismatches in the PAM region" in {
    val paddedQuery  = "GCTGACTGCATGACTATAnngrrn"
    val paddedAlign  = "|||||||||||||||||||~||.|"
    val paddedTarget = "GCTGACTGCATGACTATAC-GATT"
    val aln = GuideAlignment(s(paddedQuery), "chr1", 1, 23, '+', 100, Cigar("19M1I4M"), paddedQuery, paddedAlign, paddedTarget)

    aln.guideMismatches  shouldBe 0
    aln.guideGapBases    shouldBe 0
    aln.guideMmsPlusGaps shouldBe 0
    aln.pamMismatches    shouldBe 1
    aln.pamGapBases      shouldBe 1
    aln.pamMmsPlusGaps   shouldBe 2
    aln.mismatches       shouldBe 1
    aln.gapBases         shouldBe 1
    aln.edits            shouldBe 2
    aln.guideStartOffset shouldBe 1
    aln.guideEndOffset   shouldBe 18
  }

  it should "handle multi-base gaps correctly" in {
    val paddedQuery  = "GCTGAC---TGCATGACTATAnrg"
    val paddedAlign  = "||||||~~~||||~~|||||||||"
    val paddedTarget = "GCTGACGGGTGCA--ACTATACGG"
    val aln = GuideAlignment(s(paddedQuery), "chr1", 1, 22, '-', 100, Cigar("6M3D4M2I9M"), paddedQuery, paddedAlign, paddedTarget)

    aln.guideMismatches  shouldBe 0
    aln.guideGapBases    shouldBe 5
    aln.guideMmsPlusGaps shouldBe 5
    aln.pamMismatches    shouldBe 0
    aln.pamGapBases      shouldBe 0
    aln.pamMmsPlusGaps   shouldBe 0
    aln.mismatches       shouldBe 0
    aln.gapBases         shouldBe 5
    aln.edits            shouldBe 5
    aln.guideStartOffset shouldBe 4
    aln.guideEndOffset   shouldBe 22
  }

  it should "handle leading and trailing deletions" in { // this is an odd case, and may not come up in real world usage
    val paddedQuery  = "---GCTGACTGCATGACTATAnrg--"
    val paddedAlign  = "~~~|||||||||||||||||||||~~"
    val paddedTarget = "TGTGCTGACTGCATGACTATACGGCC"
    val aln = GuideAlignment(s(paddedQuery), "chr1", 1, 26, '+', 100, Cigar("3D21M2D"), paddedQuery, paddedAlign, paddedTarget)

    aln.guideMismatches  shouldBe 0
    aln.guideGapBases    shouldBe 3
    aln.guideMmsPlusGaps shouldBe 3
    aln.pamMismatches    shouldBe 0
    aln.pamGapBases      shouldBe 2
    aln.pamMmsPlusGaps   shouldBe 2
    aln.mismatches       shouldBe 0
    aln.gapBases         shouldBe 5
    aln.edits            shouldBe 5
    aln.guideStartOffset shouldBe 4
    aln.guideEndOffset   shouldBe 21
  }

  it should "count a gap between the guide and PAM as being in the guide" in {
    val paddedQuery  = "GCTGACTGCATGACTATA--nrg"
    val paddedAlign  = "||||||||||||||||||~~|||"
    val paddedTarget = "GCTGACTGCATGACTATATTCGG"
    val aln = GuideAlignment(s(paddedQuery), "chr1", 1, 23, '+', 100, Cigar("18M2D3M"), paddedQuery, paddedAlign, paddedTarget)

    aln.guideMismatches  shouldBe 0
    aln.guideGapBases    shouldBe 2
    aln.guideMmsPlusGaps shouldBe 2
    aln.pamMismatches    shouldBe 0
    aln.pamGapBases      shouldBe 0
    aln.pamMmsPlusGaps   shouldBe 0
    aln.mismatches       shouldBe 0
    aln.gapBases         shouldBe 2
    aln.edits            shouldBe 2
    aln.guideStartOffset shouldBe 1
    aln.guideEndOffset   shouldBe 18
  }
}
