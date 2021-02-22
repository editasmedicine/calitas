package com.editasmedicine.aligner

import com.editasmedicine.aligner.SequentialGuideAligner.{Guide, GuideAlignmentScorer}
import com.editasmedicine.commons.testing.UnitSpec
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.ReferenceSetBuilder
import com.fulcrumgenomics.util.Sequences
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import htsjdk.samtools.util.SequenceUtil

class SequentialGuideAlignerTest extends UnitSpec {
    private val ref = {
    val builder = new ReferenceSetBuilder()
    builder.add("chr1") // 25 lines x 100 bases each = 2500b
      .add("AAATAGACCTTTCCCATTTATAACTTATTTGTAAAATGATTTCTATTATAAACATAACATATACATTGTATAACAATTAGAAAACCTGTCTGTTTTGATG") // 1-100
      .add("GATCTCAAGATTTAAGAAGGCTTAGACTTCAGCTATAAGATGCACATGCCACTGTGGGAGGCCGAGGCGGGCAGATCACGAGGTCAGGAGTTCTAGACCA") // 101-200
      .add("GCCTGACCAACATGGTGAAACCCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCATGGCAGCAGACACCTGTAATCCCAGTTATTCGGGAGGCTGA") // 201-300
      .add("GGCAGGAGAATTGCTTGAATGCAGGAGGCAGAGGTTGCAGTGAGCCGAGACGGCGCCACTGCACTCCAGCCTGGGCAACAGAGCAGATGGAGACCATCCT") // 301-400
      .add("GACCAACATGATGAAACTCTGTCTCTACTAAAAATACAAAAATTAGCTGGGCATGGTGGCGTGCACCTACTAGTCCCAGCTACTCGGGAGGCTGAGGCAG") // 401-500
      .add("GAGAATTGCTTGAACCCAGGAGGCGGAGGTTTCAGTGAGCCGATACCGCGCCATTGCACTCCAGCCTGGGCAACAGAGCGAGACTGTGTCTCAAAAAAAA") // 501-600
      .add("AAAAAAAAAGGAGATGCACATGTTTAAGTCTATTTCAGGCGGTTAGCTGGTGGATTGCTACAATTCCTCTGTAAGTTTAAAAAATCATGTAAGTGCTGTT") // 601-700
      .add("TTGGAGTACTGTAATAACTCTTGAGATGTAGAACACATCTGCAAAATGAGGGTAGTATAAAAGAGACGAGGGGATGAGGGTAATACATAAGAAATAGGGG") // 701-800
      .add("AAAGGACAAGAACAGGTAAATTAAACTTCAAGTACTATTTTTGCTATTGCTGTCTACACTCAACTAGCAAGGAAAAAGCCTTGCTTCTGCTCTGCGGGTT") // 801-900
      .add("TTCTTCGGGTTTAACTTGACCAAGCAAAACAGACCATCTGGGATTAACTTTTTCCTTTTCACTGTAGGTCACAGGCTCTACGTGTAGGGTGTTGGCCACC") // 901-1000
      .add("TGTTCTTCCACCATCTCTACCTCCACCTCCTCCTTTGTGGCCACAGCAATGTCACAGCCCATACATGGGGGAGGGGAGCATTCAGGAACTCGGAGGCAGA") // 1001-1100
      .add("TGCATTTTTTTCCAAACACAATAACCTCAAACAGTGGTCTCTAAGCACTTTCCTATGCTCTTCCAAAACGTGACCTCCCCTCTTACTCACACATCCCCTA") // 1101-1200
      .add("CACACGGAAAAGGACCACTATCCGTCCAGCCTGCGCTCGAGGGAGAAGTTTATACCTTCGTCCTAGAGATGCCAAATGCAGCAGGGAAGGCTGGACCGAG") // 1201-1300
      .add("GCAGCCGAGTGCTGGAAAGGGAGGCAAGAGGTGCGGGAGCGGGGAGAGGGGGAGGGGAGGCCGGGGCGCCGCGGGAGTAACCTCCACCGCACCCCACCGC") // 1301-1400
      .add("TCCGAGGGGCAGCCGGCCCGGCCCGAGTTTCTCCCCAGAAGCCTCCAGCCGCGGCTCTCGGGGAGGAGGAAGGAAGGGGTTCCCCGTCCAGGAAGCAGCA") // 1401-1500
      .add("CCAGCGGCGACCGCCTCCAGCCTCACCCTCCTCAGCCCCGCACCGCCCATTCCTCACTCCCCGCGCCGCCGCGTCCGCGCGCCTCCCCCCTGCAGACCCC") // 1501-1600
      .add("TCTCACCCAGCCCGCCCCGACCCCGCGCCCGCGCCCCCCACCCGCCCCTCCGGGGACCCCTAATTCATTCACTCGCCGCCGGCCCCGCCCGGCGCCGGCA") // 1601-1700
      .add("AAGAGGGTCGGGACCCGGGCAGGGGCCCAGGAGGGGTGGTCCGCTCCGTACCTCTCTCCCGCACCTGGGAGCCGCTGAGCCTCTGGCCCCGCCGCCGCCT") // 1701-1800
      .add("TCAGTGCCTGCGCCGCGCTCGCTCCCAGTCCGAAATGGCGGGGGCCGGGAGTACTGGCCGAGCCGCCGCCACCTTCGCCGCCGCCACTGCCGCCGCCGCT") // 1801-1900
      .add("GCTGCCTCCGCCGCCGCGGCCGCCGCCTAGGAAAATCGAGCTCCGAGCACACCGATGAGTTCGGGGCCGGGCGGCCGCAGAGGGCAGAGCTATCGATGCG") // 1901-2000
      .add("TTCCGCGCTCGATTCTTCTTCAGACGGGCGTACGAGAGGGAGCGGCTGAGGGCGGTGTGGGAAGAGGGAAGAGGGGGAGGCAGCGAGCGCCGGCGGGGAG") // 2001-2100
      .add("AAGGAGGGGGCCGGGCCGGGCCGGCGGGGGAGGAGCGGGGGCCGGGCCGGCGGAGGAAGGGGTGGCTGGGGCGGTCTAGGGTGGCGAGCCGGGCCGGCTG") // 2101-2200
      .add("GAGAGCGGGTCTGGGCGGCGCCTTGGCGGGAGGAGGGACTGCCGGACCCACGCGGCGGCCCGCCCCCTGCCTAGCCGCAAGGCTGTCCCCGCAGCCGCCA") // 2201-2300
      .add("ATTCTGACCCGGAGCGGGACCGGACCGCGGCGGGCTGTGCGGATGCCACCAGGGAGACGCCGCGAGCGGCCACGCCGCCCCGCTGACCGGTCTCCACAGA") // 2301-2400
    builder.add("chr2")
      .add("GATACaaCTCGTACTGTCAGT")
      .add("GATACGTCTCGTACTGTCAtT")
    val path = builder.toTempFile()
    ReferenceSequenceFileFactory.getReferenceSequenceFile(path)
  }

  val aligner = new SequentialGuideAligner(Some(ref))

  /** Filters an alignment string to remove pads. */
  def s(s: String): String = s.filter(_.isLetter)

  "SequentialGuideAligner.align(query, target)" should "find a perfect alignment of a pam-less guide on the F strand" in {
    val query  =     "AACCAACC"
    val target = "TTTTAACCAACCGGGG"
    val alns   = new SequentialGuideAligner().align(Guide(query), target.getBytes, maxGuideDiffs=0, maxPamDiffs=0, maxGapsBetweenGuideAndPam=0, maxTotalDiffs=0)

    alns should have size 1
    alns.head.strand      shouldBe '+'
    alns.head.startOffset shouldBe 4
    alns.head.endOffset   shouldBe 12
    alns.head.guideStartOffset shouldBe 4
    alns.head.guideEndOffset   shouldBe 12
    alns.head.cigar.toString shouldBe "8="
    alns.head.paddedGuide  shouldBe "AACCAACC"
    alns.head.paddedTarget shouldBe "AACCAACC"
  }

  it should "find a perfect alignment of a pam-less guide on the R strand" in {
    val query  =   "GGTTGGTT"
    val target = "TTAACCAACCGGGG"
    val alns   = new SequentialGuideAligner().align(Guide(query), target.getBytes, maxGuideDiffs=0, maxPamDiffs=0, maxGapsBetweenGuideAndPam=0, maxTotalDiffs=0)

    alns should have size 1
    alns.head.strand      shouldBe '-'
    alns.head.startOffset shouldBe 2
    alns.head.endOffset   shouldBe 10
    alns.head.guideStartOffset shouldBe 2
    alns.head.guideEndOffset   shouldBe 10
    alns.head.cigar.toString shouldBe "8="
    alns.head.paddedGuide  shouldBe "GGTTGGTT"
    alns.head.paddedTarget shouldBe "GGTTGGTT"
  }

  it should "correctly describe an R strand alignment with a mismatch in it" in {
    val query  = "GGTTGGTT"
    val target = "AGCCAACC"
    val alns   = new SequentialGuideAligner().align(Guide(query), target.getBytes, maxGuideDiffs=1, maxPamDiffs=0, maxGapsBetweenGuideAndPam=0, maxTotalDiffs=1)

    alns should have size 1
    alns.head.strand      shouldBe '-'
    alns.head.startOffset shouldBe 0
    alns.head.endOffset   shouldBe 8
    alns.head.guideStartOffset shouldBe 0
    alns.head.guideEndOffset   shouldBe 8
    alns.head.cigar.toString shouldBe "6=1X1="
    alns.head.paddedGuide  shouldBe "GGTTGGTT"
    alns.head.paddedTarget shouldBe "GGTTGGCT"
  }

  it should "extend an alignment to include a PAM on the 3' end on the F strand" in {
    val query =   "AACCAACCAACCnrg"
    val target= "CCAACCAACCAACCGAGGGGGG"
    val alns   = new SequentialGuideAligner().align(Guide(query), target.getBytes, maxGuideDiffs=0, maxPamDiffs=0, maxGapsBetweenGuideAndPam=1, maxTotalDiffs=1)
    alns should have size 1
    alns.head.strand      shouldBe '+'
    alns.head.startOffset shouldBe 2
    alns.head.endOffset   shouldBe 17
    alns.head.guideStartOffset shouldBe 2
    alns.head.guideEndOffset   shouldBe 14
    alns.head.cigar.toString shouldBe "15="
    alns.head.paddedGuide  shouldBe "AACCAACCAACCnrg"
    alns.head.paddedTarget shouldBe "AACCAACCAACCGAG"
  }

  it should "extend an alignment to include a PAM on the 3' end on the R strand" in {
    val query =   "AACCAACCAACCnrg"
    val target= "CCCTGGGTTGGTTGGTTGGGGGG"
    val alns   = new SequentialGuideAligner().align(Guide(query), target.getBytes, maxGuideDiffs=0, maxPamDiffs=0, maxGapsBetweenGuideAndPam=1, maxTotalDiffs=1)
    alns should have size 1
    alns.head.strand      shouldBe '-'
    alns.head.startOffset shouldBe 2
    alns.head.endOffset   shouldBe 17
    alns.head.guideStartOffset shouldBe 5
    alns.head.guideEndOffset   shouldBe 17
    alns.head.cigar.toString shouldBe "15="
    alns.head.paddedGuide  shouldBe "AACCAACCAACCnrg"
    alns.head.paddedTarget shouldBe "AACCAACCAACCCAG"
  }

  it should "extend an alignment to include a PAM on the 5' end on the F strand" in {
    val query =   "tttvAACCAACCAACC"
    val target= "CCTTTGAACCAACCAACCGAGG"
    val alns   = new SequentialGuideAligner().align(Guide(query), target.getBytes, maxGuideDiffs=0, maxPamDiffs=0, maxGapsBetweenGuideAndPam=1, maxTotalDiffs=1)
    alns should have size 1
    alns.head.strand      shouldBe '+'
    alns.head.startOffset shouldBe 2
    alns.head.endOffset   shouldBe 18
    alns.head.guideStartOffset shouldBe 6
    alns.head.guideEndOffset   shouldBe 18
    alns.head.cigar.toString shouldBe "16="
    alns.head.paddedGuide  shouldBe "tttvAACCAACCAACC"
    alns.head.paddedTarget shouldBe "TTTGAACCAACCAACC"
  }

  it should "extend an alignment to include a PAM on the 5' end on the R strand" in {
    val query = "tttvAACCAACCAACC"
    val target= "CC" + Sequences.revcomp(query.replace("tttv", "TTTG")) + "GAGG"
    val alns   = new SequentialGuideAligner().align(Guide(query), target.getBytes, maxGuideDiffs=0, maxPamDiffs=0, maxGapsBetweenGuideAndPam=1, maxTotalDiffs=1)
    alns should have size 1
    alns.head.strand      shouldBe '-'
    alns.head.startOffset shouldBe 2
    alns.head.endOffset   shouldBe 18
    alns.head.guideStartOffset shouldBe 2
    alns.head.guideEndOffset   shouldBe 14
    alns.head.cigar.toString shouldBe "16="
    alns.head.paddedGuide  shouldBe "tttvAACCAACCAACC"
    alns.head.paddedTarget shouldBe "TTTGAACCAACCAACC"
  }

  it should "extend an alignment to include a PAM on the 5' end with a mismatch in the guide on the F strand" in {
    val query =   "tttvAACCAACCAACC"
    val target= "CCTTTGAACCAACCAAGCGAGG"
    val alns   = new SequentialGuideAligner().align(Guide(query), target.getBytes, maxGuideDiffs=1, maxPamDiffs=0, maxGapsBetweenGuideAndPam=0, maxTotalDiffs=1)
    alns should have size 1
    alns.head.strand      shouldBe '+'
    alns.head.startOffset shouldBe 2
    alns.head.endOffset   shouldBe 18
    alns.head.guideStartOffset shouldBe 6
    alns.head.guideEndOffset   shouldBe 18
    alns.head.cigar.toString shouldBe "14=1X1="
    alns.head.paddedGuide  shouldBe "tttvAACCAACCAACC"
    alns.head.paddedTarget shouldBe "TTTGAACCAACCAAGC"
  }

  it should "extend an alignment to include a PAM on the 5' end with a mismatch in the guide on the R strand" in {
    val query = "tttvAACCAACCAACC"
    val target= "CC" + Sequences.revcomp("TTTGAACCAACCAAGC") + "GAGG"
    val alns   = new SequentialGuideAligner().align(Guide(query), target.getBytes, maxGuideDiffs=1, maxPamDiffs=0, maxGapsBetweenGuideAndPam=0, maxTotalDiffs=1)
    alns should have size 1
    alns.head.strand      shouldBe '-'
    alns.head.startOffset shouldBe 2
    alns.head.endOffset   shouldBe 18
    alns.head.guideStartOffset shouldBe 2
    alns.head.guideEndOffset   shouldBe 14
    alns.head.cigar.toString shouldBe "14=1X1="
    alns.head.paddedGuide  shouldBe "tttvAACCAACCAACC"
    alns.head.paddedTarget shouldBe "TTTGAACCAACCAAGC"
  }

  it should "respect the targetOffset when given" in {
    val guide1  = "gggTTTTT"
    val guide2  = "TTTTTggg"
    val target1 = "AGAGAGAGAGGGTTTTTGGGAGAGAGAGAGAGAG"
    val target2 = "AGAGAGAGACCCAAAAACCCAGAGAGAGAGAGAG"

    val aligner = new SequentialGuideAligner()

    val r1 = aligner.align(Guide(guide1), target1.getBytes, maxGuideDiffs=0, maxPamDiffs=0, maxGapsBetweenGuideAndPam=0, maxTotalDiffs=0, targetOffset=1000).head
    r1.startOffset      shouldBe 1009
    r1.endOffset        shouldBe 1017
    r1.guideStartOffset shouldBe 1012
    r1.guideEndOffset   shouldBe 1017

    val r2 = aligner.align(Guide(guide2), target1.getBytes, maxGuideDiffs=0, maxPamDiffs=0, maxGapsBetweenGuideAndPam=0, maxTotalDiffs=0, targetOffset=1000).head
    r2.startOffset      shouldBe 1012
    r2.endOffset        shouldBe 1020
    r1.guideStartOffset shouldBe 1012
    r1.guideEndOffset   shouldBe 1017

    val r3 = aligner.align(Guide(guide1), target2.getBytes, maxGuideDiffs=0, maxPamDiffs=0, maxGapsBetweenGuideAndPam=0, maxTotalDiffs=0, targetOffset=1000).head
    r3.startOffset      shouldBe 1012
    r3.endOffset        shouldBe 1020
    r1.guideStartOffset shouldBe 1012
    r1.guideEndOffset   shouldBe 1017

    val r4 = aligner.align(Guide(guide2), target2.getBytes, maxGuideDiffs=0, maxPamDiffs=0, maxGapsBetweenGuideAndPam=0, maxTotalDiffs=0, targetOffset=1000).head
    r4.startOffset      shouldBe 1009
    r4.endOffset        shouldBe 1017
    r1.guideStartOffset shouldBe 1012
    r1.guideEndOffset   shouldBe 1017
  }

  it should "score alignments equally when reverse complemented" in {
    val query = "AATTCcgg"
    Seq("AATTCCGG", "AGTTCCGG", "AAATTCCGG", "AATTCCGAG", "AATTCCTG").foreach { target =>
      val f = aligner.alignBest(Guide(query), target.getBytes)
      val r = aligner.alignBest(Guide(Sequences.revcomp(query)), Sequences.revcomp(target).getBytes)
      r.score shouldBe f.score
      r.guideMismatches shouldBe f.guideMismatches
      r.guideGapBases   shouldBe f.guideGapBases
      r.pamMismatches   shouldBe f.pamMismatches
      r.pamGapBases     shouldBe f.pamGapBases
    }
  }

  it should "penalize N bases in the reference" in {
    val query   = "AACCGGTTnrg"
    val target  = "nnnnnnnnnnn"
    val result  = aligner.alignBest(Guide(query), target.getBytes)
    result.score shouldBe (8 * aligner.scorer.mismatchScore) + (3 * aligner.scorer.pamMismatchScore)
  }

  it should "find an alignment with the max guide diffs, when some diffs are indels" in {
    val query   = "yttnAGGAAACTTCTGGCAGGACC"
    val target  = "GTTAGTTCCAGATCTTGAGGAAGCTATCCCAGGACCCTGTCGCCACAGCCA"
    val results = aligner.align(Guide(query), target.getBytes, maxGuideDiffs=5, maxGapsBetweenGuideAndPam=1, maxPamDiffs=1, maxTotalDiffs=7, maxOverlap=10)
    results.size shouldBe 1
    results.head.startOffset shouldBe 13
  }

  it should "pick the PAM with the best scoring alignment when given multiple PAMs" in {
    val query = Guide("AACCGGTTACGTnrg", Seq("ntg"))
    val target =      "AACCGGTTACGTTTG"
    val result = aligner.alignBest(query, target.getBytes)
    result.guide shouldBe "AACCGGTTACGTntg"
    result.pamMmsPlusGaps shouldBe 0
  }

  it should "prefer a longer PAM when they match equally well" in {
    val query = Guide("AACCGGTTACGTnnn", Seq("nnnn", "nn"))
    val target =      "AACCGGTTACGTAAAAAAA"
    val result = aligner.alignBest(query, target.getBytes)
    result.guide shouldBe "AACCGGTTACGTnnnn"
  }

  it should "prefer a longer PAM with a gap vs. a shorter PAM without a gap" in {
    val query = Guide("AACCGGTTACGTacc", Seq("cccc"))
    val target =      "AACCGGTTACGTACCCC"
    val result = aligner.alignBest(query, target.getBytes)
    result.guide shouldBe "AACCGGTTACGTcccc"
    result.cigar.toString shouldBe "12=1D4="
  }


  "SequentialGuideAligner.align(chrom, pos)" should "align a perfect sequence without gaps or mismatches in the right spot on the F strand" in {
    val query  = ref.getSubsequenceAt("chr1", 50, 69).getBaseString
    val result = aligner.alignToRefBest(Guide(query), "chr1", 65)

    result.chrom  shouldBe "chr1"
    result.startOffset  shouldBe 49
    result.endOffset    shouldBe 69
    result.strand shouldBe '+'
    result.paddedGuide shouldBe result.paddedTarget
    result.paddedAlignment.forall(_ == '|') shouldBe true
    result.score should be >= 0
  }

  it should "handle Us the same as Ts" in {
    val tQuery  = ref.getSubsequenceAt("chr1", 50, 69).getBaseString
    val uQuery  = tQuery.replace('T', 'U')
    uQuery should not be tQuery

    val tResult = aligner.alignToRefBest(Guide(tQuery), "chr1", 65)
    val uResult = aligner.alignToRefBest(Guide(uQuery), "chr1", 65)
    uResult.score shouldBe tResult.score
    uResult.paddedAlignment shouldBe tResult.paddedAlignment
  }

  it should "align a perfect sequence without gaps or mismatches in the right spot on the R strand" in {
    val query   = SequenceUtil.reverseComplement(ref.getSubsequenceAt("chr1", 50, 69).getBaseString)
    val result  = aligner.alignToRefBest(Guide(query), "chr1", 65)

    result.chrom  shouldBe "chr1"
    result.startOffset  shouldBe 49
    result.endOffset    shouldBe 69
    result.strand shouldBe '-'
    result.paddedAlignment.forall(_ == '|') shouldBe true
    result.score should be >= 0
  }

  it should "align with a mismatch on the F strand" in {
    val query  = "GAGAATTGtTTGAACCCAGGnGG" // start of 5th line == 501-523 (1-based)
    val aligns = "||||||||.||||||||||||||"
    val result = aligner.alignToRefBest(Guide(query.toUpperCase), "chr1", 515)

    result.chrom shouldBe "chr1"
    result.startOffset shouldBe 500
    result.endOffset   shouldBe 523
    result.strand shouldBe '+'
    result.paddedAlignment shouldBe aligns
    result.mismatches shouldBe 1
  }

  it should "handle ambiguity codes correctly in the PAM" in {
    // ref       "TCAGTGCCTGCGCCGCGCTCGCTCCCAGTCCGAAA" // start of 18th line == 1801-1835 (1-based)
    val query  = "TCAGTGCCTGCGCCGCGCTCGCTCCCnrycwshdm"
    val aligns = "||||||||||||||||||||||||||||||.||||"
    val result = aligner.alignToRefBest(Guide(query), "chr1", 1820)

    result.chrom shouldBe "chr1"
    result.startOffset shouldBe 1800
    result.endOffset   shouldBe 1835
    result.guideStartOffset shouldBe 1800
    result.guideEndOffset   shouldBe 1826
    result.strand shouldBe '+'
    result.paddedAlignment shouldBe aligns
    result.mismatches shouldBe 1
  }

  it should "align with a two bulges on the R strand" in {
    val query  = "AGGCTGG-GGCGGTCGCtCGCNGG" // revcomp of start of 16th line == 1501-1523 (1-based)
    val aligns = "|||||||~|||||||||~||||||"
    val result = aligner.alignToRefBest(Guide(query.filter(_.isLetter).toUpperCase), "chr1", 1510)

    result.chrom shouldBe "chr1"
    result.startOffset shouldBe 1500
    result.endOffset   shouldBe 1523
    result.strand shouldBe '-'
    result.paddedAlignment shouldBe aligns
  }

  it should "prefer an alignment with two mismatches in the guide vs. one mismatch in the PAM" in {
    val query = "GATACGTCTCGTACTGTnrg"
    val result = aligner.alignToRefBest(Guide(query), "chr2", 22)
    result.chrom shouldBe "chr2"
    result.startOffset shouldBe 0
    result.endOffset shouldBe 20
    result.gapBases shouldBe 0
    result.mismatches shouldBe 2
  }

  it should "prefer mismatches to a bulges in the genome" in {
    val query =  "GATACGTCTCGTACTGTnrg"
    val target = query.replace("GATA", "GATT").replace("nrg", "AAG") + "TTTTT" + query.replace("TCTC", "TCTCC").replace("nrg", "AAG")
    val result = aligner.alignBest(Guide(query), target.getBytes)
    result.startOffset shouldBe 0
    result.mismatches shouldBe 1
    result.gapBases shouldBe 0
  }

  it should "prefer bulges in the genome to bulges in the guide" in {
    val query =  "GATACGTCTCGTACTGTnrg"
    val target = query.replace("TCTC", "TCTCC").replace("nrg", "AAG") + "NNNNN" + query.replace("TCTC", "TCT").replace("nrg", "AAG")
    val result = aligner.alignBest(Guide(query), target.getBytes)
    result.startOffset shouldBe 0
    result.mismatches shouldBe 0
    result.gapBases shouldBe 1
  }

  it should "enforce the max diffs limit separately from the other maxes" in {
    val query   = "GATACGTCTCGTACTGTnrg"
    val target1 = "GAaACGTtTCGTACTGTaac".toUpperCase.getBytes  // 2 diffs in guide, 1 in PAM
    val guide  = Guide(query)

    val r1 = aligner.align(guide, target1, maxGuideDiffs=2, maxGapsBetweenGuideAndPam=0, maxPamDiffs=1, maxTotalDiffs=3)
    r1.size shouldBe 1

    val r2 = aligner.align(guide, target1, maxGuideDiffs=2, maxGapsBetweenGuideAndPam=0, maxPamDiffs=1, maxTotalDiffs=2)
    r2.size shouldBe 0
  }
}
