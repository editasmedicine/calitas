package com.editasmedicine.aligner

import com.editasmedicine.aligner.SearchReferencesWithVariants.VariantSet
import com.editasmedicine.commons.testing.UnitSpec
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.{ReferenceSetBuilder, VcfBuilder}
import com.fulcrumgenomics.util.Metric
import com.fulcrumgenomics.vcf.api
import com.fulcrumgenomics.vcf.api.Variant
import htsjdk.samtools.reference.ReferenceSequence

class SearchReferenceTest extends UnitSpec {
  private val Guide      = "ACGTACATGCTCGATACGACGnngrrn"
  private val Perfect    = "ACGTACATGCTCGATACGACGccgaat".toUpperCase
  private val Mismatched = "ACGcACAcGCcCGAcACGACGccgaat".toUpperCase

  private val Fasta = {
    val builder = new ReferenceSetBuilder()
    val chr1 = builder.add("chr1")
    chr1.add("N", 5000)
    chr1.add("AATAT", 1000)
    chr1.add("N", 5000)

    val chr2 = builder.add("chr2")
    chr2.add("N", 3000)
    chr2.add(Perfect)
    chr2.add("GT", 500)
    chr2.add(Mismatched)
    chr2.add("CA", 500)
    chr2.add("N", 3000)

    builder.toTempFile()
  }

  /** Constructs a variant for insertion into a sequence. */
  private def v(chrom: String, pos: Int, id: String, alleles: String): Variant = {
    val builder = VcfBuilder(samples=Seq())
    builder.add(chrom, pos, id, alleles.split("/").toSeq)
    builder.iterator.next()
  }


  "SearchReference.windowIterator" should "iterate over a reference" in {
    val builder = new ReferenceSetBuilder()
    builder.add("chr1").add("ACGTC", 5000)
    val ref = builder.toTempFile()
    val iterator = SearchReference.windowIterator(ref, 451, 426, None)
    iterator.count(_ => true)
  }

  "SearchReference" should "run end to end and produce some results" in {
    val out = makeTempFile("hits.", ".txt")
    new SearchReference(guide=Guide, guideId="a", ref=Fasta, output=out, threads=1).execute()
    val hits = Metric.read[ReferenceHit](out)
    hits should have size 2

    hits.foreach(_.chromosome shouldBe "chr2")
    hits(0).coordinate_start shouldBe 3000
    hits(0).total_mm_plus_gaps shouldBe 0
    hits(1).coordinate_start shouldBe 4000 + Perfect.length
    hits(1).total_mm_plus_gaps shouldBe 4
  }

  it should "find hits with a PAM-less guide" in {
    val out = makeTempFile("hits.", ".txt")
    new SearchReference(guide=Guide.filter(_.isUpper), guideId="a", ref=Fasta, output=out, threads=1).execute()
    val hits = Metric.read[ReferenceHit](out)
    hits should have size 2
  }

  it should "find hits on adjacent short contigs" in {
    val refBuilder = new ReferenceSetBuilder()
    "CTGCCCTGACTCCCACCATGC"

    "GTGACTTGAAGTCTCAGTATA"
    refBuilder.add("ref").add("GTGCGTGACTTGAAGTCTCAGTATACCTTGCCACACGTTGCAGGTTGCCC")
    refBuilder.add("alt").add("GTGCGTGACTTGAAGTCTCAGTATgaaaTTGCCACACGTTGCAGGTTGCCC")
    val ref = refBuilder.toTempFile()

    val out = makeTempFile("hits.", ".txt")
    new SearchReference(guide="GTGACTTGAAGTCTCAGTATA", guideId="a", ref=ref, output=out, threads=1).execute()
    val hits = Metric.read[ReferenceHit](out)

    hits should have size 2
    hits(0).chromosome shouldBe "ref"
    hits(0).coordinate_start shouldBe 4
    hits(0).padded_alignment shouldBe "|||||||||||||||||||||"
    
    hits(1).chromosome shouldBe "alt"
    hits(1).coordinate_start shouldBe 4
    hits(1).padded_alignment shouldBe "||||||||||||||||||||."
  }

  it should "pull the appropriate flanking sequences for ref and variant windows" in {
    val query = "GCGTCACGGTCGAGCGATTGnrg"

    // GCGTCACGGTCGAGCGATTGggg
    val refBuilder = new ReferenceSetBuilder()
    val chr1 = refBuilder.add("chr1")

    //        000000000111111111122222222223333333333444444444455555555556666666666777777777788888888889999999999
    //        1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    chr1.add("ACACACACACACACACACACACACACACACACACACACAgcgtcacggtcgagcgattggggAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".toUpperCase) // perfect + strand hit
    chr1.add("ACACACACACACACACACACACACACACACACACACACAccccaatcgctcgaccgtgacgcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".toUpperCase) // perfect - strand hit
    chr1.add("ACACACACACACACACACACACACACACACACACACACAcacggtcgagcgattggggAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".toUpperCase) // + strand starting in ins
    chr1.add("ACACACACACACACACACACACACACACACACACACACAaatcgctcgaccgtgacgcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".toUpperCase) // - strand starting in ins

    val vcfBuilder = VcfBuilder(samples=Seq())
    vcfBuilder.add("chr1", 239, "insGAGGCGT", Seq("A", "AGAGGCGT"))
    vcfBuilder.add("chr1", 339, "insTCGCCCC", Seq("A", "ATCGCCCC"))

    val ref = refBuilder.toTempFile()
    val vcf = vcfBuilder.toTempFile()
    val out = makeTempFile("results.", ".txt")
    new SearchReference(guide=query, guideId="test", ref=ref, variants=Some(vcf), output=out, maxGapsBetweenGuideAndPam=0, maxGuideDiffs=0).execute()

    val hits = Metric.read[ReferenceHit](out)
    hits should have size 4

    // NOTE:
    //  - 10bp flanks are around the protospace and thus will include any sequence the PAM is matched to
    //  -  8bp flanks are around the whole alignment, and thus won't include the PAM sequence
    val Seq(h1, h2, h3, h4) = hits
    h1.coordinate_start shouldBe 39
    h1.padded_extra_8_bases_5_prime shouldBe "CACACACA"
    h1.padded_extra_8_bases_3_prime shouldBe "AAAAAAAA"
    h1.ten_bases_5_prime            shouldBe "CACACACACA"
    h1.ten_bases_3_prime            shouldBe "GGGAAAAAAA"

    h2.coordinate_start shouldBe 142
    h2.padded_extra_8_bases_5_prime shouldBe "TTTTTTTT"
    h2.padded_extra_8_bases_3_prime shouldBe "TGTGTGTG"
    h2.ten_bases_5_prime            shouldBe "TTTTTTTTTT"
    h2.ten_bases_3_prime            shouldBe "GGGTGTGTGT"

    h3.coordinate_start shouldBe 238
    h3.padded_extra_8_bases_5_prime shouldBe "ACACAGAG"
    h3.padded_extra_8_bases_3_prime shouldBe "AAAAAAAA"
    h3.ten_bases_5_prime            shouldBe "ACACACAGAG"
    h3.ten_bases_3_prime            shouldBe "GGGAAAAAAA"

    h4.coordinate_start shouldBe 338
    h4.padded_extra_8_bases_5_prime shouldBe "TTTTTTTT"
    h4.padded_extra_8_bases_3_prime shouldBe "CGATGTGT"
    h4.ten_bases_5_prime            shouldBe "TTTTTTTTTT"
    h4.ten_bases_3_prime            shouldBe "GGGCGATGTG"
  }


  "SearchReference.alleleCombos" should "return an element per allele for a single variant" in {
    SearchReferencesWithVariants.alleleCombos(Seq(2)) shouldBe Array(Array(0), Array(1))
    SearchReferencesWithVariants.alleleCombos(Seq(3)) shouldBe Array(Array(0), Array(1), Array(2))
  }

  it should "return all the pairwise combinations for two variants" in {
    SearchReferencesWithVariants.alleleCombos(Seq(2, 2)) shouldBe Array(Array(0, 0), Array(0, 1), Array(1, 0), Array(1, 1))
    SearchReferencesWithVariants.alleleCombos(Seq(3, 2)) shouldBe Array(Array(0, 0), Array(0, 1), Array(1, 0), Array(1, 1), Array(2, 0), Array(2, 1))
  }

  it should "return all the combinations for three variants" in {
    SearchReferencesWithVariants.alleleCombos(Seq(3, 2, 3)) shouldBe Array(
      Array(0, 0, 0),
      Array(0, 0, 1),
      Array(0, 0, 2),
      Array(0, 1, 0),
      Array(0, 1, 1),
      Array(0, 1, 2),
      Array(1, 0, 0),
      Array(1, 0, 1),
      Array(1, 0, 2),
      Array(1, 1, 0),
      Array(1, 1, 1),
      Array(1, 1, 2),
      Array(2, 0, 0),
      Array(2, 0, 1),
      Array(2, 0, 2),
      Array(2, 1, 0),
      Array(2, 1, 1),
      Array(2, 1, 2),
    )
  }

  "SearchReferenceWithVariants.buildVariantWindow" should "generate a variant window with a single SNP in it" in {
    val window = SearchReferencesWithVariants.buildVariantWindow(
      set     = VariantSet(variants=IndexedSeq(v("chr1", 20, "rs123", "C/G")), alleles=IndexedSeq(1)),
      ref     = new ReferenceSequence("chr1", 0, "CTAGACTGACTGACTAGCACTAGCCGCTTTATATATGCTATGGGACACCG".getBytes),
      padding = 15                            //  12345678901234567890123456789012345678901234567890
    )
                                    // 012345678901234567890123456789012
    new String(window.bases) shouldBe "ACTGACTGACTAGCAgTAGCCGCTTTATATA".toUpperCase
    window.cigar.toString    shouldBe "31M"
    window.refOffsetAtBaseOffset(0,  preceding=true) shouldBe 4
    window.refOffsetAtBaseOffset(15, preceding=true) shouldBe 19
    window.refOffsetAtBaseOffset(20, preceding=true) shouldBe 24
    window.refOffsetAtBaseOffset(31, preceding=true) shouldBe 35
  }

  it should "generate a variant window with a single insertion in it" in {
    val window = SearchReferencesWithVariants.buildVariantWindow(
      set     = VariantSet(variants=IndexedSeq(v("chr1", 20, "rs123", "C/CGT")), alleles=IndexedSeq(1)),
      ref     = new ReferenceSequence("chr1", 0, "CTAGACTGACTGACTAGCACTAGCCGCTTTATATATGCTATGGGACACCG".getBytes),
      padding = 15                            //  12345678901234567890123456789012345678901234567890
    )
                                    // 012345678901234567890123456789012
    new String(window.bases) shouldBe "ACTGACTGACTAGCAcgtTAGCCGCTTTATATA".toUpperCase
    window.cigar.toString    shouldBe "16M2I15M"
    window.refOffsetAtBaseOffset(0,  preceding=true) shouldBe 4
    window.refOffsetAtBaseOffset(14, preceding=true) shouldBe 18
    window.refOffsetAtBaseOffset(15, preceding=true) shouldBe 19
    window.refOffsetAtBaseOffset(16, preceding=true) shouldBe 19
    window.refOffsetAtBaseOffset(17, preceding=true) shouldBe 19
    window.refOffsetAtBaseOffset(15, preceding=false) shouldBe 19
    window.refOffsetAtBaseOffset(16, preceding=false) shouldBe 20
    window.refOffsetAtBaseOffset(17, preceding=false) shouldBe 20
  }

  it should "generate a variant window with a single deletion in it" in {
    val window = SearchReferencesWithVariants.buildVariantWindow(
      set     = VariantSet(variants=IndexedSeq(v("chr1", 20, "rs123", "CTA/C")), alleles=IndexedSeq(1)),
      ref     = new ReferenceSequence("chr1", 0, "CTAGACTGACTGACTAGCACTAGCCGCTTTATATATGCTATGGGACACCG".getBytes),
      padding = 15                            //  12345678901234567890123456789012345678901234567890
    )

                                    // 0123456789012345678901234567890
    new String(window.bases) shouldBe "ACTGACTGACTAGCAcGCCGCTTTATATATG".toUpperCase
    window.cigar.toString    shouldBe "16M2D15M"
    window.refOffsetAtBaseOffset(0,  preceding=true) shouldBe 4
    window.refOffsetAtBaseOffset(15, preceding=true) shouldBe 19
    window.refOffsetAtBaseOffset(16, preceding=true) shouldBe 22
  }

  it should "generate a variant window with a mutiple variants in it" in {
    val variants = IndexedSeq(
      v("chr1", 10, "snp", "C/T"),
      v("chr1", 20, "ins", "C/CG"),
      v("chr1", 30, "del", "TAT/T")
    )

    val window = SearchReferencesWithVariants.buildVariantWindow(
      set     = VariantSet(variants=variants, alleles=IndexedSeq(1, 1, 1)),
      ref     = new ReferenceSequence("chr1", 0, "CTAGACTGACTGACTAGCACTAGCCGCTTTATATATGCTAGGCGCTACTGAATGCTATAGCTCTGAGACTGGGACACCG".getBytes),
      padding = 15                            //  1234567890123456789012345678901234567890123456789012345678901234567890123456789
    )

    new String(window.bases) shouldBe "CTAGACTGAtTGACTAGCAcgTAGCCGCTTtATATGCTAGGCGCTA".toUpperCase
    window.cigar.toString    shouldBe "20M1I10M2D15M"
  }

  "SearchReferencesWithVariants.alleleCombos" should "generate a set of 1 for a single variant" in {
    val vs   = IndexedSeq(v("chr1", 20, "snp", "A/C"))
    val sets = SearchReferencesWithVariants.alleleCombos(vs, 10)

    sets.size shouldBe 1
    sets.head shouldBe VariantSet(vs, IndexedSeq(1))
  }

  it should "generate an entry for each allele of a single variant" in {
    val vs   = IndexedSeq(v("chr1", 20, "snp", "A/C/G/T"))
    val sets = SearchReferencesWithVariants.alleleCombos(vs, 10)

    sets should contain theSameElementsAs Seq(
      VariantSet(vs, IndexedSeq(1)),
      VariantSet(vs, IndexedSeq(2)),
      VariantSet(vs, IndexedSeq(3))
    )
  }

  it should "generate all combinations for multiple variants" in {
    // Note when a variant has the ref allele it is dropped from the variant set
    val a = v("chr1", 20, "a", "A/C")
    val b = v("chr1", 25, "b", "C/T")
    val c = v("chr1", 30, "c", "G/A")
    val sets = SearchReferencesWithVariants.alleleCombos(IndexedSeq(a, b, c), 10)

    sets should contain theSameElementsAs Seq(
      VariantSet(IndexedSeq(a), IndexedSeq(1)),
      VariantSet(IndexedSeq(b), IndexedSeq(1)),
      VariantSet(IndexedSeq(c), IndexedSeq(1)),
      VariantSet(IndexedSeq(a, b), IndexedSeq(1, 1)),
      VariantSet(IndexedSeq(a, c), IndexedSeq(1, 1)),
      VariantSet(IndexedSeq(b, c), IndexedSeq(1, 1)),
      VariantSet(IndexedSeq(a, b, c), IndexedSeq(1, 1, 1)),
    )
  }

  it should "just return a set for the first variant if the number of variants is too high" in {
    val vs = IndexedSeq(
      v("chr1", 20, "a", "A/C"),
      v("chr1", 25, "b", "C/T"),
      v("chr1", 30, "c", "G/A"),
    )

    SearchReferencesWithVariants.alleleCombos(vs, 2) should have size 1
    SearchReferencesWithVariants.alleleCombos(vs, 3) should have size 7
  }
}
