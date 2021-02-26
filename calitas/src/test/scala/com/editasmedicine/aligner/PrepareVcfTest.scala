package com.editasmedicine.aligner

import com.editasmedicine.commons.testing.UnitSpec
import com.fulcrumgenomics.testing.VcfBuilder
import com.fulcrumgenomics.testing.VcfBuilder.{DefaultHeader, Gt}
import com.fulcrumgenomics.vcf.api.{Variant, VcfCount, VcfFieldType, VcfFilterHeader, VcfInfoHeader, VcfSource}

class PrepareVcfTest extends UnitSpec {
  "PrepareVcf" should "correctly strip out genotypes in the input VCF" in {
    val samples = IndexedSeq("sample1", "sample2")
    val header  = VcfBuilder.DefaultHeader.copy(
      infos   = DefaultHeader.infos ++ Seq(VcfInfoHeader("AF", VcfCount.OnePerAltAllele, VcfFieldType.Float, "ALT allele frequency")),
      filters = DefaultHeader.filters ++ Seq(VcfFilterHeader("PASS", "Passes all filters.")),
      samples = samples
    )
    val builder = VcfBuilder(header=header)

    Range(0, 10).foreach { i =>
      builder.add(
        chrom   = "chr1",
        pos     = 1000 * (i+1),
        alleles = Seq("A", "C"),
        info    = Map("AF" -> IndexedSeq(0.5)),
        filters = Variant.PassingFilters,
        gts     = Seq(Gt("sample1", "0/1"), Gt("sample2", "./."))
      )
    }

    val vcfIn  = builder.toTempFile()
    val vcfOut = makeTempFile("prepared", ".vcf.gz")

    new PrepareVcf(input=Seq(vcfIn), output=vcfOut).execute()

    val source = VcfSource(vcfOut)
    source.header.samples.size shouldBe 0
    val variants = source.iterator.toIndexedSeq
    variants.size shouldBe 10
    variants.forall(_.genotypes.size === 0) shouldBe true
  }
}
