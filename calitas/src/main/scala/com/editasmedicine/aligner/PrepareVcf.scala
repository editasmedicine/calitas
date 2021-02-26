package com.editasmedicine.aligner

import com.editasmedicine.commons.clp.{ClpGroups, EditasTool}
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.ProgressLogger
import com.fulcrumgenomics.vcf.api.Allele.SimpleAllele
import com.fulcrumgenomics.vcf.api.{AlleleSet, ArrayAttr, Variant, VcfContigHeader, VcfGeneralHeader, VcfSource, VcfWriter}
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor

import scala.collection.immutable.ListMap

object PrepareVcf {
  val ChromsToFix = (Range.inclusive(1, 22).map(_.toString) ++ Seq("X", "Y")).toSet
}

@clp(group=ClpGroups.Alignment, description=
  """
    |Prepares a VCF for optimal use by SearchReference.  Does the following:
    |  - Removes variants and alleles below an allele frequency threshold
    |  - Remove all INFO fields except for allele-frequency
    |  - Removes any genotypes that are present
    |  - Optionally fixes up the contig lines in the VCF
    |  - Optionally adds 'chr' prefixes to chromosome names 1-22, X and Y
    |
    |Multiple input VCFs can be given but they must be disjoint (i.e. not have variants over the same
    |regions of the genome), and wil be merged to create a single output VCF.
    |
    |
  """)
class PrepareVcf
( @arg(flag='i', doc="One or more input VCFs") val input: Seq[PathToVcf],
  @arg(flag='o', doc="The output VCF to create.") val output: PathToVcf,
  @arg(flag='f', doc="The minimum allele frequency of variants to retain.") val minAf: Double = 0.01,
  @arg(flag='d', doc="An optional sequence dictionary to use to override contig lines.") val dict: Option[FilePath] = None,
  @arg(flag='c', doc="If true, add 'chr' to chroms 1-22, X and Y.") val addChrPrefix: Boolean = true
) extends EditasTool {
  import PrepareVcf._

  override def execute(): Unit = {
    // Form up the header of the output VCF
    val header = {
      val tmp = VcfSource(input.head)
      val h   = tmp.header
      tmp.close()

      val outHeader = dict match {
        case None => h
        case Some(path) =>
          val d = SAMSequenceDictionaryExtractor.extractDictionary(path)
          h.copy(
            contigs = d.getSequences.map(s => VcfContigHeader(s.getSequenceIndex, s.getSequenceName, Some(s.getSequenceLength), Option(s.getAssembly))).toIndexedSeq,
            others = h.others.filterNot(_.headerType == "reference") ++ Seq(VcfGeneralHeader("reference", d.getSequence(0).getAssembly))
          )
      }

      outHeader.copy(samples = IndexedSeq.empty)
    }

    val out = VcfWriter(output, header)
    val progress = ProgressLogger(logger, noun="variants", verb="wrote")

    input.foreach { vcfIn =>
      logger.info(s"Processing file $vcfIn")
      val in = VcfSource(vcfIn)

      in
        .filter(v => v.filters == Variant.PassingFilters)
        .filter(v => v[ArrayAttr[Float]]("AF").exists(_ >= minAf))
        .filter(v => v.alleles.forall(a => a.isInstanceOf[SimpleAllele]))
        .foreach { v =>
          val altsAndAfs = v.alleles.alts.zip(v[ArrayAttr[Float]]("AF")).filter(_._2 >= minAf)
          val fixed = v.copy(
            chrom     = if (addChrPrefix) fixChrom(v.chrom) else v.chrom,
            alleles   = AlleleSet(v.alleles.ref, altsAndAfs.map(_._1)),
            genotypes = Map.empty,
            attrs     = ListMap("AF" -> ArrayAttr(altsAndAfs.map(_._2)))
          )

          out += fixed
          progress.record(v.chrom, v.pos)
        }

      in.safelyClose()
    }

    out.close()
  }

  /** Prepends "chr" to chromosome names that need it. */
  def fixChrom(chrom: String): String = if (ChromsToFix.contains(chrom)) s"chr${chrom}" else chrom
}
