package com.editasmedicine.aligner

import java.io.FileInputStream
import java.security.{DigestInputStream, MessageDigest}
import java.text.SimpleDateFormat
import java.util.{Date, TimeZone}

import com.editasmedicine.aligner.SearchReferencesWithVariants.VariantAllele
import com.editasmedicine.aligner.SequentialGuideAligner.Guide
import com.editasmedicine.commons.EditasMetric
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.util.{Metric, Sequences}
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.reference.{ReferenceSequenceFile, ReferenceSequenceFileFactory}
import htsjdk.samtools.util.CoordMath


/**
  * Represents an alignment of a guide (with or without PAM) to a reference genome.
  *
  * Note: `gaps` in this context refers to gapped bases.  I.e. a single 2bp gap is counted as 2, not 1.
  */
trait StandardAlignmentOutput {
  /** Name of guide, typically this will be the RSQ####. */
  val guide_id: String
  /** the sequence of the guide used, unpadded (no gaps/bulges)*/
  val unpadded_guide_sequence: String
  /** eg hg38, mm6; for variants use the genome which was used to build the vcf */
  val genome_build: String
  /** Chromosome for target sequence alignment (eg: chr3) */
  val chromosome: String
  /** start of the unpadded target sequence in the genome, using 0-based open ended, does not include PAM */
  val coordinate_start: Int
  /** end of the unpadded target sequence in the genome,using 0-based open ended, does not include PAM */
  val coordinate_end: Int
  /** “+” or “-”, the target sequence strand */
  val strand: String
  /** the unpadded target sequence (as DNA) as found in the genome, without gaps/bulges and does not include PAM */
  val unpadded_target_sequence: String
  /** 10 bases before 5' (5 prime) end of target sequence in the genome, it may include the PAM; note that this is strand sensitive so 5' on the '+' strand is left of target sequence and on the '-' strand is to the right of the target sequence */
  val ten_bases_5_prime: String
  /** 10 bases before 3' (3 prime) end of target sequence in the genome, it may include the PAM; note that this is strand sensitive so 3' on the '+' strand is right of target sequence and on the '-' strand is to the left of the target sequence */
  val ten_bases_3_prime: String
  /** pam used in the alignment (eg: NNGRRN) */
  val pam_used: Option[String]
  /** if there is variant ID associated with the alignment then report it here, can be null (eg: rs1234567);  separate multiple variants with a semicolon */
  val variant_id: Option[String]
  /** */
  val variant_description: Option[String]
  /** variant vcf file name + version if available + hash(md5) */
  val variant_vcf: Option[String]
  /** global allele frequency for a variant. given as a fraction (eg 0-1) [real]; for multiple variants take the minimum of the frequency (this represents the maximal allele frequency) */
  val allele_frequency: Option[Double]
  /** final alignment score (including PAM) [int] */
  val score: Int
  /** Mismatches in the guide region [int] */
  val guide_mm: Int
  /** Gaps in the guide region [int] */
  val guide_gaps: Int
  /** Total gaps and mismatches in the guide region [int] */
  val guide_mm_plus_gaps: Int
  /** Mismatches in the PAM [int] */
  val pam_mm: Int
  /** total count of mismatches and gaps across the both Guide and Pam regions [int] */
  val total_mm_plus_gaps: Int
  /** Guide + PAM sequence with padding for mismatches and bulges/gaps; if multiple PAMs used use the PAM which was selected in the alignment */
  val padded_guide: String
  /** Visual representation of the guide-target alignment: | for matches, . for mismatches, ~ for gaps */
  val padded_alignment: String
  /** Target sequence +PAM sequence with padding for mismatches and bulges/gaps */
  val padded_target: String
  /** an additional 8 bases on the 5' side of the padded_target (past the PAM if applicable), used to help visualize the context of the alignment */
  val padded_extra_8_bases_5_prime: String
  /** an additional 8 bases on the 3' side of the padded_target (past the PAM if applicable), used to help visualize the context of the alignment */
  val padded_extra_8_bases_3_prime: String
  /** Cigar representation of Guide sequence alignment */
  val cigar: String
  /** Length of guide sequence not including PAM, nor gaps/bulges; this is the length(guide_sequence_used) */
  val unpadded_guide_sequence_length: Int
  /** Length of target sequence not including PAM, nor gaps/bulges; his is the length(bases) */
  val unpadded_target_sequence_length: Int
  /** aligner name (glsearch36/SearchReference) */
  val aligner: String = "SearchReference"
  /** version number of build id for aligner */
  val aligner_version: String = EditasMetric.Version
  /** commas separated list of pams used in the search */
  val aligner_search_pam: String
  /**  a semicolon or comma separated list of parameters used for this alignment (may need to discuss) */
  val aligner_other_parameters: String
  /** a date_time stamp for when the alignment run was started (UTC) in this format: Fri Sep 27 08:58:29 MDT 2019  */
  val time_stamp: String
}


/**
  * Instantiation of the `StandardAlignmentOutput` trait.
  */
case class ReferenceHit(override val guide_id: String,
                        override val unpadded_guide_sequence: String,
                        override val genome_build: String,
                        override val chromosome: String,
                        override val coordinate_start: Int,
                        override val coordinate_end: Int,
                        override val strand: String,
                        override val unpadded_target_sequence: String,
                        override val ten_bases_5_prime: String,
                        override val ten_bases_3_prime: String,
                        override val pam_used: Option[String],
                        override val variant_id: Option[String],
                        override val variant_description: Option[String],
                        override val variant_vcf: Option[String],
                        override val allele_frequency: Option[Double],
                        override val score: Int,
                        override val guide_mm: Int,
                        override val guide_gaps: Int,
                        override val guide_mm_plus_gaps: Int,
                        override val pam_mm: Int,
                        override val total_mm_plus_gaps: Int,
                        override val padded_guide: String,
                        override val padded_alignment: String,
                        override val padded_target: String,
                        override val padded_extra_8_bases_5_prime: String,
                        override val padded_extra_8_bases_3_prime: String,
                        override val cigar: String,
                        override val unpadded_guide_sequence_length: Int,
                        override val unpadded_target_sequence_length: Int,
                        override val aligner: String,
                        override val aligner_version: String = EditasMetric.Version,
                        override val aligner_search_pam: String,
                        override val aligner_other_parameters: String,
                        override val time_stamp: String) extends StandardAlignmentOutput with Metric {

  /** Computes the end of the alignment from the start position plus the cigar. */
  lazy val end: Int = {
    val len = Cigar(cigar).lengthOnTarget
    CoordMath.getEnd(coordinate_start, len)
  }

  /** Computes the overlap between the alignments of two reference hits on the genome. */
  def overlap(that: ReferenceHit): Int = {
    if (that.chromosome != this.chromosome) 0
    else math.max(0, math.min(this.end, that.end) - math.max(this.coordinate_start, that.coordinate_start))
  }
}

object ReferenceHit {
  object Builder {
    /**
      * Generates a new ReferenceHit.Builder object that can be used to build ReferenceHit instances.
      * Does most of the heavy lifting here so that the resulting Builder can be `.copy()`'d with new
      * guides and guide IDs without reloading the reference or re-generating the MD5 of the variant VCF.
      *
      * @param guideId the ID for the guide being used
      * @param guide the Guide object with the guide sequence and PAM(s)
      * @param ref the path to the reference genome fasta file
      * @param vcf an optional path to a VCF of variants that was used during alignment
      * @param alignerId the name or ID of the alignment software being used
      * @param arguments any relevant arguments to the aligner
      */
    def apply(guideId: String,
              guide: Guide,
              ref: PathToFasta,
              vcf: Option[PathToVcf],
              alignerId: String,
              arguments: String): Builder = {

      val refFile          = ReferenceSequenceFileFactory.getReferenceSequenceFile(ref)
      val timestamp        = {
        val fmt = new SimpleDateFormat("EEE MMM dd HH:mm:ss z yyyy")
        fmt.setTimeZone(TimeZone.getTimeZone("UTC"))
        fmt.format(new Date())
      }

      val variantVcf = vcf.map { path =>
        val buffer = new Array[Byte](64 * 1024)
        val digest = MessageDigest.getInstance("MD5")
        val stream = new DigestInputStream(new FileInputStream(path.toFile), digest)
        while (stream.read(buffer) != -1) { () }
        stream.close()
        val md5 = digest.digest.map("%02x".format(_)).mkString
        s"${path.getFileName}:$md5"
      }

      new Builder(
        guideId = guideId,
        guide   = guide,
        refFile = refFile,
        vcfId   = variantVcf,
        alignerId = alignerId,
        timestamp = timestamp,
        arguments = arguments
      )
    }
  }

  /** A builder object to make it easier to make ReferenceHit objects. */
  case class Builder private (guideId: String,
                              guide: Guide,
                              refFile: ReferenceSequenceFile,
                              vcfId: Option[String],
                              alignerId: String,
                              timestamp: String,
                              arguments: String
                             ) {

    private val alignerSearchPam = (guide.pams5Prime ++ guide.pams3Prime).mkString(",")
    private val genomeBuild      = refFile.getSequenceDictionary.getSequences.iterator().flatMap(s => Option(s.getAssembly)).find(_ => true).getOrElse("unknown")

    def build(aln: GuideAlignment, variants: Seq[VariantAllele] = Seq.empty): ReferenceHit = {
      val vs = variants.filter(v => v.pos-1 >= aln.startOffset && v.pos-1 <= aln.endOffset)

      lazy val tenLeft    = fetchBases(refFile, aln.chrom, aln.guideStartOffset+1-10, aln.guideStartOffset,  rc=aln.isNegativeStrand)
      lazy val tenRight   = fetchBases(refFile, aln.chrom, aln.guideEndOffset+1,      aln.guideEndOffset+10, rc=aln.isNegativeStrand)
      lazy val eightLeft  = fetchBases(refFile, aln.chrom, aln.startOffset+1-8,       aln.startOffset,       rc=aln.isNegativeStrand)
      lazy val eightRight = fetchBases(refFile, aln.chrom, aln.endOffset+1,           aln.endOffset+8,       rc=aln.isNegativeStrand)

      ReferenceHit(
        guide_id                        = guideId,
        unpadded_guide_sequence         = guide.guide,
        genome_build                    = if (vs.isEmpty) genomeBuild else s"${genomeBuild}+variants",
        chromosome                      = aln.chrom,
        coordinate_start                = aln.guideStartOffset,
        coordinate_end                  = aln.guideEndOffset,
        strand                          = aln.strand.toString,
        unpadded_target_sequence        = aln.unpaddedTargetWithoutPam,
        ten_bases_5_prime               = aln.leftOfGuide10bp.getOrElse(if (aln.isPositiveStrand) tenLeft else tenRight),
        ten_bases_3_prime               = aln.rightOfGuide10bp.getOrElse(if (aln.isPositiveStrand) tenRight else tenLeft),
        pam_used                        = Option(aln.guide.filter(_.isLower)).filter(_.nonEmpty),
        variant_id                      = if (vs.nonEmpty) Some(vs.map(_.id).mkString(";")) else None,
        variant_description             = if (vs.nonEmpty) Some(vs.map(_.displayString).mkString(";")) else None,
        variant_vcf                     = if (vs.nonEmpty) vcfId else None,
        allele_frequency                = if (vs.nonEmpty) Some(vs.minBy(_.af).af) else None,
        score                           = aln.score,
        guide_mm                        = aln.guideMismatches,
        guide_gaps                      = aln.guideGapBases,
        guide_mm_plus_gaps              = aln.guideMmsPlusGaps,
        pam_mm                          = aln.pamMismatches,
        total_mm_plus_gaps              = aln.edits,
        padded_guide                    = aln.paddedGuide,
        padded_alignment                = aln.paddedAlignment,
        padded_target                   = aln.paddedTarget,
        padded_extra_8_bases_5_prime    = aln.leftOfFullAln8bp.getOrElse(if (aln.isPositiveStrand) eightLeft  else eightRight),
        padded_extra_8_bases_3_prime    = aln.rightOfFullAln8bp.getOrElse(if (aln.isPositiveStrand) eightRight else eightLeft),
        cigar                           = aln.cigar.toString,
        unpadded_guide_sequence_length  = guide.guide.length,
        unpadded_target_sequence_length = aln.unpaddedTargetWithoutPam.length,
        aligner                         = alignerId,
        aligner_version                 = EditasMetric.Version,
        aligner_search_pam              = alignerSearchPam,
        aligner_other_parameters        = arguments,
        time_stamp                      = timestamp
      )
    }
  }

  /**
    * Fetches bases from a sequence, padding with Ns if the requested sequence is too close to the end.
    * Start and end should be provided 1-based inclusive.
    */
  private def fetchBases(ref: ReferenceSequenceFile, chrom: String, start: Int, end: Int, rc: Boolean) = {
    val adjustedStart = math.max(1, start)
    val adjustedEnd   = math.min(ref.getSequenceDictionary.getSequence(chrom).getSequenceLength, end)
    val bases = ("N" * (adjustedStart - start)) + ref.getSubsequenceAt(chrom, adjustedStart, adjustedEnd).getBaseString + ("N" * (end - adjustedEnd))
    if (rc) Sequences.revcomp(bases).toUpperCase else bases.toUpperCase
  }

  /** Sorts a sequence of reference hits by chromosome and position. If a sequence dictionary
    * or a reference that has a sequence dictionary is provided then the chromosome ordering will
    * match the reference.  Otherwise the chromosome ordering will be lexicographic.
    *
    * @param hits the seq of hits to be sorted
    * @param ref an optional reference
    * @return the sorted hits
    */
  def sort(hits: Seq[ReferenceHit], ref: Either[SAMSequenceDictionary,PathToFasta]): Seq[ReferenceHit] = {
    val dict = ref match {
      case Left(sd) => Some(sd)
      case Right(r) => Option(ReferenceSequenceFileFactory.getReferenceSequenceFile(r).getSequenceDictionary)
    }

    // TODO: add another criteria for minimizing variants
    dict match {
      case Some(sd) => hits.sortBy(h => (sd.getSequenceIndex(h.chromosome), h.coordinate_start, h.strand, -h.score))
      case None     => hits.sortBy(h => (h.chromosome, h.coordinate_start, h.strand, -h.score))
    }
  }
}
