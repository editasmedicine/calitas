package com.editasmedicine.aligner

import com.editasmedicine.aligner.SequentialGuideAligner.{Defaults, Guide}
import com.editasmedicine.commons.clp.{ClpGroups, EditasTool}
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.Io

@clp(group=ClpGroups.Alignment, description=
  """
    |Performs pairwise alignment of sequences.  Input is a file with two sequences per
    |line, separated by whitespace.  Sequences may be composed of DNA and RNA bases and
    |ambiguity codes.  No headers are required or expected.
    |
    |Designed for performing glocal alignment of guides to larger sequences (i.e. global
    |alignment of the query to a local portion of the target.)
    |
    |The query sequence may be a mixture of upper and lower case sequence, with lower-case bases used for the
    |PAM at one end of the sequence or the other.
    |
    |The target sequences (second on each line) should be all upper case and will be made upper-case
    |if they are not already.
  """)
class PairwiseAlignSequences
(@arg(flag='i', doc="Input file of sequence pairs.") val input: FilePath,
 @arg(flag='o', doc="Output file to write.") val output: FilePath = Io.StdOut,
 @arg(flag='t', doc="Threads to use for alignments.") val threads: Int = 8,
 @arg(flag='g', doc="Maximum gap bases between guide and PAM") val maxGapsBetweenGuideAndPam: Int = Defaults.MaxGapsBetweenGuideAndPam,
 @arg(flag='O', doc="Maximum overlap allowed between alignments on the same strand.") val maxOverlap: Int = Defaults.MaxOverlap,
 @arg(flag='m', doc="Net cost of going from a match to a mismatch in the guide.") val guideMismatchNetCost: Int = Defaults.MismatchNetCost,
 @arg(flag='M', doc="Net cost of going from a match to a mismatch in the PAM.") val pamMismatchNetCost: Int = Defaults.PamMismatchNetCost,
 @arg(flag='b', doc="Net cost of a 1bp gap in the genomic/target sequence.") val genomeGapNetCost: Int = Defaults.GenomeGapNetCost,
 @arg(flag='B', doc="Net cost of a 1bp gap in the guide.") val guideGapNetCost: Int = Defaults.GuideGapNetCost,

) extends EditasTool {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  case class Task(query: String, target: String)

  override def execute(): Unit = {
    logger.info("Reading input file.")
    val tasks = Io.readLines(input).map(_.trim).filter(_.nonEmpty).map(_.split("""\s+""")).map { fields =>
      require(fields.length == 2, s"Line found with ${fields.length} fields: ${fields.mkString(" ")}")
      Task(fields(0), fields(1).toUpperCase)
    }

    logger.info("Beginning alignments.")
    val aligner  = new SequentialGuideAligner(
      mismatchNetCost    = guideMismatchNetCost,
      pamMismatchNetCost = pamMismatchNetCost,
      genomeGapNetCost  = genomeGapNetCost,
      guideGapNetCost = guideGapNetCost
    )

    val out = Io.toWriter(output)
    val cols = Seq("query", "target", "score", "query_start", "target_start", "cigar",
    "mismatches", "gap_bases", "padded_query", "alignment", "padded_target")
    out.write(cols.mkString("", "\t", "\n"))

    while (tasks.hasNext) {
      val batch = tasks.take(10000).toIndexedSeq
      val results = batch.parWith(threads).map { task => (task, aligner.alignBest(Guide(task.query), task.target.getBytes)) }.seq

      results.foreach { case (task, aln) =>
        val fields = Seq(
          task.query,
          task.target,
          aln.score,
          1,
          aln.startOffset,
          aln.cigar.toString(),
          aln.mismatches,
          aln.gapBases,
          aln.paddedGuide,
          aln.paddedAlignment,
          aln.paddedTarget
        )

        out.write(fields.mkString("", "\t", "\n"))
      }
    }

    out.close()
  }
}
