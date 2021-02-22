package com.editasmedicine.commons.clp

import com.fulcrumgenomics.sopt.cmdline.ClpGroup

/**
  * Grouping of command line programs that are displayed together in the usage.
  */
object ClpGroups {
  class _Alignment extends ClpGroup {
    override val name: String = "Alignment"
    override val description: String = "Tools for aligning sequences."
  }

  class _AssayDesign extends ClpGroup {
    override val name: String = "Assay Design"
    override val description: String = "Tools for performing assay designs (e.g. PCR Primer design.)"
  }

  class _Blt extends ClpGroup {
    override val name: String = "BLT"
    override val description: String = "Tools for analyzing BLT data"
  }

  class _Digenome extends ClpGroup {
    override val name: String = "Digenome"
    override val description: String = "Tools for analyzing digenome-seq data"
  }

  class _GuideRnaQc extends ClpGroup {
    override val name: String = "Guide RNA QC"
    override val description: String = "Tools for QC'ing Guide RNA sequencing experiments"
  }

  final val Alignment   = classOf[_Alignment]
  final val AssayDesign = classOf[_AssayDesign]
  final val Blt         = classOf[_Blt]
  final val Digenome    = classOf[_Digenome]
  final val GuideRnaQc  = classOf[_GuideRnaQc]
}
