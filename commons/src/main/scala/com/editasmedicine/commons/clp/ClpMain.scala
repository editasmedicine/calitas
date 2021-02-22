package com.editasmedicine.commons.clp

import com.fulcrumgenomics.cmdline.FgBioMain

/** Main class used to hook commands up to the command line. */
class ClpMain(val toolkit: String, val packages: List[String] = "com.editasmedicine" :: Nil) extends FgBioMain {
  override protected def name: String = toolkit
  override protected def packageList: List[String] = packages
}
