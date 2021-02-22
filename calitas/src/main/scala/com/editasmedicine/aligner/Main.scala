package com.editasmedicine.aligner

import com.editasmedicine.commons.clp.ClpMain

/**
  * Main program for editas-tools that loads everything up and runs the appropriate sub-command
  */
object Main {
  /** The main method */
  def main(args: Array[String]): Unit = new ClpMain("calitas").makeItSoAndExit(args)
}
