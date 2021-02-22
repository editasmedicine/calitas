package com.editasmedicine.commons.clp

import com.editasmedicine.commons.LazyLogging
import com.fulcrumgenomics.cmdline.FgBioTool
import com.fulcrumgenomics.sopt.cmdline.ValidationException

import scala.collection.mutable.ListBuffer

/** Base trait which all command line tools should extend. */
trait EditasTool extends FgBioTool with LazyLogging {
  private val validationFailures = new ListBuffer[String]()

  /** Add a validation failure but don't fail just yet. */
  def addValidationError(message: String): Unit = this.validationFailures += message

  /** Throw a validation exception if any errors have been accumulated. */
  def assertValid(): Unit = if (validationFailures.nonEmpty) {
    throw new ValidationException(validationFailures.mkString("\n"))
  }
}
