package com.editasmedicine.commons

import java.text.SimpleDateFormat
import java.util.Date

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.util.Metric

object EditasMetric {
  private val DateFormat = new SimpleDateFormat("yyyy-MM-dd")

  /** The version of the software. */
  val Version: String = Option(getClass.getPackage.getImplementationVersion).getOrElse(s"unknown-${format(new Date)}")

  /** Formats a date into a standard yyyy-MM-dd format. */
  def format(d: Date): String = DateFormat.synchronized(DateFormat.format(d))
}

/** Base trait to ensure that all metrics include some standard information. */
trait EditasMetric extends Metric {
  def analysis_date: String
  def pipeline_version: String

  override protected def formatValue(value: Any): String = value match {
    case 0d | 0f | 0 | 0L => "0"
    case x  => super.formatValue(x)
  }
}
