package com.editasmedicine.commons

import com.fulcrumgenomics.commons.util.LogLevel

/** Thin extension of logger class. */
class Logger(clazz: Class[_]) extends com.fulcrumgenomics.commons.util.Logger(clazz) {
  /** Logs a large banner message to the screen with a border. */
  def banner(level: LogLevel, lines: String*): Unit = {
    val width = lines.map(_.length).max + 4
    val line = "#" * width
    emit(level, line)
    lines.foreach(l => emit(level, Seq(" #", l, " #")))
    emit(level, line)
  }
}

/** Trait that provides a lazily created logger. */
trait LazyLogging {
  protected lazy val logger = new Logger(getClass)
}
