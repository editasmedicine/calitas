package com.editasmedicine.commons.testing

import java.nio.file.{Files, Path}
import java.util.stream.Collectors

import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.io.PathUtil
import org.scalatest.{FlatSpec, Matchers}

/** Base class for unit tests. */
class UnitSpec extends FlatSpec with Matchers {
  /** Creates a new temp file for use in testing that will be deleted when the VM exits. */
  protected def makeTempFile(prefix: String, suffix: String) : Path = {
    val path = Files.createTempFile(prefix, suffix)
    path.toFile.deleteOnExit()
    path
  }

  /** Creates a new temporary directory for use. */
  protected def makeTempDir(prefix: String, suffix: String) : Path = {
    val path = Files.createTempFile(prefix, suffix)
    Files.delete(path)
    Files.createDirectory(path)
    path.toFile.deleteOnExit()
    path
  }

  /** Deletes a path/file, recursing down if the path is directory. */
  protected def deletePath(path: Path): Unit = {
    if (Files.isDirectory(path)) {
      val childStream = Files.list(path)
      val children = childStream.collect(Collectors.toList())
      childStream.close()
      children.foreach(deletePath)
    }

    Files.deleteIfExists(path)
  }

  /** Executes a block of code with access to a temp directory that is cleaned up after the block. */
  protected def withTempDir[A](body: DirPath => A): A = {
    val dir = makeTempDir("testing.", ".dir")
    try { body(dir) }
    finally { deletePath(dir) }
  }

  /** Finds a resource on the classpath and turns it into a Path. */
  protected def testResource(name: String, pkg: Package = getClass.getPackage) : Path = {
    val packagePath = pkg.getName.replace('.', '/')
    val url = getClass.getClassLoader.getResource(packagePath + "/" + name)
    PathUtil.pathTo(url.getPath)
  }
}
