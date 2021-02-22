import java.text.SimpleDateFormat
import java.util.Date

import com.typesafe.sbt.SbtGit.GitCommand
import sbt.Keys._
import sbt._
import sbtassembly.AssemblyKeys.assembly
import scoverage.ScoverageKeys._

import scala.sys.process.Process

////////////////////////////////////////////////////////////////////////////////////////////////
// Multi-project build file for the following projects:
// - commons: utility code and base classes used across projects
// - calitas: code for the CALITAS aligner
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// Common settings for all projects
////////////////////////////////////////////////////////////////////////////////////////////////
val htmlReportsDirectory: String = "target/test-reports"

lazy val htsjdkAndPicardExcludes = Seq(
  ExclusionRule(organization="org.apache.ant"),
  ExclusionRule(organization="gov.nih.nlm.ncbi"),
  ExclusionRule(organization="org.testng"),
  ExclusionRule(organization="com.google.cloud.genomics")
)

lazy val ToolkitVersion = {
  val date     = new SimpleDateFormat("yyyyMMdd").format(new Date())
  val hash     = Process("git rev-parse --short HEAD").lineStream.head
  val modified = Process("git status --porcelain").lineStream.nonEmpty

  date + "-" + hash + (if (modified) "-dirty" else "")
}

version in ThisBuild := ToolkitVersion

lazy val commonSettings = Seq(
  organization         := "com.editasmedicine",
  organizationName     := "Editas Medicine Inc.",
  organizationHomepage := Some(url("http://www.editasmedicine.com/")),
  homepage             := Some(url("https://bitbucket.org/editascomputationalbiology/")),
  startYear            := Some(2017),
  scalaVersion         := "2.12.8",
  scalacOptions        += "-target:jvm-1.8",
  autoAPIMappings      := true,
  version              := ToolkitVersion,
  testOptions in Test  += Tests.Argument(TestFrameworks.ScalaTest, "-h", Option(System.getenv("TEST_HTML_REPORTS")).getOrElse(htmlReportsDirectory)),
  testOptions in Test  += Tests.Argument("-oDF"),
  resolvers            += Resolver.jcenterRepo,
  resolvers            += Resolver.sonatypeRepo("public"),
  resolvers            += Resolver.mavenLocal,
  resolvers            += Resolver.file("ivy-local", new File("~/.ivy2/local")),
  resolvers            += "BroadSnapshots" at "https://broadinstitute.jfrog.io/broadinstitute/libs-snapshot-local/",
  shellPrompt          := { state => "%s| %s> ".format(GitCommand.prompt.apply(state), version.value) },
  updateOptions        := updateOptions.value.withCachedResolution(true),
  javaOptions in Test  += "-Xmx1G",
  fork in Test         := true
) ++ Defaults.coreDefaultSettings

lazy val assemblySettings = Seq(
  // Assembly related settings
  assemblyJarName in assembly := s"${name.value}.jar",
  test in assembly := {},
  assemblyMergeStrategy in assembly := {
    case PathList("log4j2.xml") => MergeStrategy.first
    case x =>
      val oldStrategy = (assemblyMergeStrategy in assembly).value
      oldStrategy(x)
  }
)

////////////////////////////////////////////////////////////////////////////////////////////////
// commons project
////////////////////////////////////////////////////////////////////////////////////////////////
lazy val commons = Project(id="commons", base=file("commons"))
  .settings(commonSettings: _*)
  .settings(description := "Utility and base classes used across projects.")
  .settings(libraryDependencies ++= Seq(
      "org.scalatest"       %% "scalatest"     % "3.0.5" % "test->*" excludeAll ExclusionRule(organization="org.junit", name="junit"),
      "com.fulcrumgenomics" %% "sopt"          % "0.7.0",
      "com.fulcrumgenomics" %% "commons"       % "0.7.0",
      "com.fulcrumgenomics" %% "fgbio"         % "1.1.0-4f11a43-SNAPSHOT" excludeAll(htsjdkAndPicardExcludes:_*),
      "com.beachape"        %% "enumeratum"    % "1.5.12",
      "com.github.samtools"  % "htsjdk"        % "2.19.0-26-g1e3f1fa-SNAPSHOT" excludeAll(htsjdkAndPicardExcludes: _*)
      ),
  )
  .disablePlugins(AssemblyPlugin)


////////////////////////////////////////////////////////////////////////////////////////////////
// Aligner project
////////////////////////////////////////////////////////////////////////////////////////////////
lazy val calitas = Project(id="calitas", base=file("calitas"))
  .settings(commonSettings: _*)
  .settings(assemblySettings: _*)
  .settings(description := "Library and command line tools for CALITAS.")
  .dependsOn(commons % "compile->compile;test->test")


lazy val root = Project(id="calitas-and-commons", base=file("."))
  .settings(commonSettings:_*)
  .aggregate(commons, calitas)
  .disablePlugins(AssemblyPlugin)
