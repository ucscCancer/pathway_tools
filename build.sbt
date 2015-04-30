name := "Pathway Tools"

version := "0.0.1"

scalaVersion := "2.10.4"



libraryDependencies += "org.openrdf.sesame" % "sesame-rio" % "2.7.10"

libraryDependencies += "org.openrdf.sesame" % "sesame-rio-ntriples" % "2.7.10"

libraryDependencies += "org.openrdf.sesame" % "sesame-rio-rdfxml" % "2.7.10"

libraryDependencies += "org.biopax.paxtools" % "pattern" % "2.0.0-SNAPSHOT"

resolvers ++= Seq(
  "Akka Repository" at "http://repo.akka.io/releases/",
  "Sonatype Snapshots" at "http://oss.sonatype.org/content/repositories/snapshots",
  "Sonatype Releases" at "http://oss.sonatype.org/content/repositories/releases",
  "Repository of BioPAX Pattern Framework" at "http://maven-repo.biopax-pattern.googlecode.com/hg/",
  "BioPax Releases" at "http://www.biopax.org/m2repo/releases/"
)



TaskKey[File]("mkrun") <<= (baseDirectory, fullClasspath in Runtime, mainClass in Runtime) map { (base, cp, main) =>
  val template = """#!/bin/sh
java -Xmx2g -classpath "%s" %s "$@"
"""
  val mainStr = main getOrElse error("No main class specified")
  val contents = template.format(cp.files.absString, mainStr)
  val out = base / "run-pathtools.sh"
  IO.write(out, contents)
  out.setExecutable(true)
  out
}
