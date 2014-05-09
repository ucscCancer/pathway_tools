name := "Pathway Tools"

version := "0.0.1"

scalaVersion := "2.10.0"

libraryDependencies += "org.openrdf.sesame" % "sesame-rio" % "2.7.10"

libraryDependencies += "org.openrdf.sesame" % "sesame-rio-ntriples" % "2.7.10"

libraryDependencies += "org.openrdf.sesame" % "sesame-rio-rdfxml" % "2.7.10"


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
