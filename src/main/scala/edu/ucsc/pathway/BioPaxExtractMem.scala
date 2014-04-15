
package edu.ucsc.pathway

import org.openrdf.rio.{RDFHandler, ParserConfig, RDFParser}
import java.util.zip.GZIPInputStream
import java.io.{OutputStream, FileOutputStream, FileInputStream, File}
import org.openrdf.rio.rdfxml.{RDFXMLWriter, RDFXMLParser}
import org.openrdf.model.{Statement,Resource}
import scala.collection.mutable.HashMap
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.HashSet

object Hasher {

  val FNV_OFFSET_BASIS = 0xcbf29ce484222325L;
  val FNV_prime = 1099511628211L;

  def fnv(value : Int, seed : Int = 0x811C9DC5) = (seed ^ value) * 0x1000193

  def URLHash(url:String) : Long = {
    /*
    hash = FNV_offset_basis
    for each octet_of_data to be hashed
      hash = hash × FNV_prime
      hash = hash XOR octet_of_data
    return hash
    */
    var hash = FNV_OFFSET_BASIS;
    url.getBytes.foreach( x => {
      hash = hash * FNV_prime
      hash = hash ^ x;}
    )
    return hash;
  }
}


class BioPax_PathwaysLoad() extends RDFHandler {
  var statements = new ArrayBuffer[(Long,Long)]()
  var pathway_list = new ArrayBuffer[Long]();
  val RDF_TYPE = "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"
  val BIOPAX_PATHWAY = "http://www.biopax.org/release/biopax-level3.owl#Pathway"

  def startRDF() = {}
  def endRDF() = {}
  def handleComment(comment:String) = {}
  def handleNamespace(prefix:String, uri:String) = {}

  def handleStatement(st:Statement) = {
    val sub = st.getSubject;
    val obj = st.getObject;
    if (obj.isInstanceOf[Resource]) {
      statements += Tuple2(Hasher.URLHash( sub.toString ), Hasher.URLHash( obj.toString ))
    }
    if (st.getPredicate().toString == RDF_TYPE && st.getObject().toString == BIOPAX_PATHWAY) {
      pathway_list += Hasher.URLHash( sub.toString );
    }
  }
}




class PathwayScanner(seeds : Array[Long], edges : Array[(Long,Long)])  {
  var member_map = new HashMap[Long,HashSet[Long]]();

  for (i <- seeds) {
    member_map(i) = new HashSet[Long]()
    member_map(i) += i
  }

  def scan() {
    var added = false
    do {
      added = false
      seeds.par.foreach( group => {
        val curset = member_map(group)
        val newset = new HashSet[Long]()
        for ( (src,dst) <- edges) {
          if (curset.contains(src) && !curset.contains(dst) && !seeds.contains(dst)) {
            added = true
            newset += dst
          }
        }
        member_map(group) ++= newset
      } )
      println("Scan Cycle")
    } while (added)
  }
}

object BioPaxExtractMem {

  def parseFile(file_path:String, handler:RDFHandler) = {
    val gzip = file_path.endsWith(".gz");

    val fis = if (gzip)
      new GZIPInputStream(new FileInputStream(file_path))
    else
      new FileInputStream(file_path);
    val parser = new RDFXMLParser();
    parser.setParserConfig( new ParserConfig(true, false, false,  RDFParser.DatatypeHandling.VERIFY) );
    val baseURL = "http://pathwaycommons.org/"
    parser.setRDFHandler(handler)
    parser.parse(fis , baseURL)
    fis.close();
  }

  def main(args: Array[String]) = {
    val file_path = args(0);
    val out_dir = args(1);
    val max_files = 500;
    val pathway_scan = new BioPax_PathwaysLoad();
    parseFile(file_path, pathway_scan);
    println("Statement Count: " + pathway_scan.statements.length)

    val scanner = new PathwayScanner(pathway_scan.pathway_list.toArray, pathway_scan.statements.toArray)

    scanner.scan()

    for ( (i,v) <- scanner.member_map) {
      println(i,v)
    }

  }

}
