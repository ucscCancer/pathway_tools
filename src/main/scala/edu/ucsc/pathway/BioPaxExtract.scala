
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
  var statements = new HashMap[Long,HashSet[Long]]()
  var pathway_list = new ArrayBuffer[Long]();
  var pathway_names = new HashMap[Long,String]();

  val RDF_TYPE = "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"
  val BIOPAX_PATHWAY = "http://www.biopax.org/release/biopax-level3.owl#Pathway"

  def startRDF() = {}
  def endRDF() = {}
  def handleComment(comment:String) = {}
  def handleNamespace(prefix:String, uri:String) = {}

  def handleStatement(st:Statement) = {
    val sub = st.getSubject
    val sub_id = Hasher.URLHash( sub.toString )
    val obj = st.getObject
    if (obj.isInstanceOf[Resource]) {
      if (!statements.contains(sub_id)) {
        statements(sub_id) = new HashSet[Long]()
      }
      statements(sub_id) += Hasher.URLHash( obj.toString )
    }
    if (st.getPredicate().toString == RDF_TYPE && st.getObject().toString == BIOPAX_PATHWAY) {
      pathway_list += Hasher.URLHash(sub.toString);
      pathway_names(sub_id) = sub.toString.split("[/#]").last
    }
  }
}




class PathwayScanner(seeds : Set[Long], edges : Map[Long,HashSet[Long]])  {

  def scan() : HashMap[Long,HashSet[Long]] = {
    var member_map = new HashMap[Long,HashSet[Long]]();
    var stop_set = new HashSet[Long]()

    for (i <- seeds) {
      member_map(i) = new HashSet[Long]()
      member_map(i) += i
    }

    var added : Long = 0
    do {
      added = 0
      seeds.par.foreach( group => {
        if (!stop_set.contains(group)) {
          val curset = member_map(group)
          val newset = curset.clone()
          for ( (src,dst) <- edges) {
            if (curset.contains(src)) {
              newset ++= dst
            }
          }
          newset --= seeds
          newset += group
          if (newset.size != curset.size) {
            //println(group, curset.size, newset.size)
            this.synchronized {
              added += 1
              member_map(group) ++= newset
            }
          } else {
            this.synchronized {
              stop_set += group
            }
          }
        }
      })
      printf("Scan Cycle, %d scanning\n", added)
    } while (added > 0)
    return member_map
  }
}

class BioPaxWriter(in: OutputStream) extends RDFXMLWriter(in) {
  namespaceTable.put("http://www.biopax.org/release/biopax-level3.owl#", "bp" )
}


class BioPax_Splitter(val group_names:Map[Long,String], val elementMap:Map[Long, HashSet[Long]], val outdir:File) extends RDFHandler {
  var pathway_list = new ArrayBuffer[Resource]()
  val RDF_TYPE = "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"
  val BIOPAX_PATHWAY = "http://www.biopax.org/release/biopax-level3.owl#Pathway"

  val output_map = new HashMap[Long,FileOutputStream]()
  val writer_map = new HashMap[Long,RDFXMLWriter]()
  val group_set = new HashMap[Long,HashSet[Long]]()

  elementMap.foreach( x => x._2.foreach( y => group_set(y) = new HashSet[Long]()))
  elementMap.foreach( x => x._2.foreach( y => group_set(y) += x._1))


  def startRDF() = {
    group_names.foreach( x => {
      val f = new FileOutputStream( new File(outdir, x._2 + ".owl") );
      output_map(x._1) = f
      val w = new BioPaxWriter(f)
      w.startRDF()
      writer_map(x._1) = w
    })
  }
  def endRDF() = {
    group_names.foreach( e => {
      writer_map(e._1).endRDF()
      output_map(e._1).close()
    })
  }
  def handleComment(comment:String) = {}
  def handleNamespace(prefix:String, uri:String) = {}

  def handleStatement(st:Statement) = {
    val b = st.getSubject();
    val b_id = Hasher.URLHash(b.toString)
    if (group_set.contains(b_id)) {
      for (group <- group_set(b_id)) {
        if (group_names.contains(group)) {
          writer_map(group).handleStatement(st);
        }
      }
    }
  }


}




object BioPaxExtract {

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
    println("Statement Count: " + pathway_scan.statements.size)

    val scanner = new PathwayScanner(pathway_scan.pathway_list.toSet, pathway_scan.statements.toMap)
    val member_sets = scanner.scan()

    println("Writing output files")
    for ( i <- 0 to member_sets.size by max_files) {
      println("Writing file batch: " + i)
      val cur_names = pathway_scan.pathway_names.slice(i, i+max_files);
      val writer = new BioPax_Splitter(cur_names.toMap, member_sets.toMap, new File(out_dir));
      parseFile(file_path, writer);
    }

  }

}