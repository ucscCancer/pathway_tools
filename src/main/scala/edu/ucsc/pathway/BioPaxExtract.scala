
package edu.ucsc.pathway

import org.openrdf.rio.{RDFHandler, ParserConfig, RDFParser}
import java.util.zip.GZIPInputStream
import java.io.{FileOutputStream, FileInputStream, File}
import org.openrdf.rio.rdfxml.{RDFXMLWriter, RDFXMLParser}
import org.openrdf.model.{Statement,Resource}
import scala.collection.mutable.HashMap
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.HashSet
import scala.collection.mutable


class BioPax_Pathways() extends RDFHandler {
  var pathway_list = new ArrayBuffer[Resource]()
  val RDF_TYPE = "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"
  val BIOPAX_PATHWAY = "http://www.biopax.org/release/biopax-level3.owl#Pathway"

  def startRDF() = {}
  def endRDF() = {}
  def handleComment(comment:String) = {}
  def handleNamespace(prefix:String, uri:String) = {}

  def handleStatement(st:Statement) = {
    val b = st.getSubject();
    if (st.getPredicate().toString == RDF_TYPE && st.getObject().toString == BIOPAX_PATHWAY) {
      pathway_list += b;
    }
  }
}

class NetworkGrouper(seeds : Array[Resource]) extends RDFHandler {
  val original_set = new HashSet[Resource]() ++ seeds;
  var member_map = new HashMap[Resource,HashSet[Int]]();
  var added = false;
  var name_map = new HashMap[String,Int]();

  seeds.zipWithIndex.foreach( x => {
    member_map(x._1) = member_map.getOrElse(x._1, new HashSet[Int]()) ++ List(x._2);
    name_map(x._1.toString().split("[/#]").last) = x._2
  } )


  def startRDF() = {
    added = false
  }
  def endRDF() = {}
  def handleComment(comment:String) = {}
  def handleNamespace(prefix:String, uri:String) = {}

  def handleStatement(st:Statement) = {
    val obj = st.getObject();
    if (obj.isInstanceOf[Resource]) {
      val subj = st.getSubject();
      if (!original_set.contains(obj.asInstanceOf[Resource])) {
        if (member_map.contains(subj)) {
          if (!member_map.contains(obj.asInstanceOf[Resource])) {
            member_map(obj.asInstanceOf[Resource]) = new HashSet[Int]() ++ member_map(subj);
          } else {
            val o_set = member_map(obj.asInstanceOf[Resource]);
            val s_set = member_map(subj);
            if (o_set.intersect(s_set) != s_set) {
              o_set ++= s_set;
              added = true;
            }
          }
        }
      }
    }
  }
}




class BioPax_Splitter(val groups:Map[String,Int], val elementMap:Map[Resource, HashSet[Int]], val outdir:File) extends RDFHandler {
  var pathway_list = new ArrayBuffer[Resource]()
  val RDF_TYPE = "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"
  val BIOPAX_PATHWAY = "http://www.biopax.org/release/biopax-level3.owl#Pathway"

  val output_map = new HashMap[Int,FileOutputStream]();
  val writer_map = new HashMap[Int,RDFXMLWriter]();

  val group_set = groups.values.toSet;

  def startRDF() = {
    groups.foreach( x => {
      val f = new FileOutputStream( new File(outdir, x._1 + ".owl") );
      output_map(x._2) = f;
      val w = new RDFXMLWriter(f);
      w.startRDF();
      writer_map(x._2) = w;
    })
  }
  def endRDF() = {
    groups.foreach( e => {
      writer_map(e._2).endRDF();
      output_map(e._2).close();
    })
  }
  def handleComment(comment:String) = {}
  def handleNamespace(prefix:String, uri:String) = {}

  def handleStatement(st:Statement) = {
    val b = st.getSubject();
    val b_sets = elementMap.get(b);
    if (b_sets.isDefined) {
      for ( group <- b_sets.get) {
        if ( group_set.contains(group)) {
          if (b_sets.get.contains(group)) {
            writer_map(group).handleStatement(st);
          }
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
    val pathway_scan = new BioPax_Pathways();
    parseFile(file_path, pathway_scan);
    println("Pathway Count: " + pathway_scan.pathway_list.length)

    val ng = new NetworkGrouper(pathway_scan.pathway_list.toArray)

    do {
      parseFile(file_path, ng);
      println("Selection Cycle: " + ng.member_map.size + " elements");
    } while (ng.added);

    println("Writing output files")
    for ( i <- 0 to ng.name_map.size by max_files) {
      println("Writing file batch: " + i)
      val cur_names = ng.name_map.slice(i, i+max_files);
      val writer = new BioPax_Splitter(cur_names.toMap, ng.member_map.toMap, new File(out_dir));
      parseFile(file_path, writer);
    }
  }


}
