package edu.ucsc.pathway

import java.io.FileInputStream
import org.biopax.paxtools.model.level3.{Complex, Catalysis, ProteinReference, BiochemicalReaction, Pathway}
import org.biopax.paxtools.pattern.constraint.ConBox

import collection.JavaConverters._

import org.biopax.paxtools.io.SimpleIOHandler
import org.biopax.paxtools.pattern.{Searcher, Pattern}


object BioPaxConvert {

  def stringClean(in:String) : String = {
    in.replace(" ", "_").replace("/", "_").replace("(", "_").replace(")", "_")
  }

  def pathway_name(pathway:Pathway) : String = {
    if (pathway.getStandardName != null)
      stringClean(pathway.getStandardName)
    else
      stringClean(pathway.getName.asScala.last)
  }

  def main(args:Array[String]) = {
    val fin = new FileInputStream(args(0))
    val handler = new SimpleIOHandler()
    val model = handler.convertFromOWL(fin)

    /*
    model.getObjects(classOf[ProteinReference]).asScala.foreach( x => {
      println("==Protein==")
      println(x.getXref)
      println(x.getEntityReferenceOf)
    })

    model.getObjects(classOf[Complex]).asScala.foreach( x => {
      println("==Complex==")
      println(x.getXref)
      println(x.getComponent)
    })


    model.getObjects(classOf[Catalysis]).asScala.foreach( x => {
      println("==Catalysis==")
      println(x.getXref)
      println(x.getController)
      println(x.getControlled)
      println(x.getCofactor)
      println(x.getInteractionType)
      println(x.getPathwayComponentOf)
    })

    model.getObjects(classOf[BiochemicalReaction]).asScala.foreach( x => {
      println("==BiochemicalReaction==")
      println(x.getXref)
      println(x.getLeft)
      println(x.getRight)
      println(x.getControlledOf)
      println(x.getInteractionType)
      println(x.getPathwayComponentOf)
    })
  */

    val pathway_list = model.getObjects(classOf[Pathway]).asScala.map( x => {
      (pathway_name(x),x)
    })

    println("Pathway Count", model.getObjects(classOf[Pathway]).asScala.size)
    println(pathway_list.size)

    var pathway_map = pathway_list.groupBy(_._1).flatMap( x => {
      if (x._2.size == 1)
        x._2.toSeq
      else
        x._2.zipWithIndex.map(y => ("%s_%d".format(y._1._1, y._2), y._1._2) )
    }).toMap


    val pathway_filter = pathway_map.toSeq.map( x => {
      (x._1, x._2, x._2.getPathwayComponent().asScala.seq)
    }).filter( x => x._3.size > 0 )

    println("mapsize:", pathway_map.size)
    println("filtersize:", pathway_filter.size)

    pathway_filter.map( x => {
      x._3.toSeq.map( y => {
        y.getModelInterface.toString
      }).mkString(",")
    }).foreach(println)

    /*
    val pattern = new Pattern(classOf[BiochemicalReaction], "react")
    pattern.add(ConBox.left(), "react", "left")
    pattern.add(ConBox.right(), "react", "right")

    Searcher.search(model, pattern).asScala.foreach( x => {
      println("Element", x._1)
      //println("Matches", x._2)
      x._2.asScala.map( y => println("Left:", y.get("left", pattern) ) )
    })
    */


  }
}
