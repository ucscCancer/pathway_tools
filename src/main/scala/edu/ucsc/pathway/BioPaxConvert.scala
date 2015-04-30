package edu.ucsc.pathway

import java.io.FileInputStream
import org.biopax.paxtools.model.level3.{Complex, Catalysis, ProteinReference, BiochemicalReaction}
import org.biopax.paxtools.pattern.constraint.ConBox

import collection.JavaConverters._

import org.biopax.paxtools.io.SimpleIOHandler
import org.biopax.paxtools.pattern.{Searcher, Pattern}


object BioPaxConvert {
  def main(args:Array[String]) = {
    val fin = new FileInputStream(args(0))
    val handler = new SimpleIOHandler()
    val model = handler.convertFromOWL(fin)

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
    })

    model.getObjects(classOf[BiochemicalReaction]).asScala.foreach( x => {
      println("==BiochemicalReaction==")
      println(x.getXref)
      println(x.getLeft)
      println(x.getRight)
      println(x.getControlledOf)
      println(x.getInteractionType)
    })






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
