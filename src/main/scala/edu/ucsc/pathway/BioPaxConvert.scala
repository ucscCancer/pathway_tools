package edu.ucsc.pathway

import java.io.{File, FileInputStream}
import java.util.zip.GZIPInputStream
import java.{lang, util}
import javax.xml.parsers.DocumentBuilderFactory
import javax.xml.transform.{OutputKeys, TransformerFactory}
import javax.xml.transform.dom.DOMSource
import javax.xml.transform.stream.StreamResult
import org.biopax.paxtools.model.BioPAXElement
import org.biopax.paxtools.model.level3._
import org.biopax.paxtools.pattern.constraint.ConBox
import org.w3c.dom.Document

import collection.JavaConverters._

import org.biopax.paxtools.io.SimpleIOHandler
import org.biopax.paxtools.pattern.{Searcher, Pattern}

class Vertex(val src: Entity, val attr : Array[(String,String)]) {
  override def toString() : String = {
    "%s(%s)".format(src.getRDFId, attr.mkString(","))
  }

  def getAttributes() : Array[(String,String)] = {
    var out = attr ++ Array( ("url", src.getRDFId) )
    if (src.getDisplayName != null)
      out = out ++ Array(("displayName", src.getDisplayName))
    out = out ++ src.getName.asScala.map( x => ("name", x) )
    if (src.getStandardName != null)
      out = out ++ Array(("standardName", src.getStandardName))
    out
  }
}

class Edge(val src: Entity, val edge_type : String, val src_id : String, val dst_id : String, val attr : Array[(String,String)] ) {
  override def toString() : String = {
    "%s : %s %s %s".format(src.getRDFId, src_id, edge_type, dst_id)
  }

  def getAttributes() : Array[(String,String)] = {
    var out = attr ++ Array( ("url", src.getRDFId) )
    if (src.getDisplayName != null)
      out = out ++ Array(("displayName", src.getDisplayName))
    out = out ++ src.getName.asScala.map( x => ("name", x) )
    if (src.getStandardName != null)
      out = out ++ Array(("standardName", src.getStandardName))
    out
  }

}

trait Subnet {
  def vertices() : Array[Vertex]
  def edges() : Array[Edge]
  def inputVertices() : Array[Vertex]
  def outputVertices() : Array[Vertex]
}

abstract class BaseSubnet(stack: Set[String]) extends Subnet {

}

class BlankSubnet extends Subnet {
  override def vertices(): Array[Vertex] = Array()

  override def outputVertices(): Array[Vertex] = Array()

  override def edges(): Array[Edge] = Array()

  override def inputVertices(): Array[Vertex] = Array()
}


object BioChemicalReactionSubnet {
  val interaction_map = Map(
    "-a>" -> "-a>",
    "ACTIVATION" -> "-a>",
    "INHIBITION" -> "-a|",
    "STATE_TRANSITION" -> "-a>",
    "DISSOCIATION" -> "-a>",
    "TRUNCATION" -> "-a|",
    "HETERODIMER_ASSOCIATION" -> "component>"  //BUG: this needs to be check, may end up pointing to other member, not the complex
  )
}

abstract class ReactionSubnet[T <: Conversion](stack:Set[String], entity:T) extends BaseSubnet(stack) {
  val left:Array[Subnet] = entity.getLeft.asScala.toArray.map( x => new PhysicalEntitySubnet(stack + entity.getRDFId, x) )
  val right:Array[Subnet] = entity.getRight.asScala.toArray.map( x => new PhysicalEntitySubnet(stack + entity.getRDFId, x) )
  val controls:Array[Subnet] = entity.getControlledOf.asScala.toArray.map( x => new ControlSubnet(stack + entity.getRDFId, x) )

  val reactionEdgeType : String
  val controlEdgeType : String

  def vertices() : Array[Vertex] = {
    left.flatMap( x => x.vertices() ) ++ right.flatMap( x => x.vertices() ) ++ controls.flatMap( x => x.vertices() )
  }

  def inputVertices() : Array[Vertex] = left.flatMap(_.vertices()) ++ controls.flatMap(_.vertices())
  def outputVertices() : Array[Vertex] = right.flatMap(_.vertices())

  def edges() : Array[Edge] = {
    val leftVerts = left.flatMap( _.outputVertices() )
    val rightVerts = right.flatMap( _.inputVertices() )

    val edges_control = controls.flatMap( x => x.outputVertices() ).flatMap( x => rightVerts.map( y => new Edge(entity, controlEdgeType, x.src.getRDFId, y.src.getRDFId, Array()) ))
    val edge_reaction = leftVerts.flatMap( x => rightVerts.map( y => new Edge(entity, reactionEdgeType, x.src.getRDFId, y.src.getRDFId, Array()) ) )

    left.flatMap( _.edges() ) ++ right.flatMap( _.edges()) ++ edges_control ++ edge_reaction
  }
}

class ConversionSubnet(stack:Set[String], entity: Conversion) extends BaseSubnet(stack) {

  val sub = if (entity.getModelInterface == classOf[BiochemicalReaction]) {
    new BioChemicalReactionSubnet(stack + entity.getRDFId, entity.asInstanceOf[BiochemicalReaction])
  } else if (entity.getModelInterface == classOf[ComplexAssembly]) {
    new ComplexAssemblySubnet(stack + entity.getRDFId, entity.asInstanceOf[ComplexAssembly])
  } else if (entity.getModelInterface == classOf[TransportWithBiochemicalReaction]) {
    new TransportWithBiochemicalReactionSubnet(stack, entity.asInstanceOf[TransportWithBiochemicalReaction])
  } else if (entity.getModelInterface == classOf[Transport]) {
    new TransportSubnet(entity.asInstanceOf[Transport])
  } else if (entity.getModelInterface == classOf[Degradation]) {
    new DegradationSubnet(stack, entity.asInstanceOf[Degradation])
  } else if (entity.getModelInterface == classOf[Conversion]) {
    //labeling something as a 'conversion' is the most generic way to refer to a reaction
    new ReactionSubnet[Conversion](stack, entity) {
      val reactionEdgeType = "-converts>"
      val controlEdgeType = "-a>"
    }
  } else {
    throw new RuntimeException("Unknown Conversion Type: " + entity.getModelInterface.toString)
  }

  override def vertices(): Array[Vertex] = sub.vertices()

  override def outputVertices(): Array[Vertex] = sub.outputVertices()

  override def edges(): Array[Edge] = sub.edges()

  override def inputVertices(): Array[Vertex] = sub.inputVertices()
}


class InteractionSubnet(stack:Set[String], entity: Interaction) extends BaseSubnet(stack) {

  val sub = if (entity.getModelInterface == classOf[TemplateReaction]) {
    new TemplateReactionSubnet(stack, entity.asInstanceOf[TemplateReaction])
  } else {
    throw new RuntimeException("Unknown Conversion Type")
  }

  override def vertices(): Array[Vertex] = sub.vertices()

  override def outputVertices(): Array[Vertex] = sub.outputVertices()

  override def edges(): Array[Edge] = sub.edges()

  override def inputVertices(): Array[Vertex] = sub.inputVertices()
}




class BioChemicalReactionSubnet(stack:Set[String], entity:BiochemicalReaction) extends BaseSubnet(stack) {

  val left:Array[Subnet] = entity.getLeft.asScala.toArray.map( x => new PhysicalEntitySubnet(stack + entity.getRDFId, x) )
  val right:Array[Subnet] = entity.getRight.asScala.toArray.map( x => new PhysicalEntitySubnet(stack + entity.getRDFId, x) )
  val controls:Array[Subnet] = entity.getControlledOf.asScala.toArray.map( x => new ControlSubnet(stack + entity.getRDFId, x) )

  def vertices() : Array[Vertex] = {
    left.flatMap( x => x.vertices() ) ++ right.flatMap( x => x.vertices() ) ++ controls.flatMap( x => x.vertices() )
  }

  def inputVertices() : Array[Vertex] = left.flatMap(_.vertices()) ++ controls.flatMap(_.vertices())
  def outputVertices() : Array[Vertex] = right.flatMap(_.vertices())

  def edges() : Array[Edge] = {

    val iset = entity.getInteractionType.asScala.flatMap( x => x.getTerm.asScala.map( y => {
      BioChemicalReactionSubnet.interaction_map(y)
    } ) )

    if (iset.size == 0) {
      iset.add("-a>")
    }

    val leftVerts = left.flatMap( _.outputVertices() )
    val rightVerts = right.flatMap( _.inputVertices() )

    val edges_control = controls.flatMap( x => x.outputVertices() ).flatMap( x => rightVerts.map( y => new Edge(entity, "-control>", x.src.getRDFId, y.src.getRDFId, Array()) ))

    val edge_reaction = iset.flatMap( itype => leftVerts.flatMap( x => rightVerts.map( y => new Edge(entity, itype, x.src.getRDFId, y.src.getRDFId, Array()) ) ) )
    left.flatMap( _.edges() ) ++ right.flatMap( _.edges()) ++ edges_control ++ edge_reaction
  }
}



class ModulationSubnet(stack:Set[String], entity: Modulation) extends BlankSubnet {
  //FIXME: implement this section
}

class DegradationSubnet(stack:Set[String], entity:Degradation) extends BaseSubnet(stack) {

  val left:Array[Subnet] = entity.getLeft.asScala.toArray.map( x => new PhysicalEntitySubnet(stack + entity.getRDFId, x) )
  val right:Array[Subnet] = entity.getRight.asScala.toArray.map( x => new PhysicalEntitySubnet(stack + entity.getRDFId, x) )
  val controls:Array[Subnet] = entity.getControlledOf.asScala.toArray.map( x => new ControlSubnet(stack + entity.getRDFId, x) )

  def vertices() : Array[Vertex] = {
    left.flatMap( x => x.vertices() ) ++ right.flatMap( x => x.vertices() ) ++ controls.flatMap( x => x.vertices() )
  }

  def inputVertices() : Array[Vertex] = left.flatMap(_.vertices()) ++ controls.flatMap(_.vertices())
  def outputVertices() : Array[Vertex] = right.flatMap(_.vertices())

  def edges() : Array[Edge] = {
    val leftVerts = left.flatMap( _.outputVertices() )
    val rightVerts = right.flatMap( _.inputVertices() )

    val edges_control = controls.flatMap( x => x.outputVertices() ).flatMap( x => rightVerts.map( y => new Edge(entity, "-control|", x.src.getRDFId, y.src.getRDFId, Array()) ))
    val edge_reaction = leftVerts.flatMap( x => rightVerts.map( y => new Edge(entity, "-a|", x.src.getRDFId, y.src.getRDFId, Array()) ) )

    left.flatMap( _.edges() ) ++ right.flatMap( _.edges()) ++ edges_control ++ edge_reaction
  }
}

class TransportWithBiochemicalReactionSubnet(stack:Set[String], entity:TransportWithBiochemicalReaction) extends ReactionSubnet(stack, entity) {
  override val reactionEdgeType: String = "-transport>"
  override val controlEdgeType: String = "-control>"
}

class ComplexAssemblySubnet(stack:Set[String], entity:ComplexAssembly) extends ReactionSubnet[ComplexAssembly](stack, entity) {
  val reactionEdgeType: String = "component>"
  val controlEdgeType: String = "-a>"
}

class ControlSubnet(stack:Set[String], control: Control) extends BaseSubnet(stack) {

  val sub_controller : Array[Subnet] = control.getController.asScala.toArray.map( x => {
    if ( x.isInstanceOf[PhysicalEntity] ) {
      new PhysicalEntitySubnet(stack + control.getRDFId, x.asInstanceOf[PhysicalEntity])
    } else if (x.getModelInterface == classOf[Pathway]) {
      new BlankSubnet() //FIXME: Should be abstract node here
    } else {
      throw new RuntimeException("Unknown controller: " + x.getModelInterface.toString)
    }
  })

  val sub_controlled =
    control.getControlled.asScala.toArray.filter( x => !stack.contains(x.getRDFId) ).map( x => {
      if ( x.isInstanceOf[PhysicalEntity] ) {
        new PhysicalEntitySubnet(stack + control.getRDFId, x.asInstanceOf[PhysicalEntity])
      } else if (x.isInstanceOf[Conversion]) {
        new ConversionSubnet(stack, x.asInstanceOf[Conversion])
      } else if (x.isInstanceOf[Interaction]) {
        new InteractionSubnet(stack + control.getRDFId, x.asInstanceOf[Interaction])
      } else if (x.getModelInterface == classOf[Pathway]) {
        new BlankSubnet() //FIXME: Should be abstract node here
      } else {
        throw new RuntimeException("Unknown controlled: " + x.getModelInterface.toString)
      }
    } )

  def inputVertices() : Array[Vertex] = sub_controller.flatMap(_.vertices())
  def outputVertices() : Array[Vertex] = sub_controlled.flatMap(_.outputVertices())

  def vertices() : Array[Vertex] = {
    inputVertices() ++ outputVertices()
  }

  def edges() : Array[Edge] = {
    inputVertices().flatMap( x => outputVertices().map(y => new Edge(control, "-a>", x.src.getRDFId, y.src.getRDFId, Array())))
  }
}

class TransportSubnet(entity: Transport) extends BlankSubnet {

}

class CatalysisSubnet(stack:Set[String], entity: Catalysis) extends BaseSubnet(stack) {

  val cofactor_subnet = entity.getCofactor.asScala.map( x => new PhysicalEntitySubnet(stack + entity.getRDFId, x) ).toArray

  override def vertices(): Array[Vertex] = {
    inputVertices() ++ outputVertices()
  }

  override def outputVertices(): Array[Vertex] = Array()

  override def edges(): Array[Edge] = Array()

  override def inputVertices(): Array[Vertex] = {
    cofactor_subnet.flatMap(_.vertices())
  }
}

class TemplateReactionSubnet(stack:Set[String], entity: TemplateReaction) extends BaseSubnet(stack) {

  val product_subnet = entity.getProduct.asScala.toArray.map( x => new PhysicalEntitySubnet(stack + entity.getRDFId, x) )
  val control_subnet =
    entity.getControlledOf.asScala.toArray
      .filter( x => !stack.contains(x.getRDFId))
      .map( x => new ControlSubnet(stack + entity.getRDFId, x) )


  override def vertices(): Array[Vertex] = {
    inputVertices() ++ outputVertices()
  }
  override def outputVertices(): Array[Vertex] = {
    product_subnet.flatMap(_.vertices())
  }

  override def edges(): Array[Edge] = {
    inputVertices().flatMap( x => outputVertices().map( y => new Edge(entity, "->", x.src.getRDFId, y.src.getRDFId, Array()) ) )
  }

  override def inputVertices(): Array[Vertex] = {
    control_subnet.flatMap(_.vertices())
  }
}

class TemplateReactionRegulationSubnet(entity: TemplateReactionRegulation) extends Subnet {
  override def vertices(): Array[Vertex] = {
    Array()
  }

  override def outputVertices(): Array[Vertex] = Array()

  override def edges(): Array[Edge] = Array()

  override def inputVertices(): Array[Vertex] = Array()

}

class ProteinSubnet(entity: Protein) extends Subnet {

  def inputVertices() : Array[Vertex] = vertices()
  def outputVertices() : Array[Vertex] = vertices()

  def vertices() :Array[Vertex] = {
    val xrefs = entity.getXref.asScala.map( x => ("xref", "%s:%s".format(x.getDb, x.getId)) )
    val attr = if (entity.getEntityReference == null)
      Array("type" -> "protein") ++ xrefs
    else
      Array("type" -> "protein") ++ xrefs ++ entity.getEntityReference.getXref.asScala.map( x => ("xref", "%s:%s".format(x.getDb, x.getId)) )

    return Array(new Vertex(entity, attr))
  }

  def edges() : Array[Edge] = {
    Array()
  }
}

class RnaSubnet(entity: Rna) extends Subnet {
  def inputVertices() : Array[Vertex] = vertices()
  def outputVertices() : Array[Vertex] = vertices()

  def vertices() :Array[Vertex] = {
    val xrefs = entity.getXref.asScala.map( x => ("xref", "%s:%s".format(x.getDb, x.getId)) )
    val attr = if (entity.getEntityReference == null)
      Array("type" -> "rna") ++ xrefs
    else
      Array("type" -> "rna") ++ xrefs ++ entity.getEntityReference.getXref.asScala.map( x => ("xref", "%s:%s".format(x.getDb, x.getId)) )

    return Array(new Vertex(entity, attr))
  }

  def edges() : Array[Edge] = {
    Array()
  }
}

class SmallMoleculeSubnet(entity: SmallMolecule) extends Subnet {

  def inputVertices() : Array[Vertex] = vertices()
  def outputVertices() : Array[Vertex] = vertices()

  def vertices() :Array[Vertex] = {
    return Array(new Vertex(entity, Array("type" -> "small molecule")))
  }

  def edges() : Array[Edge] = {
    Array()
  }
}

class ComplexSubnet(stack:Set[String], entity: Complex) extends BaseSubnet(stack) {
  val components = entity.getComponent.asScala.toArray.map( x => new PhysicalEntitySubnet(stack + entity.getRDFId, x) )

  def inputVertices() : Array[Vertex] = components.flatMap( x => x.vertices() )
  def outputVertices() : Array[Vertex] = Array( new Vertex(entity, Array( "type" -> "complex")) )

  def vertices() : Array[Vertex] = {
    inputVertices() ++ outputVertices()
  }

  def edges() : Array[Edge] = {
    components.flatMap( x => x.outputVertices() ).map( x => {
      new Edge(entity, "component>", x.src.getRDFId, entity.getRDFId, Array() )
    }) ++ components.flatMap( x => x.edges() )
  }
}

class PhysicalEntitySubnet(stack:Set[String], entity:PhysicalEntity) extends BaseSubnet(stack) {

  val sub : Subnet = if (entity.getModelInterface == classOf[Protein]) {
    new ProteinSubnet(entity.asInstanceOf[Protein])
  } else if (entity.getModelInterface == classOf[Complex]) {
    new ComplexSubnet(stack + entity.getRDFId, entity.asInstanceOf[Complex])
  } else if (entity.getModelInterface == classOf[SmallMolecule]) {
    new SmallMoleculeSubnet(entity.asInstanceOf[SmallMolecule])
  } else if (entity.getModelInterface == classOf[PhysicalEntity]) {
    new BlankSubnet
  } else if (entity.getModelInterface == classOf[DnaRegion]) {
    new BlankSubnet
  } else if (entity.getModelInterface == classOf[RnaRegion]) {
    new BlankSubnet
  } else if (entity.getModelInterface == classOf[Rna]) {
    new RnaSubnet(entity.asInstanceOf[Rna])
  } else if (entity.getModelInterface == classOf[Dna]) {
    new BlankSubnet
  } else if (entity.getModelInterface == classOf[Interaction]) {
    new BlankSubnet
  } else {
    throw new RuntimeException("Unknown Component: " + entity.getModelInterface.toString)
  }

  def inputVertices() : Array[Vertex] = sub.inputVertices()
  def outputVertices() : Array[Vertex] = sub.outputVertices()

  def vertices() : Array[Vertex] = {
    sub.vertices()
  }

  def edges() : Array[Edge] = {
    sub.edges()
  }
}


class PathwaySubnet(stack:Set[String], entity: Pathway) extends BaseSubnet(stack) {

  val sub : Array[Subnet] = entity.getPathwayComponent.asScala.toArray.map( y => {
    if (y.getModelInterface == classOf[BiochemicalReaction] ) {
      new BioChemicalReactionSubnet(stack + entity.getRDFId, y.asInstanceOf[BiochemicalReaction])
    } else if (y.getModelInterface == classOf[Degradation]) {
      val deg = y.asInstanceOf[Degradation]
      new DegradationSubnet(stack + entity.getRDFId, deg)
    } else if (y.getModelInterface == classOf[TemplateReaction]) {
      new TemplateReactionSubnet(stack + entity.getRDFId, y.asInstanceOf[TemplateReaction])
    } else if (y.getModelInterface == classOf[TemplateReactionRegulation]) {
      new TemplateReactionRegulationSubnet(y.asInstanceOf[TemplateReactionRegulation])
    } else if (y.getModelInterface == classOf[Control]) {
      new ControlSubnet(stack + entity.getRDFId, y.asInstanceOf[Control])
    } else if (y.getModelInterface == classOf[ComplexAssembly]) {
      new ComplexAssemblySubnet(stack + entity.getRDFId, y.asInstanceOf[ComplexAssembly])
    } else if (y.getModelInterface == classOf[Pathway]) {
      new BlankSubnet()
    } else if (y.getModelInterface == classOf[Catalysis]) {
      new CatalysisSubnet(stack + entity.getRDFId, y.asInstanceOf[Catalysis])
    } else if (y.getModelInterface == classOf[TransportWithBiochemicalReaction]) {
      new TransportWithBiochemicalReactionSubnet(stack + entity.getRDFId, y.asInstanceOf[TransportWithBiochemicalReaction])
    } else if (y.getModelInterface == classOf[Transport]) {
      new TransportSubnet(y.asInstanceOf[Transport])
    } else if (y.getModelInterface == classOf[Modulation]) {
      new ModulationSubnet(stack + entity.getRDFId, y.asInstanceOf[Modulation])
    } else if (y.getModelInterface == classOf[Conversion]) {
      new ConversionSubnet(stack + entity.getRDFId, y.asInstanceOf[Conversion])
    } else if (y.getModelInterface == classOf[Interaction]) {
      new BlankSubnet
    } else {
      throw new RuntimeException("Unknown Type: " + y.getModelInterface.toString)
    }
  })

  def inputVertices() : Array[Vertex] = sub.flatMap( _.inputVertices() )
  def outputVertices() : Array[Vertex] = sub.flatMap( _.outputVertices() )

  def vertices() : Array[Vertex] = {
    sub.flatMap( x => x.vertices() )
  }

  def edges() : Array[Edge] = {
    sub.flatMap( x => x.edges() )
  }
}

object BioPaxConvert {

  def stringClean(in:String) : String = {
    in.replace(" ", "_").replace("/", "_").replace("(", "_").replace(")", "_").replace("$", "_")
  }

  def pathway_name(pathway:Pathway) : String = {
    if (pathway.getStandardName != null)
      stringClean(pathway.getStandardName)
    else
      stringClean(pathway.getName.asScala.last)
  }

  def main(args:Array[String]) = {
    val fin = if (args(0).endsWith(".gz")) {
      new GZIPInputStream(new FileInputStream(args(0)))
    } else {
      new FileInputStream(args(0))
    }
    val outdir = args(1)
    val handler = new SimpleIOHandler()
    val model = handler.convertFromOWL(fin)

    //find all the pathways
    val pathway_list = model.getObjects(classOf[Pathway]).asScala.map( x => {
      (pathway_name(x),x)
    }).filter( x => x._2.getPathwayComponent().asScala.seq.size > 0 ) //filter out the pathways with no components

    //give each of them a unique name
    var pathway_map = pathway_list.groupBy(_._1).flatMap( x => {
      if (x._2.size == 1)
        x._2.toSeq
      else
        x._2.zipWithIndex.map(y => ("%s_%d".format(y._1._1, y._2), y._1._2) )
    }).toMap

    pathway_map.foreach( x => {
      println("Processing", x._1)
    })

    val subnets : Map[String, Subnet] = pathway_map.map( x => {
      (x._1, new PathwaySubnet(Set[String](), x._2))
    }).filter( x => {
      x._2.vertices().length > 0
    })

    subnets.foreach( x => {
      println(x._1, x._2.vertices().length)
    })

    subnets.foreach( x => {
      val path = new File(outdir, x._1 + ".xgmml")
      write_xgmml(path.getAbsolutePath,
        x._2.vertices().groupBy(x => x.src.getRDFId).map( x => x._2.head ).toArray,
        x._2.edges().groupBy(x => x.src.getRDFId).map( x => x._2.head).toArray
      )
    })


  }


  def write_xgmml(path: String, vertices : Array[Vertex], edges : Array[Edge]) = {
    val docFactory = DocumentBuilderFactory.newInstance();
    val docBuilder = docFactory.newDocumentBuilder();

    val idmap = vertices.zipWithIndex.map( x => (x._1.src.getRDFId, x._2) ).toMap

    // root elements
    val doc = docBuilder.newDocument();

    val graph_node = doc.createElement("graph")
    graph_node.setAttribute("xmlns", "http://www.cs.rpi.edu/XGMML")

    graph_node.setAttribute("xmlns:dc", "http://purl.org/dc/elements/1.1/")
    graph_node.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink" )
    graph_node.setAttribute("xmlns:rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#" )
    graph_node.setAttribute("xmlns:cy", "http://www.cytoscape.org" )
    graph_node.setAttribute("directed", "1")
    doc.appendChild(graph_node)

    vertices.foreach( x => {
      val node = doc.createElement("node")

      node.setAttribute("id", idmap(x.src.getRDFId).toString)
      x.getAttributes().groupBy(_._1).foreach( y => {
        val o = doc.createElement("att")
        if (y._2.size == 1) {
          o.setAttribute("type", "string")
          o.setAttribute("name", y._1)
          o.setAttribute("value", y._2(0)._2)
        } else {
          o.setAttribute("type", "list")
          o.setAttribute("name", y._1)
          y._2.foreach( z => {
            val e = doc.createElement("att")
            e.setAttribute("type", "string")
            e.setAttribute("name", y._1)
            e.setAttribute("value", z._2)
            o.appendChild(e)
          })
        }
        node.appendChild(o)
      } )
      graph_node.appendChild(node)
    })

    edges.foreach( x => {
      val node = doc.createElement("edge")
      node.setAttribute("source", idmap(x.src_id).toString)
      node.setAttribute("target", idmap(x.dst_id).toString)
      val o = doc.createElement("att")
      o.setAttribute("type", "string")
      o.setAttribute("name", "type")
      o.setAttribute("value", x.edge_type)
      node.appendChild(o)

      x.getAttributes().groupBy(_._1).foreach( y => {
        val o = doc.createElement("att")
        if (y._2.size == 1) {
          o.setAttribute("type", "string")
          o.setAttribute("name", y._1)
          o.setAttribute("value", y._2(0)._2)
        } else {
          o.setAttribute("type", "list")
          o.setAttribute("name", y._1)
          y._2.foreach( z => {
            val e = doc.createElement("att")
            e.setAttribute("type", "string")
            e.setAttribute("name", y._1)
            e.setAttribute("value", z._2)
            o.appendChild(e)
          })
        }
        node.appendChild(o)
      } )
      graph_node.appendChild(node)
    })


    // write the content into xml file
    val transformerFactory = TransformerFactory.newInstance();
    val transformer = transformerFactory.newTransformer();
    transformer.setOutputProperty(OutputKeys.INDENT, "yes");
    transformer.setOutputProperty(OutputKeys.METHOD, "xml");
    transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount","4");

    val source = new DOMSource(doc);

    val result = new StreamResult(new File(path));
    // Output to console for testing
    // StreamResult result = new StreamResult(System.out);
    transformer.transform(source, result)
  }

}

