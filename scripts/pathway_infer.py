#!/usr/bin/env python

import sys
import re
import yaml
import array
import math
import csv
from StringIO import StringIO

from network_toolkit import convert, dai_util

from argparse import ArgumentParser

try:
    import dai
except ImportError:
    dai = None
    sys.stderr.write("libdai not found\n")

class Dogma:
    """
    DogmaSet 

    Holds different dogma templates that can be applied to different 
    types of pathway elements
    Right now, hard coded to dogma for protein and other
    """
    def __init__(self, handle):
        txt = handle.read()
        self.data = yaml.load(txt)

        #dogma template map, tuple of src-dst stored by edge type
        self.d_map = {}

        for dogma in self.data['dogma_template']:
            self.d_map[dogma] = {}
            for tmp in self.data['dogma_template'][dogma]:
                if len(tmp) == 3:
                    self.d_map[dogma][ tmp[2] ] = tmp[0], tmp[1]
        
        #dogma elements, set of all dogma element types
        self.el_set = {}

        for dogma in self.data['dogma_template']:
            self.el_set[dogma] = []
            for line in self.data['dogma_template'][dogma]:
                if line[0] not in self.el_set[dogma]:
                    self.el_set[dogma].append(line[0])
                if len(line) > 1:
                    if line[1] not in self.el_set[dogma]:
                        self.el_set[dogma].append(line[1])

        self.i_map = {}
        for line in self.data['interactions']:
            self.i_map[line[2]] = (line[0], line[1])

    def get_dogma_template_map(self, name):
        return self.d_map[name]
    
    def get_dogma_template_nodes(self, name):
        return self.el_set[name]

    def get_interaction_map(self):
        return self.i_map

    def get_factor_table(self, edge_type, dst_template, dst_nodetype):
        ft = None
        if dst_template in self.data['factors_fixed'][edge_type]:
            ft = self.data['factors_fixed'][edge_type][dst_nodetype]
        else:
            ft = self.data['factors_fixed'][edge_type]['*'][dst_nodetype]
        return self.data['factor_tables'][ft]


class ExpandedPathway:
    """
    ExpandedPathway

    The ExpandedPathway (rename?) represents a pathway that has been expanded out to 
    a factor graph. It contains factor variable mappings, network connections, and 
    CPT definitions
    """
    def __init__(self, factor_graph):
        self.factor_graph = factor_graph

class CPTPair:
    """
    CPTPair reprents the source information that can be used for 
    CPT generation
    """
    def __init__(self, var_a, dim_a, var_b, dim_b, ab_table):
        self.var_a = var_a
        self.var_b = var_b
        self.dim_a = dim_a
        self.dim_b = dim_b
        self.ab_table = ab_table


class PairCalculator(dai_util.FactorCalculator):
    """
    FactorCalculator

    This class represents a set of variables and CPTPair classes.
    It can be used to generate a CPT from the variables (either all of them 
    or a subset).
    """
    def __init__(self, variables, pair_table_list):
        """
        
        """

        self.variables = {}
        for v in variables:
            self.variables[ v.variable_id ] = v
        self.pair_table_list = pair_table_list
        self.pair_tables = {}
        for pt in pair_table_list:
            if pt.var_a not in self.pair_tables:
                self.pair_tables[pt.var_a] = {}
            self.pair_tables[pt.var_a][pt.var_b] = pt


    def calculate(self, v_set):
        prob = 1.0
        for i in v_set:
            if i in self.pair_tables:
                for j in v_set:
                    if j in self.pair_tables[i]:
                        prob *= self.pair_tables[i][j].ab_table[v_set[i]][v_set[j]]
        return prob

class PositiveFactorCalculator(dai_util.FactorCalculator):
    MATRIX = [
        [0.90, 0.05, 0.05],
        [0.05, 0.90, 0.05],
        [0.05, 0.05, 0.90]
    ]
    def __init__(self, src, dst):
        self.src = src
        self.dst = dst

    def calculate(self, v_set):
        return self.MATRIX[v_set[self.src]][v_set[self.dst]]


class Pathway:
    """
    Pathway

    This is the file loader for the Paradigm pathway format. It uses two data structures to store 
    the pathway data

    el_map : element map, element type stored by element name

    c_map : child_map, map of dst named to maps of src names, inner value is edge type: p_map[dst][src] = type


    """

    def __init__(self, handle):
        #element map: element type stored by element name
        self.el_map = {}
        #child_map: map of dst named to maps of src names, inner value is edge type: p_map[dst][src] = type
        self.c_map = {}
        for line in handle:
            tmp = line.rstrip().split("\t")
            if len(tmp) == 2:
                self.el_map[tmp[1]] = tmp[0]
                if tmp[1] not in self.c_map:
                    self.c_map[tmp[1]] = {}
            elif len(tmp) == 3:
                self.c_map[tmp[1]][tmp[0]] = tmp[2]
    
    def expand(self, dogma, state_count=3):
        """
        This method takes a dogma object and an interaction map object to 
        create an expanded probability graph
        """

        fg = dai_util.FactorGraph()

        #hold node enumeration value (starting at 0)

        #store edges for child reduce and emit later in code
        factor_edge = []
        for node in self.c_map:
            #what is the current dogma?
            cdog = dogma.get_dogma_template_map(self.el_map[node])
            #for ever node in dogma add a factor graph node
            for t_node in dogma.get_dogma_template_nodes( self.el_map[node] ):
                pv = fg.var_map.get_or_add_variable( label=node, elem_type=t_node, dim=state_count)
                for edge_type in cdog:
                    if cdog[edge_type][1] == t_node:
                        factor_edge.append( ( (node, cdog[edge_type][0]), (node, cdog[edge_type][1]), edge_type ) )

            #for every parent, add an edge, store the edge type (for CPD calc later)
            for parent in self.c_map[node]:
                e_type = self.c_map[node][parent]
                i_type = dogma.get_interaction_map()[e_type]

                factor_edge.append( ((parent, i_type[0]), (node, i_type[1]), e_type ) )
            
        #reduce edge graph to child mapping
        child_map = {}
        for edge in factor_edge:
            par_var = fg.var_map.get_variable_by_label( edge[0][0], edge[0][1] )
            if par_var is None:
                raise Exception("Missing node %s %s" % (edge[0][0], edge[0][1]))
            child_var = fg.var_map.get_variable_by_label( edge[1][0], edge[1][1] )
            if child_var is None:
                raise Exception("Missing node %s %s" % (edge[1][0], edge[1][1]))
            child_id = child_var.variable_id
            if child_id not in child_map:
                child_map[child_id] = [ (par_var.variable_id, edge[2]) ]
            else:
                child_map[child_id].append( (par_var.variable_id, edge[2]) )
        
        for child_id in child_map:
            #now we need to yield out CPDs
            child_variable = fg.var_map.get_variable_by_id(child_id)
            template_type = self.el_map[child_variable.variable_name]
            child_type = child_variable.variable_type

            variable_list = [i[0] for i in child_map[child_id]]
            variable_list.append(child_id)
            edge_map = {i[0] : i[1] for i in child_map[child_id]}
            #yield child_id, child_map[child_id]
            cpt_variables = []
            pair_probs = []
            for v_id in variable_list:
                cpt_variables.append( fg.var_map.get_variable_by_id(v_id) )
                if v_id != child_id:
                    edge = edge_map[v_id]
                    ft = dogma.get_factor_table(edge_type=edge, dst_template=template_type, dst_nodetype=child_type)
                    pt = CPTPair( v_id, 3, child_id, 3, ft )
                    pair_probs.append(pt)
            cpt = PairCalculator(cpt_variables, pair_probs)
            fg.append_cpt(dai_util.CPTGenerator(cpt_variables, cpt, child_variable.variable_name + ":" + child_variable.variable_type, "inputs"))
        return fg

class EvidenceMatrix:
    """
    array.array based float matrix class
    """
    def __init__(self):
        self.corner_name = "probe"
        self.data = None
        self.nrows = None
        self.ncols = None
        self.rowmap = None
        self.colmap = None

    def read(self, handle):
        header = None
        for line in handle:
            row = line.rstrip().split("\t")
            if header is None:
                header = row
                self.data = array.array("f")
                self.colmap = {}
                self.rowmap = {}
                self.ncols = len(row) - 1
                self.nrows = 0
                for i, c in enumerate(row[1:]):
                    self.colmap[c] = i
            else:
                if len(row) - 1 != self.ncols:
                    raise DataException("Misformed matrix")
                if row[0] not in self.rowmap:
                    self.rowmap[row[0]] = len(self.rowmap)
                    a = []
                    for v in row[1:]:
                        try:
                            a.append(float(v))
                        except ValueError:
                            a.append(float('Nan'))
                    self.data.extend(a)
                    self.nrows += 1

    def init_blank(self, rows, cols, type="f", null_value=float('nan')):
        self.data = array.array(type)
        self.colmap = {}
        for i,c in enumerate(cols):
            self.colmap[c] = i
        self.rowmap = {}
        for i,r in enumerate(rows):
            self.rowmap[r] = i
        self.ncols = len(cols)
        self.nrows = len(rows)
        for i in range(self.nrows):
            self.data.extend([null_value] * self.ncols)

    def has_row(self, row_name):
        return row_name in self.rowmap

    def has_col(self, col_name):
        return col_name in self.colmap

    def values(self):
        return self.data

    def get_value(self, row_name, col_name):
        return self.data[ self.rowmap[row_name] * self.ncols + self.colmap[col_name] ]

    def set_value(self, row_name, col_name, value):
        self.data[ self.rowmap[row_name] * self.ncols + self.colmap[col_name] ] = value
    
    def get_row(self, row_name):
        if row_name not in self.rowmap:
            raise KeyError
        out = {}
        for c in self.colmap:
            out[c] = self.data[ self.rowmap[row_name] * self.ncols + self.colmap[c] ]
        return out

    def get_col(self, col_name):
        if col_name not in self.colmap:
            raise KeyError
        out = {}
        for r in self.rowmap:
            out[r] = self.data[ self.rowmap[r] * self.ncols + self.colmap[col_name] ]
        return out


    def set_row(self, row_name, row_data):
        if row_name not in self.rowmap:
            raise KeyError
        row_offset = self.rowmap[row_name] * self.ncols
        for c in row_data:
            if c in self.colmap:
                self.data[ row_offset + self.colmap[c] ] = row_data[c] 

    def get_cols(self):
        if self.colmap is None:
            return None
        return self.colmap.keys()

    def get_rows(self):
        if self.rowmap is None:
            return None
        return self.rowmap.keys()
    
    def write(self, handle, missing='NA', format_str="%.5f"):
        write = csv.writer(handle, delimiter="\t", lineterminator='\n')
        col_list = self.get_cols()
        
        write.writerow([self.corner_name] + col_list)
        for rowName in self.rowmap:
            out = [rowName]
            row = self.get_row(rowName)
            for col in col_list:
                val = row[col]
                if val is None or math.isnan(val):
                    val = missing
                else:
                    val = format_str % (val)
                out.append(val)
            write.writerow(out)
    
    def set_nan(self, value=0.0):
        for i in range(len(self.data)):
            if math.isnan(self.data[i]):
                self.data[i] = value

    def descritize(self, bounds, labels):
        out = EvidenceMatrix()
        out.init_blank(cols=self.colmap, rows=self.rowmap, type='i', null_value=0)
        out.rowmap = self.rowmap
        out.colmap = self.colmap

        for i, v in enumerate(self.data):
            o = float('nan')
            if v <= bounds[0]:
                o = labels[0]
            if v > bounds[-1]:
                o = labels[-1]
            for j in range(1, len(bounds)):
                if v > bounds[j-1] and v <= bounds[j]:
                    o = labels[j]
            out.data[i] = o
        return out

                
reEvidenceArg = re.compile(r'^([^:]+):(.*)$')
if __name__ == "__main__":
    parser = ArgumentParser()    
    parser.add_argument("-d", "--dogma", 
        dest="dogma", help="DogmaFile", default=None)
    parser.add_argument("-p", "--pathway", 
        dest="pathway", help="PathwayFile", default=None)
    parser.add_argument("-e", "--evidence", nargs=2, action="append", default=[])
    parser.add_argument("--report", action="store_true", help="Print out factor graph and quit", default=False)
    parser.add_argument("-s", "--sample", help="Sample Name", default=None)

    
    args = parser.parse_args()
    
    handle = open(args.dogma)
    dogma = Dogma(handle)
    handle.close()
    
    handle = open(args.pathway)
    pathway = Pathway(handle)
    handle.close()  
    
    e_map = {}
    for e_type, e_path in args.evidence:
        mat = EvidenceMatrix()
        handle = open(e_path)
        mat.read(handle)
        handle.close()
        e_map[e_type] = mat.descritize([0.33, 0.66], [0,1,2])

    if args.report:
        expanded_pathway = pathway.expand(dogma)
        for line in expanded_pathway.describe():
            print line
    else:
        if not dai:
            sys.stderr.write("Cannot proceed without libdai")
            sys.exit(1)
        else:
            vecfac = dai.VecFactor()

            expanded_pathway = pathway.expand(dogma)

            clamp_map = {}
            if args.sample is not None:
                for e_type in e_map:
                    col = e_map[e_type].get_col(args.sample)
                    for probe in col:
                        v = expanded_pathway.var_map.get_variable_by_label(probe, e_type)
                        if v is not None:
                            v_obs = expanded_pathway.var_map.get_or_add_variable(probe, e_type + ":obs", 3)
                            clamp_map[v_obs.variable_id] = col[probe]
                            expanded_pathway.append_cpt( dai_util.CPTGenerator([v, v_obs], PositiveFactorCalculator(v_obs.variable_id, v.variable_id), v_obs.variable_name, "obs_connect" ) )

            dai_fg = expanded_pathway.generate_dai_factor_graph()
            lb_inf = dai_fg.get_inf_bp(verbose=True)

            lb_inf.init()
            lb_inf.run()

            lb_inf_post = lb_inf.clone()

            for v in clamp_map:
                lb_inf_post.clamp(v, clamp_map[v])

            lb_inf_post.init()
            lb_inf_post.run()
                        
            for variable, dai_variable in dai_fg.variables():
                factor_pre = lb_inf.belief( dai_variable )
                factor_post = lb_inf_post.belief( dai_variable )
                for i in range(factor_pre.nrStates()):
                    print variable.label, variable.elem_type, i, factor_pre[i], factor_post[i]
            print lb_inf_post._iters
            



