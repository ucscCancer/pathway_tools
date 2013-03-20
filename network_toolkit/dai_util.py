import dai
import network_toolkit.dai_conf

class FactorCalculator:

    def calculate(self, *args, **kw):
        raise Exception("Unimplemented method 'calculate'")


class Variable:
    def __init__(self, label, elem_type, variable, dim):
        self.label = label
        self.elem_type = elem_type
        self.dim = dim
        self.variable = variable


class VariableMap:
    def __init__(self):
        #from label/type tuple to PathwayVariable object
        self.variable_map = {}
        #from variable ID to PathwayVariable object
        self.varid_map = {}
        self.var_count = 0

    def get_or_add_variable(self, label, elem_type, dim):
        tup = (label, elem_type)
        if tup in self.variable_map:
            return self.variable_map[tup]
        var = Variable(label, elem_type, self.var_count, dim)
        self.var_count += 1
        self.insert(var)
        return var

    def insert(self, pv):
        self.variable_map[ (pv.label, pv.elem_type) ] = pv
        self.varid_map[ pv.variable ] = pv

    def get_variable_by_id(self, v):
        return self.varid_map[v]

    def get_variable_by_label(self, label, elem_type):
        try:
            return self.variable_map[ (label, elem_type) ]
        except KeyError:
            return None

    def list_variable(self, label, elem_type):
        for v_label, v_elem_type in self.variable_map:
            if label is None or v_label == label:
                if elem_type is None or v_elem_type == elem_type:            
                    yield self.variable_map[ (v_label, v_elem_type) ]


    def __iter__(self):
        for i in self.varid_map:
            yield self.varid_map[i]


class MultiDimMatrix(object):
    """
    From http://stackoverflow.com/questions/508657/multidimensional-array-in-python/508677#508677
    """
    def __init__(self, *dims):
        self._shortcuts = [i for i in self._create_shortcuts(dims)]
        self._li = [None] * (self._shortcuts.pop())
        self._shortcuts.reverse()

    def _create_shortcuts(self, dims):
        dimList = list(dims)
        dimList.reverse()
        number = 1
        yield 1
        for i in dimList:
            number *= i
            yield number

    def _flat_index(self, index):
        if len(index) != len(self._shortcuts):
            raise TypeError()

        flatIndex = 0
        for i, num in enumerate(index):
            flatIndex += num * self._shortcuts[i]
        return flatIndex

    def __getitem__(self, index):
        return self._li[self._flat_index(index)]

    def __setitem__(self, index, value):
        self._li[self._flat_index(index)] = value

    def __str__(self):
        return "%s" % (self._li)


class CPT:
    """
    CPT 

    Conditional Probablity Table

    This class stores a set of variables (and their dimensions) as well
    as the multidimensional matrix that represents the their probablities.

    Methods for normalizing the table against different variables should go here
    """
    def __init__(self, variables):
        self.variables = variables
        self.varmap = {}
        dims = []
        for i, v in enumerate(variables):
            self.varmap[v] = i
            dims.append(variables[v])
        self._table = MultiDimMatrix(*dims)

    def __str__(self):
        num_vars = len(self.variables)
        fac_states = [ 0 ] * num_vars
        return "%s" % ("\n".join(str(i) for i in self._table._li))

    def set_value(self, value, factors):
        idx = []
        for i in self.variables:
            idx.append( factors[i] )
        self._table.__setitem__(idx, value)

    def factors(self, variable_set=None):
        if variable_set is None:
            variable_set = self.variables.keys()
        num_vars = len(variable_set)
        fac_states = [ 0 ] * num_vars
        fac_index = 0
        fac_values = []

        out = []
        done = False
        while not done and fac_states[-1] < self.variables[variable_set[-1]]:
            out.append( self._table.__getitem__(fac_states) )
            inc_index = 0
            fac_states[inc_index] += 1
            while fac_states[inc_index] >= self.variables[variable_set[inc_index]]:
                fac_states[inc_index] = 0
                inc_index += 1
                if inc_index > len(fac_states)-1:
                    done = True
                    break
                fac_states[inc_index] += 1
        return out


class CPTGenerator:
    """
    CPT 

    Conditional Probablity Tabel

    This class stores a set of variables (and their dimensions) as well
    as the multidimensional matrix that represents the their probablities.

    Methods for normalizing the table against different variables should go here
    """
    def __init__(self, variables, calculator, name="CPT"):
        self.variables = variables
        self.calculator = calculator
        self.name = name

    def __str__(self):
        return self.name
        #num_vars = len(self.variables)
        #fac_states = [ 0 ] * num_vars
        #return "%s" % ("\n".join(str(i) for i in self.table._li))

    def get_variable_ids(self):
        for var in self.variables:
            yield var.variable

    def get_variable_by_id(self, id):
        for var in self.variables:
            if var.variable == id:
                return var
        return None

    def generate(self, variable_set=None):
        """
        Generate a CPT based on the variable_set (by default all the variables that where used
        to intialize the CPTGenerator), and the CPTPair rules that apply.
        """
        if variable_set is None:
            variable_set = sorted(self.get_variable_ids())

        num_vars = len(variable_set)
        fac_states = [ 0 ] * num_vars
        fac_index = 0
        fac_values = []

        var_dim_map = {}
        var_dims = []
        for i in variable_set:
            d = self.get_variable_by_id(i).dim
            var_dim_map[i] = d
            var_dims.append(d)

        cpt = CPT(var_dim_map)
        done = False
        while not done and fac_states[-1] < var_dims[-1]:
            c_map = {}
            for i, v in enumerate(variable_set):
                c_map[v] = fac_states[i]
            val = self.calculator.calculate(c_map)
            if val is None:
                raise Exception("Null value from factor calculator")
            cpt.set_value(val, c_map)
            inc_index = 0
            fac_states[inc_index] += 1
            while fac_states[inc_index] >= self.get_variable_by_id(variable_set[inc_index]).dim:
                fac_states[inc_index] = 0
                inc_index += 1
                if inc_index > len(fac_states)-1:
                    done = True
                    break
                fac_states[inc_index] += 1
        return cpt



class FactorGraph:
    """
    ExpandedPathway

    The ExpandedPathway (rename?) represents a pathway that has been expanded out to 
    a factor graph. It contains factor variable mappings, network connections, and 
    CPT definitions
    """
    def __init__(self, variable_map=None, cpt_list=None):
        self.var_map = variable_map
        self.cpt_list = cpt_list

        if self.var_map is None:
            self.var_map = VariableMap()
        if self.cpt_list is None:
            self.cpt_list = []

    def append_cpt(self, cpt):
        self.cpt_list.append(cpt)

    def describe(self):

        yield "# Variables"
        for path_var in self.var_map:
            yield "# %d\t%s\t%s" % (path_var.variable, path_var.label, path_var.elem_type)

        yield "# Factor Graphs"

        for cpt in self.cpt_list:
            node_order = []
            node_dim   = []

            for n in sorted(cpt.get_variable_ids()):
                node_order.append(n)
                node_dim.append(cpt.get_variable_by_id(n).dim)

            yield "#graph to child " + str(cpt.name)
            yield " ".join([str(i) for i in node_order])
            yield " ".join([str(i) for i in node_dim])
            yield cpt.generate()

    def generate_dai_factor_graph(self):
        vecfac = dai.VecFactor()

        elem_map = {}
        dai_var_map = {}
        factor_map = {}
        for elem in self.var_map:
            if elem.label not in elem_map:
                elem_map[ elem.label ] = { elem.elem_type : elem.variable }
            else:
                elem_map[ elem.label ][elem.elem_type] = elem.variable
            dai_var_map[elem.variable] = dai.Var(elem.variable, elem.dim)

        factor_i = 0
        for cpt_gen in self.cpt_list:
            elem = cpt_gen.generate()
            var_list = dai.VarSet()
            v_list = list(elem.variables)
            v_list.sort()
            for v in v_list:
                var_list.append(dai_var_map[v])
            factor = dai.Factor(var_list) 
            for i, v in enumerate(elem.factors()):
                factor[i] = v
            factor_label = cpt_gen.name
            name_offset = 1
            while factor_label in factor_map:
                factor_label = cpt_gen.name + "_" + str(name_offset)
                name_offset += 1
            factor_map[factor_label] = factor_i
            factor_i += 1
            vecfac.append(factor)
        return DaiFactorGraph(self.var_map, vecfac, dai_var_map, factor_map)

class DaiFactorGraph:
    def __init__(self, variable_map, vecfac, var_map, factor_map):
        self.variable_map = variable_map
        self.vecfac = vecfac
        self.var_map = var_map
        self.factor_map = factor_map

    def run_em(self, evidence):
        ev = dai.Evidence()

        #estimate = dai.CondProbEstimation()
        props = dai.PropertySet()
        #see: http://cs.ru.nl/~jorism/libDAI/doc/classdai_1_1CondProbEstimation.html#a548db9c240464901a5d5cc37fe0e8caa
        props["total_dim"] = total_dim 
        props["target_dim"] = VARIABLE_DIMENSION
        pe = dai.ParameterEstimation.construct("CondProbEstimation", props)

    def get_factor_graph(self):
        return dai.FactorGraph(self.vecfac)

    def get_inf(self, name, verbose=False):
        cls, config = network_toolkit.dai_conf.config_map[name]
        prop = dai.PropertySet()
        for c in config:
            prop[c] = config[c]
        if verbose:
            prop["verbose"] = "1"
        sn = dai.FactorGraph(self.vecfac)
        inf = cls(sn, prop)
        return inf

    def get_inf_jtree(self):
        sn = self.get_factor_graph()
        prop = dai.PropertySet()
        prop["inference"] = "SUMPROD"
        prop["updates"] = "HUGIN"
        prop["verbose"] = "1"
        inf = dai.JTree(sn, prop)
        return inf

    def get_inf_bp(self, logdomain=False, verbose=False):
        sn = dai.FactorGraph(self.vecfac)
        prop = dai.PropertySet()
        prop["tol"] = "1e-9"
        if logdomain:
            prop["logdomain"] = "1"
        else:
            prop["logdomain"] = "0"
        prop["updates"] = "SEQFIX"
        if verbose:
            prop["verbose"] = "1"
        else:
            prop["verbose"] = "0"
        lb_inf = dai.BP(sn, prop)
        return lb_inf

    def variables(self):
        for var in self.variable_map:
            yield var, self.var_map[var.variable]


