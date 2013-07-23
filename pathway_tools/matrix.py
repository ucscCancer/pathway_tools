
import array
import csv

class DataException(Exception):
    pass

class NamedMatrix:
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

    def read(self, handle, delim="\t"):
        header = None
        for line in handle:
            row = line.rstrip("\n\r").split(delim)
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

    def descritize(self, bounds, labels, default=None):
        out = NamedMatrix()
        out.init_blank(cols=self.colmap, rows=self.rowmap, type='i', null_value=0)
        out.rowmap = self.rowmap
        out.colmap = self.colmap

        for i, v in enumerate(self.data):
            o = default
            if v <= bounds[0]:
                o = labels[0]
            if v > bounds[-1]:
                o = labels[-1]
            for j in range(1, len(bounds)):
                if v > bounds[j-1] and v <= bounds[j]:
                    o = labels[j]
            out.data[i] = o
        return out


