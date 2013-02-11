## Data module
## Written By: Sam Ng
## Last Updated: 9/28/11
import re, sys, urllib2, os, random, math
from copy import deepcopy
import mCalculate

def log(msg, die = False):
    """logger function"""
    sys.stderr.write(msg)
    if die:
        sys.exit(1)

def openAnyFile(inf):
    """performs an open() on a file or url"""
    if inf.startswith("http"):
        stream = urllib2.urlopen(inf)
    else:
        stream = open(inf, 'r')
    return stream

def retColumns(inf, delim = "\t"):
    """returns the columns of a .tsv"""
    f = openAnyFile(inf)
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank header\n", die = True)
    line = line.rstrip("\r\n")
    return(re.split(delim, line)[1:])

def retRows(inf, delim = "\t"):
    """returns the rows of a .tsv"""
    rows = []
    f = openAnyFile(inf)
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank header\n", die = True)
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        rows.append(re.split(delim, line)[0])
    return(rows)
    
def rCRSData(inf, appendData = dict(), delim = "\t", retFeatures = False, debug = False):
    """reads .tsv into a [col][row] dictionary"""
    inData = dict()
    colFeatures = []
    rowFeatures = []
    ## copy appendData
    for col in appendData.keys():
        inData[col] = dict()
        for row in appendData[col].keys():
            inData[col][row] = [appendData[col][row]]
    ## read header
    f = openAnyFile(inf)
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank on line 1\n", die = True)
    line = line.rstrip("\r\n")
    pline = re.split(delim, line)
    if debug:
        log("%s\nLENGTH: %s\n" % (line, len(pline)))
    colFeatures = pline[1:]
    for col in colFeatures:
        if col not in inData:
            inData[col] = dict()
    ## read data
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        pline = re.split(delim, line)
        rowFeatures.append(pline[0])
        if debug:
            log("%s\nLENGTH: %s\n" % (line, len(pline)))
        if len(pline) != (1+len(colFeatures)):
            log("ERROR: length of line does not match the rest of the file\n", die = True)
        for i, col in enumerate(colFeatures):
            row = pline[0]
            if row not in inData[col]:
                inData[col][row] = []
            if pline[i+1] == "":
                inData[col][row].append("NA")
            else:            
                inData[col][row].append(pline[i+1])
    f.close()
    ## average entries
    for col in inData.keys():
        for row in inData[col].keys():
            inData[col][row] = mean(inData[col][row], null = inData[col][row][0])
    if retFeatures:
        return(inData, colFeatures, rowFeatures)
    else:
        return(inData)

def wCRSData(outf, outData, delim = "\t", useCols = None, useRows = None, null = "NA"):
    """write [col][row] dictionary to .tsv"""
    ## get colFeatures and rowFeatures
    if useCols == None:
        colFeatures = outData.keys()
    else:
        colFeatures = useCols
    if useRows == None:
        rowFeatures = []
        for col in colFeatures:
            if col in outData:
                rowFeatures = outData[col].keys()
                break
    else:
        rowFeatures = useRows
    ## write header
    f = open(outf, "w")
    f.write("id")
    for col in colFeatures:
        f.write("\t%s" % (col))
    f.write("\n")
    for row in rowFeatures:
        f.write("%s" % (row))
        for col in colFeatures:
            if col in outData:
                if row in outData[col]:
                    f.write("\t%s" % (outData[col][row]))
                else:
                    f.write("\t%s" % (null)) 
            else:
                f.write("\t%s" % (null))
        f.write("\n")
    f.close()

def rwCRSData(outf, inf, delim = "\t", useCols = None, useRows = None, null = "NA", numeric = False, enumerateRows = False):
    """read and write .tsv for lower memory usage and efficiency"""
    f = openAnyFile(inf)
    o = open(outf, "w")
    ## read header
    line = f.readline()
    if line.isspace():
        log("ERROR: no header found\n", die = True)
    line = line.rstrip("\r\n")
    dataCols = re.split(delim, line)[1:]
    ## write header
    if useCols is None:
        useCols = set(dataCols)
    outCols = list(set(useCols)&set(dataCols))+list(set(useCols)-set(dataCols))
    outCols.sort()
    o.write("id")
    for i in outCols:
        o.write("\t%s" % (i))
    o.write("\n")
    ## read and write rest of the file
    if enumerateRows:
        rowID = 0
        e = open(outf+".fmap", "w")
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        pline = re.split(delim, line)
        if useRows is not None:
            if pline[0] not in useRows:
                continue
        if enumerateRows:
            rowID += 1
            rowItem = "rid_%s" % (rowID)
            e.write("%s\t%s\n" % (pline[0], rowItem))
        else:
            rowItem = pline[0]
        ## read row
        lData = dict()
        for i in range(len(dataCols)):
            if numeric:
                try:
                    fval = float(pline[i+1])
                    lData[dataCols[i]] = "%.10f" % (fval)
                except ValueError:
                    lData[dataCols[i]] = null
            else:
                lData[dataCols[i]] = pline[i+1]
        ## write row
        o.write("%s" % (rowItem))
        for i in outCols:
            if i in lData:
                if lData[i] == "":
                    o.write("\t%s" % (null))
                else:
                    o.write("\t%s" % (lData[i]))
            else:
                o.write("\t%s" % (null))
        o.write("\n")
    if enumerateRows:
        e.close()
    f.close()
    o.close()

def remapCRSData(inData, colMap = dict(), rowMap = dict()):
    """remaps the keys of a CRS data matrix"""
    ## map entries
    outData = {}
    for col in inData.keys():
        mcol = col
        if col in colMap:
            mcol = colMap[col]
        if mcol not in outData:
            outData[mcol] = {}
        for row in inData[col].keys():
            mrow = row
            if row in rowMap:
                mrow = rowMap[row]
            if mrow not in outData[mcol]:
                outData[mcol][mrow] = [] 
            outData[mcol][mrow].append(inData[col][row])
    ## average entries
    for col in outData.keys():
        for row in outData[col].keys():
            outData[col][row] = mean(outData[col][row], null = outData[col][row][0])
    return(outData)

def rPARADIGM(inf, delim = "\t"):
    """read PARADIGM format .fa output"""
    inLikelihood = dict()
    inScore = dict()
    f = openAnyFile(inf)
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            pline = re.split("[= ]", line)
            sample = pline[1]
            inLikelihood[sample] = float(pline[3])
            inScore[sample] = dict()
        else:
            pline = re.split(delim, line)
            feature = pline[0]    
            inScore[sample][feature] = float(pline[1])
    f.close()
    return(inLikelihood, inScore)

def r2Col(inf, appendData = dict(), delim = "\t", header = False, null = ""):
    """read 2 column data"""
    inData = deepcopy(appendData)
    f = openAnyFile(inf)
    if header:
        line = f.readline()
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        pline = re.split(delim, line)
        if len(pline[1]) == 0:
            pline[1] = null
        if len(pline) != 2:
            log("ERROR: Length of data line is not 2\n", die = True)
        inData[pline[0]] = pline[1]
    f.close()
    return(inData)

def rSet(inf, header = True, delim = "\t", enumerate = False):
    """read sets file"""
    inSets = dict()                 #Dictionary with (name : set)
    f = openAnyFile(inf)
    if header:
        f.readline()
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\t\r\n")
        pline = re.split(delim, line)
        if enumerate:
            value = 1
            while "_".join([pline[0]]+[value]) in inSets:
                value += 1
            inSets["_".join([pline[0]]+[value])] = set(pline[1:])
        else:
            inSets[pline[0]] = set(pline[1:])
    f.close()
    return(inSets)

def rList(inf):
    """read 1 column list"""
    inList = list()
    f = openAnyFile(inf)
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\t\r\n")
        inList.append(line)
    f.close()
    return(inList)

def lineCount(inf):
    """returns line count"""
    f = openAnyFile(inf)
    for i, l in enumerate(f):
        pass
    f.close()
    return(i+1)

def rMeta(clinf, delim = "\t", null = True):
    """read .meta format clinical information"""
    metaLabels = dict()
    f = openAnyFile(clinf)
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank on line 1\n", die = True)
    line = line.rstrip("\r\n")
    pline = re.split(delim, line)
    samples = pline[1:]
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank on line 2\n", die = True)
    line = line.rstrip("\r\n")
    pline = re.split(delim, line)
    for i in range(len(samples)):
        if not null:
            if pline[i+1] == "NULL":
                continue
        metaLabels[samples[i]] = pline[i+1]
    return(metaLabels)

def wMeta(inf, col, method = "discrete", mparams = "-;-1;0,+;1", name = None, samples = None, directory = "."):
    """write .meta format clinical information"""
    cData = rCRSData(inf)[col]
    if name == None:
        name = re.sub(" ", "", col)
    if samples == None:
        samples = cData.keys()
    if directory.endswith("/"):
        directory = directory.rstrip("/")
    f = open("%s/%s.metadata" % (directory, name), "w")
    f.write("labels\t"+"\t".join(samples)+"\n")
    vals = []
    if method == "discrete":
        labelList = []
        for i in re.split(",", mparams):
            labelList.append(re.split(";", i))
        for i in samples:
            for label, j in enumerate(labelList):
                if i not in cData:
                    vals.append("NULL")
                    break
                elif cData[i] in j:
                    vals.append(str(label))
                    break
                elif label+1 == len(labelList):
                    vals.append("NULL")
    elif method == "quartile":
        medVal = mCalculate.median(cData[i].values())
        for j in samples:
            try:
                if float(cData[i][j]) > medVal:
                    f.write("\t1")
                else:
                    f.write("\t0")
            except ValueError:
                f.write("\tNULL")
    f.write("knownVal\t"+"\t".join(vals)+"\n")
    f.close()

def rFolds(inf, limit = None):
    """read .folds format file"""
    splitMap = dict()
    f = openAnyFile(inf)
    ## Header
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank on line 1\n", die = True)
    line = line.rstrip("\r\n")
    colFeatures = re.split("\t", line)[1:]
    ## Data
    r = 1
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        pline = re.split("\t", line)
        if len(pline) != (1+len(colFeatures)):
            log("ERROR: length of line does not match the rest of the file\n", die = True)
        splitMap[r] = dict()
        for i in range(len(colFeatures)):
            if pline[i+1] not in splitMap[r]:
                splitMap[r][pline[i+1]] = set()
            splitMap[r][pline[i+1]].update([colFeatures[i]])
        if limit is not None:
            if r >= limit:
                break
        r += 1
    f.close()
    return(splitMap)

def wFolds(outf, splitMap):
    samples = []
    m = 1
    while m in splitMap[1]:
        samples += splitMap[1][m]
        m += 1
    samples.sort()
    f = open(outf, "w")
    f.write("id\t%s\n" % ("\t".join(samples)))
    r = 1
    while r in splitMap:
        f.write("repeat_%s" % (r))
        for i in samples:
            val = -1
            for m in splitMap[r].keys():
                if i in splitMap[r][m]:
                    val = m
                    break
            f.write("\t%s" % (val))
        f.write("\n")
        r += 1
    f.close()

def createSplits(metaGroups, seed = None, nrepeats = 1, mfolds = 5):
    """create splitMap from samples"""
    if seed != None:
        random.seed(seed)
    for i in metaGroups.keys():
        if i == "NULL":
            continue
        if len(metaGroups[i]) < mfolds:
            log("ERROR: Not enough samples for mfold\n", die = True)
    splitMap = dict()
    for r in range(1, nrepeats+1):
        selectGroups = deepcopy(metaGroups)
        groups = selectGroups.keys()
        sampleMap = dict()
        while len(groups) > 0:
            for m in range(1, mfolds+1):
                sampleMap[selectGroups[groups[0]].pop(random.randint(0,len(selectGroups[groups[0]])-1))] = m
                if len(selectGroups[groups[0]]) == 0:
                    groups.pop(0)
                    if len(groups) == 0:
                        break
        splitMap[r] = reverseDict(sampleMap)
    return(splitMap)

def rMAF(inf, delim = "\t", retSamples = False):
    """read .maf format file"""
    mutData = dict()
    f = openAnyFile(inf)
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank on line 1\n", die = True)
    line = line.rstrip("\r\n")
    pline = re.split(delim, line)
    hugoCol = -1
    tumorCol = -1
    for i, j in enumerate(pline):
        if j == "Hugo_Symbol":
            hugoCol = i
        elif j == "Tumor_Sample_Barcode":
            tumorCol = i
    samples = []
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\t\r\n")
        pline = re.split(delim, line)
        if pline[hugoCol] not in mutData:
            mutData[pline[hugoCol]] = list()
        mutData[pline[hugoCol]].append(pline[tumorCol])
        if pline[tumorCol] not in samples:
            samples.append(pline[tumorCol])
    f.close()
    if retSamples:
        return(mutData, samples)
    else:
        return(mutData)

def rVCF(inf, delim = "\t"):
    """read .vcf files from directory"""
    mutSet = set()
    f = openAnyFile(inf)
    for line in f:
        if line.isspace():
            continue
        if line.startswith("#"):
            continue
        line = line.rstrip("\t\r\n")
        pline = re.split(delim, line)
        gene = re.split("[=/]", pline[7])[1]
        if gene not in mutSet:
            mutSet.update([gene])
    f.close()
    return(mutSet)

def wMutData(outf, mutData, samples, features):
    """write mutation data into a paradigm rawFile"""
    f = open(outf, "w")
    f.write("id\t%s\n" % ("\t".join(features)))
    for i in samples:
        f.write("%s" % (i))
        for j in features:
            if j not in mutData:
                f.write("\t%s" % ("NA"))
            else:
                if i in mutData[j]:
                    f.write("\t%s" % ("1"))
                else:
                    f.write("\t%s" % ("0"))
        f.write("\n")
    f.close()

def reverseCRS(inData):
    outData = dict()
    for col in inData.keys():
        for row in inData[col].keys():
            if row not in outData:
                outData[row] = dict()
            outData[row][col] = inData[col][row]
    return(outData)

def getListIndices(inItem, inList):
    """returns indices of the occurence of inItem in inList"""
    indices = []
    for i, item in enumerate(inList):
        if item == inItem:
            indices.append(i)
    return(indices)

def applyData(inData, fh):
    outData = dict()
    for i in inData.keys():
        outData[i] = dict()
        for j in inData[i].keys():
            outData[i][j] = fh(inData[i][j])
    return(outData)

def randomNoise(val, scale = 0.1):
    try:
        return(float(val)*(1+scale*random.gauss(0,1)))
    except ValueError:
        return("NA")

def floatList(inList):
    """returns only numeric elements of a list"""
    outList = []
    for i in inList:
        try:
            fval = float(i)
            if fval != fval:
                raise ValueError
            outList.append(fval)
        except ValueError:
            continue
    return(outList)

def mean(inList, null = "NA"):
    """Calculates mean"""
    fList = floatList(inList)
    if len(fList) == 0:
        mean = null
    else:
        mean = sum(fList)/len(fList)
    return (mean)

def mean_std(inList, sample = True):
    """Calculates mean and std"""
    fList = floatList(inList)
    if len(fList) == 0:
        mean = "NA"
        std = "NA"
    else:
        mean = sum(fList)/float(len(fList))
        std = 0.0
        for i in fList:
            std += (i-mean)**2
        if len(fList) > 1:
            if sample:
                std = math.sqrt(std/(len(fList)-1))
            else:
                std = math.sqrt(std/len(fList))
        else:
            std = 0.0
    return(mean, std)

def corrMaps(map1, map2):
    """Calculates pearson correlation"""        
    features = list(set(map1.keys()) & set(map2.keys()))
    list1 = []
    list2 = []
    for feature in features:
        list1.append(map1[feature])
        list2.append(map2[feature])
    mean1 = mean(list1)
    mean2 = mean(list2)
    cov = 0.0
    stdev1 = 0.0
    stdev2 = 0.0
    for i in range(len(list1)):
        try:
            fval1 = float(list1[i])
            fval2 = float(list2[i])
            cov += (fval1-mean1)*(fval2-mean2)
            stdev1 += (fval1-mean1)**2
            stdev2 += (fval2-mean2)**2
        except:
            continue
    stdev1 = math.sqrt(stdev1)
    stdev2 = math.sqrt(stdev2)
    if stdev1 == 0 or stdev2 == 0:
        value = "NA"
    else:
        value = cov/(stdev1*stdev2)
    return(value)
