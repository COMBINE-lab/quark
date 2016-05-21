#*****************************************************************************
# Take equivalence classes and make graph out of it
# Input : Eq Class
# Output: 1. A .net file with the graph
# 2. A file containing read ids with eq2readID.txt with the following
# format,
# <eq class id>
# <number of read ids>
# <id1>
# <id2>
# ...
#****************************************************************************
from functools import partial
import numpy as np
import cPickle as pickle

from joblib import Parallel, delayed
import more_itertools
import random
from scipy.spatial.distance import cosine
from collections import defaultdict
from scipy.spatial.distance import cosine
import pandas as pd
import os
import cPickle as pickle

#import psutil
from multiprocessing import cpu_count
import itertools

#p = psutil.Process(os.getpid())
#p.set_cpu_affinity(list(range(cpu_count())))

count = 0

class EquivCollection(object):
    def __init__(self):
        self.tnames = []
        self.eqClasses = {}
        self.hasNames = False

    def setNames(self, names):
        self.tnames = names
        self.hasNames = True

    def add(self, tids, count):
        if tids in self.eqClasses:
            self.eqClasses[tids] += count
        else:
            self.eqClasses[tids] = count

def readEqClass(eqfile, eqCollection, oFile):
    with open(oFile,'w') as wfile:
        with open(eqfile) as ifile:
            numTran = int(ifile.readline().rstrip())
            numEq = int(ifile.readline().rstrip())
            print("file: {}; # tran = {}; # eq = {}".format(eqfile, numTran, numEq))
            if not eqCollection.hasNames:
                tnames = []
                for i in xrange(numTran):
                    tnames.append(ifile.readline().rstrip())
                eqCollection.setNames(tnames)
            else:
                for i in xrange(numTran):
                    ifile.readline()

            for i in xrange(numEq):
                toks = map(str, ifile.readline().rstrip().split('\t'))
                nt = toks[0]
                tids = tuple(toks[1:int(nt)+1])
                count = toks[-1]
                eqCollection.add(tids, count)
                readnames = tuple(toks[int(nt)+2:-1])
                print>>wfile,i
                print>>wfile,len(readnames)
                for rid in readnames:
                    print>>wfile,rid

import itertools
from joblib import Parallel, delayed

def transformClasses(equiv, gFile):
    from collections import defaultdict
    tr2eq = defaultdict(list)
    classDict = {}
    weightedDict = {}
    for gkey,tids in enumerate(equiv.eqClasses.keys()):
        classDict[gkey] = tids
        for t in tids:
            tr2eq[t].append(gkey)

    for t,gkeys in tr2eq.iteritems():
        if len(gkeys) > 1:
            for g1,g2 in itertools.permutations(gkeys,2):
                if (g1,g2) in weightedDict:
                    weightedDict[(g1,g2)] += 1
                else:
                    weightedDict[(g1,g2)] = 1
        else:
            k = gkeys[0]
            weightedDict[(k,k)] = 1
    with open(gFile,'w') as ofWriter:
        for gpair,count in weightedDict.iteritems():
            denom = len(classDict[gpair[0]]) + len(classDict[gpair[1]])
            denom = float(len(set(classDict[gpair[0]]).union(set(classDict[gpair[1]]))))
            eW = float(count)/float(denom)
            ofWriter.write("{}\t{}\t{}\n".format(gpair[0],gpair[1],eW))


def writeEdges(chunk,eqPath,groups):
    with open(eqPath,'w') as oGraphFile:
        for i,j in chunk:
            g1 = groups[i]
            g2 = groups[j]
            eW = float(len(set(g1).intersection(set(g2))))/float(len(set(g1).union(set(g2))))
            if eW > 0.0:
                #global count
                #count += 1
                import sys
                #if count%10000 == 0:
                #    sys.stdout.write("\r\rSaw %d reads"%(count))
                #    sys.stdout.flush()
                oGraphFile.write("{}\t{}\t{}\n".format(i,j,eW))

def parallelize_func(iterable, func, chunksz, n_jobs=20, *args, **kwargs):
    """ Parallelize a function over each element of an iterable. """
    chunker = func
    chunks = more_itertools.chunked(iterable, chunksz)
    chunks_results = Parallel(n_jobs=n_jobs, verbose=50)(
        delayed(chunker)(chunk, *args, **kwargs) for chunk in chunks)
    results = more_itertools.flatten(chunks_results)
    return list(results)

def main():
    import argparse
    parser = argparse.ArgumentParser(description="blah blah")
    parser.add_argument('-i',type=str,help="eq class file")
    parser.add_argument('-g',type=str,help="Graph file")
    parser.add_argument('-o',type=str,help="File with ids")
    args = parser.parse_args()
    #eqPath = "/mnt/scratch1/hirak/RapCompressData/sailfish/sailfish_quant/aux/eq_classes.txt"
    eqPath = args.i
    equiv = EquivCollection()
    readEqClass(eqPath, equiv, args.o)
    transformClasses(equiv,args.g)
    groups = equiv.eqClasses.keys()
    liX = xrange(len(groups))
    #oGraphFile = open("/mnt/scratch1/hirak/RapCompressData/graph.net",'w')
    #comb = itertools.combinations(liX,2)
    #print ("Done with all pairs")
    #parallelize_func(comb, writeEdges,1000,20,eqPath,groups)
    #Parallel(n_jobs=30)(delayed(writeEdges)(i,j,oGraphFile,groups) for i,j in comb)


if __name__ == "__main__":
    main()
