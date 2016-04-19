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

def readEqClass(eqfile, eqCollection):
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
            toks = map(int, ifile.readline().rstrip().split('\t'))
            nt = toks[0]
            tids = tuple(toks[1:-1])
            count = toks[-1]
            eqCollection.add(tids, count)

import itertools
from joblib import Parallel, delayed

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
    eqPath = "/mnt/scratch1/hirak/RapCompressData/sailfish/sailfish_quant/aux/eq_classes.txt"
    equiv = EquivCollection()
    readEqClass(eqPath, equiv)
    groups = equiv.eqClasses.keys()
    liX = xrange(len(groups))
    oFile ="/mnt/scratch1/hirak/RapCompressData/graph.net"
    comb = itertools.combinations(liX,2)
    print ("Done with all pairs")
    parallelize_func(comb, writeEdges,1000,20,oFile,groups)
    #Parallel(n_jobs=30)(delayed(writeEdges)(i,j,oGraphFile,groups) for i,j in comb)
    parallelize_func(comb,writeEdges,)


if __name__ == "__main__":
    main()
