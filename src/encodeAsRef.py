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
from pyfaidx import Fasta
import operator
from Bio import SeqIO
import linecache

def revcomp(s):
    comp = {'A':'T','T':'A','G':'C','C':'G','a':'t','t':'a','g':'c','c':'g','N':'N'}
    r = ''.join(comp[e] for e in s[::-1])
    return r


def readEqClass(eqfile, fastaFile,fastqFile,oFile):
    with open(eqfile) as ifile:
        numTran = int(ifile.readline().rstrip())
        numEq = int(ifile.readline().rstrip())
        print("file: {}; # tran = {}; # eq = {}".format(eqfile, numTran, numEq))
        #read the transcript names
        tnames = []
        for i in xrange(numTran):
            tnames.append(ifile.readline())

        for i in xrange(numEq):
            toks = map(str, ifile.readline().rstrip().split('\t'))
            nt = toks[0]
            tids = tuple(toks[1:int(nt)+1])
            count = toks[-1]
            readnames = tuple(toks[int(nt)+2:-1])
            if len(readnames) > 1000:
                encode(readnames,fastaFile,fastqFile,oFile,tnames)
                break

def encode(readnames,fastaFile,fastqFile,oFile,tname):
    tid = int(readnames[0].split(" ")[3])
    print tid
    import linecache
    tseq = linecache.getline(fastaFile,tid*2+2).rstrip()

    #print "tsq: ",tseq[:20]
    #sort the read names
    leftreads = []
    for rid in readnames:
        headerno=int(rid.split(" ")[0].split("_")[1])
        seq = linecache.getline(fastqFile,headerno*2).rstrip()
        leftreads.append((rid.split(" ")[0],seq,int(rid.split(" ")[5]),int(rid.split(" ")[2])))
        #print rid.split(" ")[0],int(rid.split(" ")[2])

    import operator
    sortedl = sorted(leftreads,key=operator.itemgetter(2))
    with open(oFile,'w') as wFile:
        wFile.write("{}\n".format(tseq))
        for header,seq,fwd,pos in sortedl:
            wFile.write("{}".format(abs(pos)))
            seqr = revcomp(seq) if(not(fwd)) else seq
            overhang = ''
            #print abs(pos)
            #if it is unmapped only hope is to encode like previous
            if (pos == 0):
                if(noRef):
                    ref = seq
                else:
                    i = 0
                    while()
            if(pos < 0):
                overhang = seqr[0:abs(pos)]
                i = 0
                while((abs(pos)+i < len(seqr))and(tseq[i] == seqr[abs(pos)+i])):
                    #print tseq[i],seqr[abs(pos)+i],i,abs(pos)+i,len(seqr)
                    i = i + 1
                matches = i
                if i < len(seqr):
                    wFile.write("{}M{}{}\n".format(overhang,matches,seqr[i:]))
                else:
                    wFile.write("{}M{}\n".format(overhang,matches))
            else:
                i = 0
                while((i < len(seqr))and(tseq[pos+i] == seqr[i])):
                    i = i + 1
                matches = i
                if i < len(seqr):
                    wFile.write("M{}{}\n".format(matches,seqr[i:]))
                else:
                    wFile.write("M{}\n".format(matches))



def main():
    import argparse
    parser = argparse.ArgumentParser(description="blah blah")
    parser.add_argument('-i',type=str,help="eq class file")
    parser.add_argument('-tr',type=str,help="transcript file")
    parser.add_argument('-fq',type=str,help="fastq file")
    parser.add_argument('-o',type=str,help="output file")
    args = parser.parse_args()
    eqfile = args.i
    fastaFile = args.tr
    fastqFile = args.fq
    oFile = args.o
    readEqClass(eqfile,fastaFile,fastqFile,oFile)


if __name__ == "__main__":
    main()
