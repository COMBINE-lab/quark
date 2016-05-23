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
import sys
import time

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

        with open(oFile,'w'):
            pass

        eqclassnum = 0
        starttime = time.time()
        for i in xrange(numEq):
            if(eqclassnum%100 == 0):
                percent = float(eqclassnum)/numEq
                hashes = '#'* int(round(percent * 50))
                dashes = ' '* (50 - len(hashes))
                elapsedtime = float(time.time() - starttime)
                #print type(elapsedtime)
                sys.stdout.write("\rProgress [{}] {}% {} seconds elapsed".format(hashes+dashes, int(round(percent * 100)),elapsedtime))
                sys.stdout.flush()
            eqclassnum += 1
            toks = map(str, ifile.readline().rstrip().split('\t'))
            nt = toks[0]
            tids = tuple(toks[1:int(nt)+1])
            count = toks[-1]
            readnames = tuple(toks[int(nt)+2:-1])
            if len(readnames) > 1000:
                encode(readnames,fastaFile,fastqFile,oFile,tnames)
            elif len(readnames) <= 1000 and len(readnames)>2:
                encodeAsShift(readnames,fastqFile,oFile)
            else:
                with open(oFile,'a') as wFile:
                    for rid in readnames:
                        headerno=int(rid.split(" ")[0].split("_")[1])
                        import linecache
                        seq = linecache.getline(fastqFile,headerno*2).rstrip()
                        wFile.write("{}\n".format(seq))

    print "\n"



def encodeAsShift(readnames,fastqFile,oFile):
    import linecache
    leftreads = []
    for rid in readnames:
        headerno=int(rid.split(" ")[0].split("_")[1])
        seq = linecache.getline(fastqFile,headerno*2).rstrip()
        leftreads.append((rid.split(" ")[0],seq,int(rid.split(" ")[5]),int(rid.split(" ")[2])))

    import operator
    sortedl = sorted(leftreads,key=operator.itemgetter(3))
    prevseq = ''
    currseq = ''
    prevpos = 0
    currpos = 0
    with open(oFile,'a') as wFile:
        for header,seq,fwd,pos in sortedl:
            seqr = revcomp(seq) if(not(fwd)) else seq
            currseq = seqr
            currpos = pos
            wstr = ''

            if(prevseq == ''):
                prevseq = seqr
                prevpos = pos
                wFile.write("{}{}\n".format(seqr,fwd))
                continue
            else:
                shift = currpos - prevpos
                counter = 0
                match = 0
                encodeshift = False
                if(shift > 0 and shift < len(prevseq)):
                    wstr += 'S'+str(shift)
                else:
                    shift = 0

                while((shift+counter <= len(prevseq)-1) and (counter < len(currseq))):
                    if(prevseq[shift+counter] == currseq[counter]):
                        match = match + 1
                    else:
                        if(match < 3):
                            if match == 0:
                                wstr += currseq[counter]
                            else:
                                wstr += currseq[counter - match : match+1]
                        else:
                            wstr += 'M'+str(match)+currseq[counter]
                        match = 0
                    counter += 1
                if(match > 0):
                    wstr += 'M'+str(match)
                if(counter < len(currseq)):
                    wstr += currseq[counter:]
            prevseq = currseq
            prevpos = currpos
            wFile.write("{}{}\n".format(wstr,fwd))


def encode(readnames,fastaFile,fastqFile,oFile,tname):
    tid = int(readnames[0].split(" ")[3])
    #print tid
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
    sortedl = sorted(leftreads,key=operator.itemgetter(3))
    minpos = min(leftreads,key=operator.itemgetter(3))[3]
    maxpos = max(leftreads,key=operator.itemgetter(3))[3]
    if(minpos < 0):
        minpos = 0
    if(maxpos > len(tseq)-54):
        maxpos = len(tseq)


    noRef = True
    ref = ''
    prev = ''
    wstr = ''
    prevpos = 0
    currpos = 0
    prevseq = ''
    currseq = ''
    with open(oFile,'a') as wFile:
        wFile.write("{}\t{}\t{}\n".format(tseq[minpos:(maxpos-minpos)+54],minpos,maxpos))
        for header,seq,fwd,pos in sortedl:
            wFile.write("{}\t".format(pos))
            seqr = revcomp(seq) if(not(fwd)) else seq
            overhang = ''
            currseq = seqr
            currpos = pos
            wstr = ''
            #print abs(pos)
            #if it is unmapped only hope is to encode like previous
            if(pos == 0):
                if(noRef):
                    ref = seqr
                    wFile.write("{}{}\n".format(seqr,fwd))
                    noRef = False
                else:
                    i = 0
                    while((i<len(seqr))and(ref[i] == seqr[i])):
                        i = i+1
                    if(i > len(seqr)-2):
                        wFile.write("M{}{}\n".format(i,fwd))
                    else:
                        wFile.write("{}{}\n".format(seqr,fwd))
                        ref=seqr
                continue

            if(pos < 0):
                overhang = seqr[0:abs(pos)]
                i = 0
                while((abs(pos)+i < len(seqr))and(tseq[i] == seqr[abs(pos)+i])):
                    #print tseq[i],seqr[abs(pos)+i],i,abs(pos)+i,len(seqr)
                    i = i + 1
                matches = i
                if i < len(seqr):
                    wstr = overhang+"M"+str(matches)+seqr[i:]
                    #wFile.write("{}M{}{}\n".format(overhang,matches,seqr[i:]))
                else:
                    wstr = overhang + "M" + str(matches)
                    #wFile.write("{}M{}\n".format(overhang,matches))
            else:
                i = 0
                overhang = ''
                match = 0
                counter = 0
                encoded = ''
                while(pos+counter <= len(tseq)-1 and counter < len(seqr)):
                    #if tid == 8719:
                    #    print len(seqr),seqr,pos,len(tseq)
                    if(tseq[pos+counter] == seqr[counter]):
                        match += 1
                    else:
                        if(match < 3):
                            if match == 0:
                                encoded += seqr[counter]
                            else:
                                encoded += seqr[counter-match:match+1]
                        else:
                            encoded += "M" + str(match) + seqr[counter]
                        match = 0

                    counter += 1
                if(match > 0):
                    encoded += "M"+str(match)
                if(counter < len(seqr)):
                    encoded += seqr[counter:]

                score = len(encoded)
                #in this stage we have match statistics
                #also the overhang statistics
                #we can either go with shit theory or
                #we can just encode it using position
                #with respect to transcript
                if (match < len(seqr)):
                    if((score > (len(seqr)* 0.8)) and (prevseq != '')):
                        #we have to decide either to go with refbased
                        #encoding or shift encoding with respect to
                        #previous sequence
                        shift = currpos - prevpos
                        counter = 0
                        match = 0
                        encodeshift = False
                        if(shift > 0 and shift < len(prevseq)):
                            wstr += 'S'+str(shift)
                        else:
                            shift = 0

                        while((shift+counter <= len(prevseq)-1) and (counter < len(currseq))):
                            if(prevseq[shift+counter] == currseq[counter]):
                                match = match + 1
                            else:
                                if(match < 3):
                                    if match == 0:
                                        wstr += currseq[counter]
                                    else:
                                        wstr += currseq[counter - match : match+1]
                                else:
                                    wstr += 'M'+str(match)+currseq[counter]
                                match = 0
                            counter += 1
                        if(match > 0):
                            wstr += 'M'+str(match)
                    else:
                        wstr = encoded
                        #wstr = "M"+str(matches)+seqr[i:]
                        #wFile.write("M{}{}\n".format(matches,seqr[i:]))
                else:
                    wstr = encoded
                    #wFile.write("M{}\n".format(matches))

                #logic for writing L which blindly
                #encode the previous line
                if prev == '':
                    prev = wstr
                    wFile.write("{}{}\n".format(wstr,fwd))
                else:
                    if wstr == prev:
                        wFile.write("L\n")
                    else:
                        wFile.write("{}{}\n".format(wstr,fwd))
                        prev = wstr

                prevseq = currseq
                prevpos = currpos



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
