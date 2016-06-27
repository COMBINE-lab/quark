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

class island:
    start = 0
    end = 0
    def create(self,s,e):
        self.start = s
        self.end = e
    def updateboth(self,s,e):
        self.start = s
        self.end = e
    def updatestart(self,s):
        self.start = s
    def updateend(self,e):
        self.end = e
    def get(self):
        return (self.start,self.end)


def readEqClass(eqfile, fastaFile,f1,f2,oFile,threshold):
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
            if len(readnames) > 20:
                encode(readnames,fastaFile,f1,f2,oFile,tnames)
            else:
                with open(oFile,'a') as wFile:
                    for rid in readnames:
                        headerno=int(rid.split(" ")[0].split("_")[1])
                        import linecache
                        seq = linecache.getline(f1,headerno*2).rstrip()
                        seqm = linecache.getline(f2,headerno*2).rstrip()
                        wFile.write("{}:{}\n".format(seq,seqm))

    print "\n"

def discoverIsland(tseq,boundaries):
    islands = []
    islandDict = {}
    for (h,b,e) in boundaries:
        #print b,e
        if len(islands) == 0:
            i = island()
            #print "create new island with ",b,e
            if b >= 0:
                i.create(b,e)
            else:
                i.create(0,e)
            islands.append(i)
            islandDict[h] = 0
            #modified.append((h,b,e,0))
        else:
            cs,ce = islands[-1].get()
            #print cs,ce
            iid = len(islands)
            #no need to create new
            if(b < 0):
                #modified.append((h,b,e,0))
                islandDict[h] = 0
            elif(b >= cs and e <= ce):
                #modified.append((h,b,e,iid-1))
                islandDict[h] = iid-1
            #update existing
            elif(b <= ce and e > ce):
                islands[-1].updateend(e)
                #modified.append((h,b,e,iid-1))
                islandDict[h] = iid-1
            #create new
            elif(b > ce):
                i = island()
                i.create(b,e)
                islands.append(i)
                #modified.append((h,b,e,iid))
                islandDict[h] = iid
        cs,ce = islands[-1].get()
        if ce > len(tseq):
            islands[-1].updateend(len(tseq)-1)

    return (islands,islandDict)



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

def manageNeg(pos,seqr,tseq):
    overhang = seqr[0:abs(pos)]
    i = 0
    while((abs(pos)+i < len(seqr))and(tseq[i] == seqr[abs(pos)+i])):
        #print tseq[i],seqr[abs(pos)+i],i,abs(pos)+i,len(seqr)
        i = i + 1
    matches = i
    if i < len(seqr):
        leftstr = overhang+"M"+str(matches)+seqr[i:]
        #wFile.write("{}M{}{}\n".format(overhang,matches,seqr[i:]))
    else:
        leftstr = overhang + "M" + str(matches)

    return leftstr


def manageRight(pos,seqr,tseq):
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

    return encoded

def manageLeft(pos,seqr,prevseq,currseq,prevpos,currpos,tseq):
    i = 0
    overhang = ''
    match = 0
    counter = 0
    encoded = ''
    wstr = ''
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

    return wstr


def encode(readnames,fastaFile,f1,f2,oFile,tname):
    tid = int(readnames[0].split(" ")[3])
    import linecache
    tseq = linecache.getline(fastaFile,tid*2+2).rstrip()

    #print "tsq: ",tseq[:20]
    #sort the read names
    leftreads = []
    rightreads = []
    for rid in readnames:
        headerno=int(rid.split(" ")[0].split("_")[1])
        seq = linecache.getline(f1,headerno*2).rstrip()
        seqR = linecache.getline(f2,headerno*2).rstrip()
        leftreads.append((rid.split(" ")[0],seq,int(rid.split(" ")[5]),int(rid.split(" ")[2]), \
            seqR,int(rid.split(" ")[4]),int(rid.split(" ")[1])))
        #rightreads.append((rid.split(" ")[0],seqR,int(rid.split(" ")[4]),int(rid.split(" ")[1])))
        #print rid.split(" ")[0],int(rid.split(" ")[2])

    import operator
    sortedl = sorted(leftreads,key=operator.itemgetter(3))
    #discover islands
    boundaries =[(h,p,p+len(s)) for h,s,_,p,_,_,_ in sortedl]
    islands, islandDict  = discoverIsland(tseq,boundaries)

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
    lastId = -1
    currId = -1


    with open(oFile,'a') as wFile:
        for i in islands:
            #sys.stdout.write(tseq[i.start:(i.end-i.start+1)])
            wFile.write("{}\t".format(tseq[i.start:i.end+1]))
        wFile.write("\n")


    with open(oFile,'a') as wFile:
        #wFile.write("{}\t{}\t{}\n".format(tseq[minpos:(maxpos-minpos)+54],minpos,maxpos))
        for header,seq,fwd,pos,seqM,fwdM,posM in sortedl:
            #print header
            s,e = islands[islandDict[header]].get()
            relpos=0
            relposM = 0
            if pos >= 0:
                relpos = pos-s
            else:
                relpos = pos

            if posM >= 0:
                relposM = posM-s
            else:
                relposM = posM

            currId = islandDict[header]
            if lastId == -1:
                wFile.write("{}\t{}\t".format(relpos,islandDict[header]))
            elif lastId == currId:
                wFile.write("{}\t".format(relpos))
            elif lastId != currId:
                wFile.write("{}\t{}\t".format(relpos,islandDict[header]))

            lastId = currId

            lOrphan = False
            rOrphan = False

            seqr = revcomp(seq) if(not(fwd)) else seq
            seqm = revcomp(seqM) if(not(fwdM)) else seqM

            #for left end
            overhang = ''
            currseq = seqr
            currpos = pos
            wstr = ''
            #for right end
            overhangM = ''
            currseqM = seqm
            currposM = posM
            wstr = ''
            #print abs(pos)
            #if it is unmapped only hope is to encode like previous

            #start working on the right end Now it becomes more complex
            #because either of the ends might be unmapped handle both ends
            #differently then concatenate them use flags to control them
            #leftstr and rightstr will be touched once and once only

            rightstr = ''
            leftstr = ''


            if(pos == 0 or posM == 0):
                #get done with orphan read
                if(pos == 0):
                    lOrphan = True
                    if(noRef):
                        ref = seqr
                        #wFile.write("{}{}".format(seqr,fwd))
                        leftstr = "{}{}".format(seqr,fwd)
                        noRef = False
                    else:
                        i = 0
                        while((i<len(seqr))and(ref[i] == seqr[i])):
                            i = i+1
                        if(i > len(seqr)-2):
                            #wFile.write("M{}{}".format(i,fwd))
                            leftstr="M{}{}".format(i,fwd)
                        else:
                            leftstr = "{}{}".format(seqr,fwd)
                            ref=seqr
                if(posM == 0):
                    rOrphan = True
                    rightstr = "{}".format(seqm)

            if(pos < 0):
                leftstr = manageNeg(pos,seqr,tseq)
            if(posM < 0):
                rightstr = manageNeg(posM,seqm,tseq)


            if(pos > 0):
                #logic for writing L which blindly
                #encode the previous line
                #call manageLeft(pos,seqr,prevseq,currseq,tseq):
                leftstr = manageLeft(pos,seqr,prevseq,currseq,prevpos,currpos,tseq)
                if prev == '':
                    prev = leftstr
                    leftstr = leftstr+str(fwd)
                    #wFile.write("{}{}\n".format(wstr,fwd))
                else:
                    if leftstr == prev:
                        #wFile.write("{}\n".format(fwd))
                        leftstr = str(fwd)
                    else:
                        #wFile.write("{}{}\n".format(wstr,fwd))
                        leftstr = leftstr+str(fwd)
                        prev = leftstr
                prevseq = currseq
                prevpos = currpos

            if(posM > 0):
                rightstr = manageRight(posM,seqm,tseq)

            wFile.write("{}:{}\t{}{}\n".format(leftstr,relposM,rightstr,fwdM))






def main():
    import argparse
    parser = argparse.ArgumentParser(description="blah blah")
    parser.add_argument('-i',type=str,help="eq class file")
    parser.add_argument('-tr',type=str,help="transcript file")
    parser.add_argument('-f1',type=str,help="fastq file")
    parser.add_argument('-f2',type=str,help="fastq file")
    parser.add_argument('-t',type=int,help="threshold value")
    parser.add_argument('-o',type=str,help="output file")
    args = parser.parse_args()
    eqfile = args.i
    fastaFile = args.tr
    f1 = args.f1
    f2 = args.f2
    oFile = args.o
    threshold = args.t
    readEqClass(eqfile,fastaFile,f1,f2,oFile,threshold)


if __name__ == "__main__":
    main()
