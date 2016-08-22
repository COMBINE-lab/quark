import sys
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

def isInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def decodeQuark(island,left,pos):
    ind = 0
    #pos = 5
    real = ""
    #ore = left[-1]
    #left = left[:-1]
    while(ind < len(left)):
        if(left[ind] == "M"):
            digit = ""
            digitInd = ind+1
            #I have to check for numbers
            #print left[digitInd]
            while(isInt(left[digitInd])):
                #print left[digitInd],digitInd
                digit += left[digitInd]
                digitInd += 1
                if (digitInd == len(left)):
                    break
            matches = int(digit)
            matchChars = island[pos:pos+matches]
            real += matchChars
            pos += matches
            ind = digitInd
        else:
            real += left[ind]
            ind += 1
            pos += 1
    return real

import time

def decode(islandFile, quarkFile, outFile):
    lFile = outFile + "mapped.1.seq"
    rFile = outFile + "mapped.2.seq"
    wlFile = open(lFile,"w")
    wrFile = open(rFile,"w")

    with open(islandFile,"r") as iFile:
        with open(quarkFile,"r") as qFile:
            starttime = time.time()
            numEq = int(iFile.readline().rstrip())
            for eqNum in xrange(numEq):
                if(eqNum%100 == 0):
                    percent = float(eqNum)/numEq
                    hashes = '#'* int(round(percent * 50))
                    dashes = ' '* (50 - len(hashes))
                    elapsedtime = float(time.time() - starttime)
                    #print type(elapsedtime)
                    sys.stdout.write("\rProgress [{}] {}% {} seconds elapsed".format(hashes+dashes, int(round(percent * 100)),elapsedtime))
                    sys.stdout.flush()
                numIslands = int(iFile.readline().rstrip())
                islands = []
                for i in xrange(numIslands):
                    islands.append(iFile.readline().rstrip())

                numSeq = int(qFile.readline().rstrip())
                for i in xrange(numSeq):
                    lDecoded = ""
                    rDecoded = ""
                    bothSeq,metaData = qFile.readline().rstrip().split("\t")
                    #bothSeq,metaData = qFile.readline().rstrip().strip("\t")
                    #return
                    #print bothSew
                    left,right = bothSeq.split("|")
                    lIsland,lPos,rIsland,rPos = map(int, metaData.split(","))
                    #decode the left end first
                    #ignore if it maps to $
                    if(islands[lIsland] == "$"):
                        lDecoded = left
                    else:
                        #need decoding
                        ore = left[-1]
                        orebool = False if (ore == "0") else True
                        lDecoded = decodeQuark(islands[lIsland],left[:-1],lPos)
                        #if(orebool):
                            #lDecoded = revcomp(lDecoded)

                    if(islands[rIsland] == "$"):
                        rDecoded = right
                    else:
                        #need decoding
                        ore = right[-1]
                        orebool = False if (ore == "0") else True
                        rDecoded = decodeQuark(islands[rIsland],right[:-1],rPos)
                        #if(orebool):
                            #rDecoded = revcomp(rDecoded)
                    if(len(lDecoded) != 54):
                        print "\n---------"
                        print left
                        print islands[lIsland]
                        print lIsland
                        print lPos
                        print "---------"
                        break
                    wlFile.write("{}\n".format(lDecoded))
                    wrFile.write("{}\n".format(rDecoded))
    wlFile.close()
    wrFile.close()

def main():
    import argparse
    parser = argparse.ArgumentParser(description="blah blah")
    parser.add_argument('-i',type=str,help="Island File")
    parser.add_argument('-q',type=str,help="Quark File")
    parser.add_argument('-o',type=str,help="Output file")
    args = parser.parse_args()
    #eqPath = "/mnt/scratch1/hirak/RapCompressData/sailfish/sailfish_quant/aux/eq_classes.txt"
    #eqPath = args.i
    decode(args.i,args.q,args.o)
    #oGraphFile = open("/mnt/scratch1/hirak/RapCompressData/graph.net",'w')
    #comb = itertools.combinations(liX,2)
    #print ("Done with all pairs")
    #parallelize_func(comb, writeEdges,1000,20,eqPath,groups)
    #Parallel(n_jobs=30)(delayed(writeEdges)(i,j,oGraphFile,groups) for i,j in comb)


if __name__ == "__main__":
    main()
