#!/usr/bin/env python

import logging
import sys
import argparse

from Parser.Gff.GffGeneParser import GffGeneParser
from Utils.Metrics.BamMetrics import BamMetrics

class TPMCounter(object):

    def __init__(self, bamFile, gffile, logLevel='ERROR'):

        self.gffFile = gffile
        self.bamFile = bamFile
        self.logLevel =logLevel
        logging.basicConfig(level=self.logLevel)

    def __del__(self):
        pass

    def getGenomeIndex(self, lGenes):
        """Index CDS positions for each sequence and return a dict"""

        dSeq = {}
        currentSeq = None
        lSeq = []
        lTuples = []
        for gene in sorted(lGenes):
            if gene.seqid != currentSeq:
                if currentSeq != None:
                    lSeq = [1]*(lTuples[-1][1]+5000)
                    for tup in lTuples:
                        for i in range(tup[0]-1,tup[1]):
                            lSeq[i] = 0
                    dSeq[currentSeq] = lSeq
                    logging.info("Indexing CDS for seq: {}".format(currentSeq))
                    lTuples = []
                    lSeq = []
                currentSeq = gene.seqid
            for transcript in gene.lTranscripts:
                for cds in transcript.lCDS:
                    lTuples.append((cds.start, cds.end))
        lSeq = [1]*(lTuples[-1][1]+5000)
        for tup in lTuples:
            for i in range(tup[0]-1,tup[1]):
                lSeq[i] = 0
        dSeq[currentSeq] = lSeq
        logging.info("Indexing CDS for seq: {}".format(currentSeq))
        return dSeq

                    
    def run(self, minNbFrags, stranded, countFile=None):
        """run"""
        

        # get transcript len
        GffParser = GffGeneParser(self.gffFile) 
        genomeIndex = self.getGenomeIndex(GffParser.getAllGenes())
        
        # create genome Index
       
        # get mean fragment len
        lMeans = []
        meanFragLength = 0
        dNbFragPerTranscripts = {}
        iBamMetrics = BamMetrics(self.bamFile)
        for gene in GffParser.getAllGenes():
            for t in gene.lTranscripts:
                #if t.id =="lm_SuperContig_0_v2_lmctg_0053_v2_egn4_orf_Lema_T005660.1":
                mean, nbfrags = iBamMetrics.getMeanFragmentLengthFromSplicedMapping(chr=t.seqid, start=t.start,end=t.end, chrIndex=genomeIndex[t.seqid], annotIsReverse=t.isOnReverseStrand(), stranded=stranded)
                dNbFragPerTranscripts[t.id] = nbfrags
                logging.info("{}: nb frags: {}, mean frag length: {}".format(t.id, nbfrags, mean))
                if nbfrags >= minNbFrags:
                    lMeans.append(mean)
        meanFragLength = sum(lMeans)/len(lMeans)
        logging.info("mean fragment length:{}".format(meanFragLength))


        # read countFile
        if countFile:
            dNbFragPerTranscripts = dict(zip(dNbFragPerTranscripts.keys(),[0]*len(dNbFragPerTranscripts.keys())))
            dNbFragPerTranscripts = self.readCountFile(countFile ,dNbFragPerTranscripts)

        # compute effective len
        dLengthTranscripts = {}
        dEffectiveLenghtTranscripts = {}
        for gene in GffParser.getAllGenes():
            for t in gene.lTranscripts:
                effectiveLen = max(1, t.getCDSTotalLength() - meanFragLength + 1)
                dLengthTranscripts[t.id] = t.getCDSTotalLength()
                dEffectiveLenghtTranscripts[t.id] = effectiveLen

        sumRatioCountEffectiveLength = sum ([ dNbFragPerTranscripts[t]/dEffectiveLenghtTranscripts[t] for t in dLengthTranscripts ])
        #print sumRatioCountEffectiveLength
        print "ID\tlength\tEffectiveLength\tcounts\tTPM"
        for t in dLengthTranscripts:
            TPM = (dNbFragPerTranscripts[t] / dEffectiveLenghtTranscripts[t])*(1/sumRatioCountEffectiveLength) * 1.e6 
            print "{}\t{}\t{}\t{}\t{}".format(t,dLengthTranscripts[t],dEffectiveLenghtTranscripts[t],dNbFragPerTranscripts[t],TPM)
        

    def readCountFile(self, countFile, dNbFragPerTanscript):
        """read count File"""
        print dNbFragPerTanscript
        try:
            with open(countFile, 'r') as f:
                for line in f:
                    values = line.rstrip().split("\t")
                    if values[0] not in dNbFragPerTanscript:
                        raise Exception('Read count file; Transcript: {} not in annotation file'.format(values[0]))
                    else:
                        dNbFragPerTanscript[values[0]] = float(values[1])
            return dNbFragPerTanscript
        except Exception as e:
            logging.error(e)
            sys.exit(1)
 

    def export(self):
        pass

if __name__ == '__main__':

    program = 'count_TPM.py'
    version = 0.1
    description = 'Count TPM for each transcript from a BAM file or a file with counts \
                   like FeatureCounts program. It requires an annotation file in gff \
                   format, a BAM file and a counts file (like FeatureCounts) \
                   if you want to compute TPM from a precompute counts file. \
                   The program computes the mean fragment length for each transcript\
                   and the mean for all. Then, it computes the effective\
                   length and the TPM'
    parser = argparse.ArgumentParser(prog=program)
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("BamFile", help="BamFile", type=str)
    parser.add_argument("AnnotFile", help="Annotation file in gff3 format", type=str)

    parser.add_argument("-m","--minNbFrags", help="Minimum number of fragment per trancript to \
                        compute mean fragment length, [default=3]", type=int,
                        default=3) 
    parser.add_argument("-s","--stranded", help="Define which fragments take into account\
                        for computing mean fragment length and counts. If your reads are \
                        not oriented, set to \"no\" as default, if not use \"yes\" or \
                        \"reverse\". [default=no]", type=str, default="no")
    parser.add_argument("-c","--countFile", help="file with counts (as Htseq-count, ...) \
                        instead of internal count computation. Format: at least 2 columns \
                        with no header")
    parser.add_argument("--version", action='version', version='{} {}'.format(program,version))
    parser.add_argument("-v", "--verbosity", type=int, choices=[1,2,3],
                        help="increase output verbosity 1=error, 2=info, 3=debug")

    args = parser.parse_args()

    logLevel='ERROR'
    if args.verbosity == 1:
        logLevel = 'ERROR'
    if args.verbosity == 2:
        logLevel = 'INFO'
    if args.verbosity == 3:
        logLevel = 'DEBUG'
    logging.getLogger().setLevel(logLevel)


    i = TPMCounter(args.BamFile, args.AnnotFile, logLevel)
    i.run(args.minNbFrags, args.stranded, countFile=args.countFile)
