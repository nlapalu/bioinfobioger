#!/usr/bin/env python

import numpy as np
import pysam

class BamMetrics(object):

    def  __init__(self, bamFile, logLevel='ERROR'):

        try:
            self.bamFileH = pysam.AlignmentFile(bamFile, 'r', check_sq=True)
        except:
            raise

    def getMeanFragmentLength(self, chr='', start=None, end=None):
        """
            return mean fragment length

            :param self: object ref
            :param chr: chromosome name
            :param start: chromosome start 1-based indexing
            :param end: chromosome end 1-based indexing

            If chr, start, end not specified, return mean fragment across all
            the file.

        """

        lFragLength = []
        lFragLenM = []

        if (chr and start and end):
            readIterator = self.bamFileH.fetch(reference=chr, start=start, end=end)
        else:
            readIterator = self.bamFileH.fetch()

        for read in readIterator:
            # process paired reads
            if read.is_paired:
                if read.is_read1:
                    if (read.template_length < 0):
                        fragStart = read.next_reference_start
                        fragEnd = read.next_reference_start-read.template_length -1
                    else:
                        fragStart = read.reference_start
                        fragEnd = read.reference_start + read.template_length -1 
            # process single
            else:
                fragStart = read.reference_start
                fragEnd = read.reference_start + read.query_alignment_length -1
 
            lFragLength.append(fragEnd-fragStart+1)
            if len(lFragLength) == 1000000:
                lFragLenM.append(sum(lFragLength)/1000000)
                lFragLength = []
                
        return sum(lFragLenM)/len(lFragLenM)


    def getMeanFragmentLengthFromSplicedMapping(self, chr='', start=None, end=None, chrIndex=None, annotIsReverse=None, stranded="no"):
        """
            return mean fragment length for spliced mapping

        """

        lFragLength = []
        lFragLenM = []
        lFragCoords = []
        nbFragments = 0
        fragStart = 0
        fragEnd = 0

        if (chr and start and end):
            readIterator = self.bamFileH.fetch(reference=chr, start=start, end=end)
        else:
            readIterator = self.bamFileH.fetch()

        for read in readIterator:
            # process paired reads
            if read.is_paired:
                if read.is_read1:
                    if (stranded == 'no') or (stranded == 'yes' and read.is_reverse == annotIsReverse) or (stranded == 'reverse' and read.is_reverse != annotIsReverse):
                        if (read.template_length < 0):
                            fragStart = read.next_reference_start
                            fragEnd = read.next_reference_start-read.template_length -1
                        else:
                            fragStart = read.reference_start
                            fragEnd = read.reference_start + read.template_length -1
                    else:
                        continue 
                else:
                    continue
            # process single
            else:
                if (stranded == 'no') or (stranded == 'yes' and read.is_reverse == annotIsReverse) or (stranded == 'reverse' and read.is_reverse != annotIsReverse):
                    #print "simple"
                    fragStart = read.reference_start
                    fragEnd = read.reference_start + read.query_alignment_length -1


            #nbBasesToExclude = self.getNbOfNonConsideredBases(chr, fragStart, fragEnd, chrIndex)
            #print "nb bases nbBasesToExclude: {}".format(nbBasesToExclude)

            #print "start:{}, end:{}, len: {}, nbtoex:{}".format(fragStart,fragEnd,fragEnd-fragStart+1,  nbBasesToExclude)

            #lFragLength.append(fragEnd-fragStart+1-nbBasesToExclude)
            #lFragLength.append(fragEnd-fragStart+1)
            lFragCoords.append((fragStart,fragEnd))


            #if len(lFragLength) == 1000000:
            if len(lFragCoords) == 1000000:
                lFragLength =  self.getFragLengthWithSplicing(lFragCoords, chrIndex)
                lFragLenM.append(sum(lFragLength)/1000000)
                lFragLength = []
                nbFragments += 1000000

        #nbFragments += len(lFragLength)
        nbFragments += len(lFragCoords)
        if nbFragments > 0: 
            lFragLength =  self.getFragLengthWithSplicing(lFragCoords, chrIndex)
            lFragLenM.append(sum(lFragLength)/len(lFragLength))
            return sum(lFragLenM)/float(len(lFragLenM)), nbFragments
        else:
            return 0,0


    def getFragLengthWithSplicing(self, lFragCoords, chrIndex):
        lFragLen = []
        for frag in lFragCoords:
            #lFragLen.append(frag[1]-frag[0]+1-sum(chrIndex[frag[0]-1:frag[1]]))
            #lFragLen.append(frag[1]-frag[0]+1-2)
            lFragLen.append(frag[1]-frag[0]+1-sum([chrIndex[x] for x in range(frag[0]-1,frag[1])]))
        return lFragLen 

    def getNbOfNonConsideredBases(self, chr, start, end, chrIndex): 
        """ ... """

        return sum(chrIndex[start-1:end])
