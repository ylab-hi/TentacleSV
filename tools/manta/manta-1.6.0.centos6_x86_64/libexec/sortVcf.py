#!/usr/bin/env python2
#
# Manta - Structural Variant and Indel Caller
# Copyright (c) 2013-2019 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

"""
sort input vcf
"""

import os, sys
import re


def isInfoKey(string,key) :
    match=re.search("(^|;)%s(;|$)" % (key) ,string)
    return match is not None


def getKeyVal(string,key) :
    match=re.search("%s=([^;]*);?" % (key) ,string)
    if match is None : return None
    return match.group(1);


VCF_CHROM = 0
VCF_POS = 1
VCF_ID = 2
VCF_REF = 3
VCF_ALT = 4
VCF_QUAL = 5
VCF_FILTER = 6
VCF_INFO = 7



class VcfRecord :
    def __init__(self, line) :
        self.line = line
        w=line.strip().split('\t')
        self.chrom=w[VCF_CHROM]
        self.pos=int(w[VCF_POS])
        self.id = w[VCF_ID]
        self.ref=w[VCF_REF]
        self.alt=w[VCF_ALT]
        self.qual=w[VCF_QUAL]
        self.isPass=(w[VCF_FILTER] == "PASS")

        self.invState=None
        inv3 = isInfoKey(w[VCF_INFO],"INV3")
        inv5 = isInfoKey(w[VCF_INFO],"INV5")
        assert(not (inv3 and inv5))
        if inv3: self.invState = "INV3"
        if inv5: self.invState = "INV5"

        self.endPos=self.pos+len(w[VCF_REF])-1
        val = getKeyVal(w[VCF_INFO],"END")
        if val is not None :
            self.endPos = int(val)
        self.mateId = None
        val = getKeyVal(w[VCF_INFO],"MATEID")
        if val is not None :
            self.mateId = val


class Constants :

    import re

    contigpat = re.compile("^##contig=<ID=([^,>]*)[,>]")


def processFile(vcfFile, isFirst, chromOrder, header, recList) :
    """
    read in a vcf file
    """

    import re

    for line in open(vcfFile) :
        if line[0] == "#" :
            if not isFirst : continue
            header.append(line)
            match = re.match(Constants.contigpat,line)
            if match is not None :
                chromOrder.append(match.group(1))
        else :
            recList.append(VcfRecord(line))



def listInputVcfs(vcfListFile,args) :
    for arg in args :
        yield arg
    if vcfListFile is None : return
    for vcfFile in open(vcfListFile) :
        yield vcfFile.strip()



def getOptions() :

    from optparse import OptionParser

    parser = OptionParser(usage="%prog [options] [input_vcf [input_vcf..]] > output_vcf")

    parser.add_option("-f", dest="vcfListFile",
                      help="File listing input vcf files, one file per line. These will be used in addition to any provided directly on the command-line")

    parser.add_option("-a", dest="isPrintAll", action="store_true", default=False,
                      help="print out all vcf records without filtering")

    (options,args) = parser.parse_args()

    if len(args) == 0 and not options.vcfListFile:
        parser.print_help()
        sys.exit(2)

    # validate input:
    if options.vcfListFile is not None :
        if not os.path.exists(options.vcfListFile) :
            raise Exception("Can't find vcf list file: " + options.vcfListFile)

    for vcfFile in listInputVcfs(options.vcfListFile,args) :
        if not os.path.isfile(vcfFile) :
            raise Exception("Can't find input vcf file: " +vcfFile)

    return (options,args)



def resolveRec(recEqualSet, recList, removedMateIds) :
    """
    determine which of a set of 'equal' vcf records is the best

    right now best is a record with PASS in the filter field, and
    secondarily having the highest quality

    Keep track of MateIDs of removed variants
    """

    if not recEqualSet: return

    bestIndex=0
    bestQual=0.
    bestIsPass=False
    bestIsAssembled=False
    for (index,rec) in enumerate(recEqualSet) :
        try:
            rec.qual = float(rec.qual)
        except ValueError:
            rec.qual = 0.

        assert rec.qual >= 0.

        isNewPass=((not bestIsPass) and rec.isPass)
        isHighQual=((bestIsPass == rec.isPass) and (rec.qual > bestQual))
        isNewAssembled=((not bestIsAssembled) and (rec.alt[0] != '<'))
        if (isNewPass or isHighQual or isNewAssembled) :
            bestIndex = index
            bestQual = rec.qual
            bestIsPass = rec.isPass
            bestIsAssembled = (rec.alt[0] != '<')

    recList.append(recEqualSet[bestIndex])
    for index,rec in enumerate(recEqualSet):
        if index != bestIndex and rec.mateId is not None:
            removedMateIds.add(rec.mateId)



def main() :

    outfp = sys.stdout

    (options,args) = getOptions()

    header=[]
    recList=[]
    chromOrder=[]

    isFirst=True
    for vcfFile in listInputVcfs(options.vcfListFile,args) :
        processFile(vcfFile, isFirst, chromOrder, header, recList)
        isFirst = False

    def vcfRecSortKey(x) :
        """
        sort vcf records for final output

        Fancy chromosome sort rules:
        if contig records are found in the vcf header, then sort chroms in that order
        for any chrom names not found in the header, sort them in lexicographical order after the
        found chrom names
        """

        try :
            headerOrder = chromOrder.index(x.chrom)
        except ValueError :
            headerOrder = len(chromOrder)

        return (headerOrder, x.chrom, x.pos, x.endPos, x.ref, x.alt, x.id)

    recList.sort(key = vcfRecSortKey)

    for line in header :
        outfp.write(line)

    def isEqualRec(rec1,rec2) :
        if (rec1 is None) or (rec2 is None) : return False

        if rec1[0] != rec2[0]: return False # chrom
        if rec1[1] != rec2[1]: return False # pos
        if rec1[2] != rec2[2]: return False # ref
        if rec1[4] != rec2[4]: return False # endPos
        if rec1[5] != rec2[5]: return False # invState

        # special handling to find duplications when alt is the only difference:
        if rec1[3] != rec2[3]:
            if rec1[3] != "<INS>" and rec2[3] != "<INS>":
                return False

            def matchTest(rec) :
                if rec[0] == "<" : return False
                if len(rec) < 80 : return False
                return True

            if rec1[3] == "<INS>" :
                return matchTest(rec2[3])
            if rec2[3] == "<INS>" :
                return matchTest(rec1[3])

        return True

    idsToRemove = set() #IDs for BNDs we want to remove since their mates have been removed
    if not options.isPrintAll:
        recList2 = []
        recEqualSet = []
        lastRec = None
        for vcfrec in recList :
            if vcfrec.id in idsToRemove: continue #Remove a variant if its mate is already removed. Ensures that consistent BND pairs are selected
            rec = (vcfrec.chrom, vcfrec.pos, vcfrec.ref, vcfrec.alt, vcfrec.endPos, vcfrec.invState)
            if not isEqualRec(rec,lastRec) :
                resolveRec(recEqualSet,recList2, idsToRemove)
                recEqualSet = []
            recEqualSet.append(vcfrec)
            lastRec = rec
        resolveRec(recEqualSet,recList2, idsToRemove)
        recList = recList2

    for vcfrec in recList :
        if vcfrec.id in idsToRemove: continue # Remove BNDs that were kept originally but their mates were removed later in the loop
        outfp.write(vcfrec.line)


main()

