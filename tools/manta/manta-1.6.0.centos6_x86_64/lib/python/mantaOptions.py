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
Workflow configuration options shared by multiple
configuration scripts.
"""

import os,sys

scriptDir=os.path.abspath(os.path.dirname(__file__))
scriptName=os.path.basename(__file__)

sys.path.append(scriptDir)

from checkChromSet import getFastaInfo, getTabixChromSet
from configureOptions import ConfigureWorkflowOptions
from configureUtil import assertOptionExists, checkFixTabixIndexedFileOption, joinFile, OptParseException, \
                          validateFixExistingFileArg
from workflowUtil import exeFile, parseGenomeRegion



def cleanLocals(locals_dict) :
    """
    When passed a locals() dictionary, clean out all of the hidden keys and return
    """

    return dict((k,v) for (k,v) in locals_dict.items() if not k.startswith("__") and k != "self")



class MantaWorkflowOptionsBase(ConfigureWorkflowOptions) :

    def addWorkflowGroupOptions(self,group) :
        group.add_option("--referenceFasta",type="string",metavar="FILE",
                         help="samtools-indexed reference fasta file [required]")
        group.add_option("--runDir", type="string", metavar="DIR",
                         help="Name of directory to be created where all workflow scripts and output will be written. "
                              "Each analysis requires a separate directory. (default: %default)")
        group.add_option("--callRegions", dest="callRegionsBed", metavar="FILE",
                         help="Optionally provide a bgzip-compressed/tabix-indexed BED file containing the set of regions to call. "
                              "No VCF output will be provided outside of these regions. The full genome will still be used "
                              "to estimate statistics from the input (such as expected fragment size distribution). "
                              "Only one BED file may be specified. (default: call the entire genome)")


    def addExtendedGroupOptions(self,group) :
        group.add_option("--scanSizeMb", dest="scanSizeMb", type="int", metavar="INT",
                         help="Maximum sequence region size (in megabases) scanned by each task during "
                         "SV Locus graph generation. (default: %default)")
        group.add_option("--region", type="string",dest="regionStrList",metavar="REGION", action="append",
                         help="Limit the analysis to a region of the genome for debugging purposes. "
                              "If this argument is provided multiple times all specified regions will "
                              "be analyzed together. All regions must be non-overlapping to get a "
                              "meaningful result. Examples: '--region chr20' (whole chromosome), "
                              "'--region chr2:100-2000 --region chr3:2500-3000' (two regions)'. If this "
                              "option is specified (one or more times) together with the --callRegions BED file, then "
                              "all region arguments will be intersected with the callRegions BED track.")
        group.add_option("--callMemMb",dest="callMemMbOverride",type="int",metavar="INT",
                         help="Set default task memory requirement (in megabytes) for common tasks. This may benefit "
                              "an analysis of unusual depth, chimera rate, etc.. 'Common' tasks refers to most "
                              "compute intensive scatter-phase tasks of graph creation and candidate generation.")

        ConfigureWorkflowOptions.addExtendedGroupOptions(self,group)


    def getOptionDefaults(self) :
        """
        Set option defaults.

        Every local variable in this method becomes part of the default hash
        """

        configCommandLine=sys.argv

        libexecDir=os.path.abspath(os.path.join(scriptDir,"../../libexec"))
        assert os.path.isdir(libexecDir)

        bgzipBin=joinFile(libexecDir,exeFile("bgzip"))
        htsfileBin=joinFile(libexecDir,exeFile("htsfile"))
        tabixBin=joinFile(libexecDir,exeFile("tabix"))
        samtoolsBin=joinFile(libexecDir,exeFile("samtools"))

        mantaStatsBin=joinFile(libexecDir,exeFile("GetAlignmentStats"))
        mantaMergeStatsBin=joinFile(libexecDir,exeFile("MergeAlignmentStats"))
        getChromDepthBin=joinFile(libexecDir,exeFile("GetChromDepth"))
        mantaGraphBin=joinFile(libexecDir,exeFile("EstimateSVLoci"))
        mantaGraphMergeBin=joinFile(libexecDir,exeFile("MergeSVLoci"))
        mantaGraphCheckBin=joinFile(libexecDir,exeFile("CheckSVLoci"))
        mantaHyGenBin=joinFile(libexecDir,exeFile("GenerateSVCandidates"))
        mantaGraphStatsBin=joinFile(libexecDir,exeFile("SummarizeSVLoci"))
        mantaStatsSummaryBin=joinFile(libexecDir,exeFile("SummarizeAlignmentStats"))

        mergeChromDepth=joinFile(libexecDir,"mergeChromDepth.py")
        mantaSortVcf=joinFile(libexecDir,"sortVcf.py")
        mantaExtraSmallVcf=joinFile(libexecDir,"extractSmallIndelCandidates.py")
        mantaPloidyFilter=joinFile(libexecDir,"ploidyFilter.py")
        mantaSortEdgeLogs=joinFile(libexecDir,"sortEdgeLogs.py")
        catScript=joinFile(libexecDir,"cat.py")
        vcfCmdlineSwapper=joinFile(libexecDir,"vcfCmdlineSwapper.py")
        mantaSortBam=joinFile(libexecDir,"sortBam.py")
        mantaMergeBam=joinFile(libexecDir,"mergeBam.py")

        workflowScriptName = "runWorkflow.py"

        # default memory request per process-type
        #
        # where different values are provided for SGE and local runs note:
        #  1. for SGE the memory limits must be greater than the highest memory use ever
        #      expected in a production run. The consequence of exceeding this limit is a failed
        #      run.
        #   2. for localhost the memory usage should be at least above the highest mean memory
        #       use ever expected in a production run. The consequence of exceeding the mean is
        #       a slow run due to swapping.
        #
        estimateMemMb=1.5*1024
        mergeMemMb=4*1024
        hyGenSGEMemMb=4*1024
        hyGenLocalMemMb=1.5*1024

        # extended options
        scanSizeMb = 12
        callRegionsBed = None
        regionStrList = None
        callMemMbOverride = None

        return cleanLocals(locals())


    def validateAndSanitizeOptions(self,options) :

        assertOptionExists(options.runDir,"run directory")
        options.runDir = os.path.abspath(options.runDir)

        workflowScriptPath = os.path.join(options.runDir, options.workflowScriptName)
        if os.path.exists(workflowScriptPath):
            raise OptParseException("Run directory already contains workflow script file '%s'. Each analysis must be configured in a separate directory." % (workflowScriptPath))

        # check reference fasta file exists
        assertOptionExists(options.referenceFasta,"reference fasta file")
        options.referenceFasta=validateFixExistingFileArg(options.referenceFasta,"reference")

        # check for reference fasta index file:
        faiFile=options.referenceFasta + ".fai"
        if not os.path.isfile(faiFile) :
            raise OptParseException("Can't find expected fasta index file: '%s'" % (faiFile))

        # check for bed file of call regions and its index file
        options.callRegionsBed = checkFixTabixIndexedFileOption(options.callRegionsBed, "call-regions bed")

        if (options.regionStrList is None) or (len(options.regionStrList) == 0) :
            options.genomeRegionList = None
        else :
            options.genomeRegionList = [parseGenomeRegion(r) for r in options.regionStrList]

        # validate chromosome names appearing in region tags and callRegions bed file
        if (options.callRegionsBed is not None) or (options.genomeRegionList is not None) :
            refChromInfo = getFastaInfo(options.referenceFasta)
            if options.callRegionsBed is not None :
                for chrom in getTabixChromSet(options.tabixBin, options.callRegionsBed) :
                    if chrom not in refChromInfo :
                        raise OptParseException("Chromosome label '%s', in call regions bed file '%s', not found in reference genome." %
                                                (chrom, options.callRegionsBed))

            if options.genomeRegionList is not None :
                for (genomeRegionIndex, genomeRegion) in enumerate(options.genomeRegionList) :
                    chrom = genomeRegion["chrom"]
                    if chrom not in refChromInfo :
                        raise OptParseException("Chromosome label '%s', parsed from region argument '%s', not found in reference genome." %
                                                (chrom, options.regionStrList[genomeRegionIndex]))
