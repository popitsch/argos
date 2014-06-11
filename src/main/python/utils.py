'''
Created on Mar 20, 2013

@author: niko
'''

if __name__ == '__main__':
    pass

import subprocess
from subprocess import *
from argparse import ArgumentTypeError
import os
import sys
import logging
import signal


# FIXME: there is a problem with called subprocessed that do not terminate! This will hang this code!
# example-cmd:
# 'bowtie2', '-t', '10', '-x', '/project/ngs-work/meta/reference/genomes/hg19_human/bowtie2_path/hg19', '-1', '../FASTQ/Probe7_CAGATC_L004_R1_001.fastq.gz', '-2', '../FASTQ/Probe7_CAGATC_L004_R2_001.fastq.gz', '-S', '../FASTQ/Probe7_CAGATC_L004_R1_001.fastq.gz.bt2.sam'
# This will result in a non-terminating PERL job. IN turn, the python program also stalls. 
# run a commandline task
def runTask(cmd, shell=False):
    logging.info(cmd)
    if shell:
        out = check_output(" ".join(cmd), shell=True, stderr=subprocess.STDOUT)
    else: 
        out = check_output(cmd, stderr=subprocess.STDOUT)
    return out
        
def files_exist(files):
    if (type(files) is list) :
        for f in files:
            if not os.path.exists(f):
                return False
    else:
        if not os.path.exists(files):
            return False
    return True

# use in argument parser
def existing_file(files):
    if files_exist(files):
        return files
    raise ArgumentTypeError("Not all files exist")

# remove a (list of) file(s) (if it/they exists)
def removeFile(files):
    if (type(files) is list) :
        for f in files:
            if os.path.exists(f):
                os.remove(f)
    else:
        if os.path.exists(files):
            os.remove(files)

# check whether a file exists and exit if not        
def checkFile(files):
    if (type(files) is list) :
        for f in files:
            checkFile(f)
    else :
        if not os.path.exists(files):
            print "Error: file", files, "was not found! Exiting..."
            sys.exit(1)       
                
def pipelineStep(inputfile, outFile, cmd, shell=False, stdout=None, append=False):
    try:
        if inputfile is not None:
            if (type(inputfile) is list):
                for f in inputfile:
                    checkFile(f)
            else:
                checkFile(inputfile) 
        
        if stdout is not None:
            if shell is False:
                raise ArgumentTypeError("When using the parameter stdout, the shell parameter must be set to True!")
            else: 
                if append:
                    cmd.append(">>")
                else:
                    cmd.append(">")
                cmd.append(stdout)                                                                                                                                                                                                                                                                 

                
        # Overwrite output file?!
        if outFile is not None:
            if (type(outFile) is list) :
                for f in outFile:
                    if files_exist(f):
                        logging.warn("Overwriting file %s", f)
                else:
                    if files_exist(outFile):
                        logging.warn("Overwriting file %s", outFile)

        out = runTask(cmd, shell)         

        logging.debug(out)        
        if outFile is not None:
            checkFile(outFile)
        return True
    except CalledProcessError as e:
        logging.debug(e.output)
        logging.error("ERROR %s - removing outputfile %s", e, outFile)
        if outFile is not None:
            if (type(outFile) is list) :
                for f in outFile:
                    removeFile([f]) 
            else:
                removeFile([outFile]) 
        return False        

# sample <number> of reads randomly from inFile (.fa/.fq) and writes to outFile. <number> can either be a absolute number (e.g. 1 000 000) or a percentage (e.g. 0.1 for 10 %)
def sampleFromFastX(inFile, outFile, number, seed=100):     
    return pipelineStep(inFile, outFile, ["seqtk", "sample", "-s", str(seed), inFile, str(number)], shell=True, stdout=outFile)


#Replaces the file extension of inFile to with <newExtension> and adds a suffix
#Example replaceExtension("reads.fq", ".sam", suffix="_namg") => reads_ngm.sam
def replaceExtension(inFile, newExtension, suffix=""):    
    return os.path.splitext(inFile)[0] + suffix + newExtension        
        
#Converts a bam file to fq using bam2fastq
def bam2fq(inFile, outFile, delinFile=False, override=False):
    success = True
    
    if((not files_exist(outFile[::-1].replace("#"[::-1], "_1"[::-1], 1)[::-1]) and not files_exist(outFile[::-1].replace("#", "", 1)[::-1])) or override):
        success = pipelineStep(inFile, outFile, ["bam2fastq", "-f --aligned --unaligned --no-filtered", "-o", outFile, inFile], shell=True)
        if(success and delinFile):
            removeFile(inFile)
    else:
        logging.info("bam2fq skipped. File " + outFile[::-1].replace("#"[::-1], "_1"[::-1], 1)[::-1] + " already exists.")
    return success

def runFASTQC(inFile, outFolder, threads=1, additionalParameters="", override=False):
    success = True
    if ( not files_exist(outFolder)):
        os.makedirs(outFolder)
    outFile = os.path.join(outFolder, os.path.basename(inFile).split(".")[0] + "_fastqc.zip" )
    if(not files_exist(outFile) or override):
        success = pipelineStep(inFile, None, ["fastqc", "-t", str(threads), "-o", outFolder, additionalParameters, inFile], shell=True)
    return success

def indexBam(inFileBam, override=False):
    success = True
    idxFile = inFileBam + ".bai"
    if(not files_exist(idxFile) or override):
        success = success and pipelineStep(inFileBam, idxFile, ["samtools", "index", inFileBam], shell=True)
    return success
    
def sam2bam(inFile, outFile, index=True, sort=True, delinFile=False, override=False, onlyUnique=False, onlyProperPaired=False, filterMQ=0):
    if(onlyUnique and filterMQ == 0):
        filterMQ = 1;
        
    success = True    
    if(not files_exist(outFile) or override):        
        cmd = ["samtools", "view", "-Sb", inFile, "-o", outFile]
        if filterMQ > 0:
            cmd+=["-q", str(filterMQ)]
        if onlyProperPaired:
            cmd+=["-f", "2"]
        success = success and pipelineStep(inFile, outFile, cmd, shell=True)
        
        if(sort):         
            tmp = outFile + "_tmp"
            os.rename(outFile, tmp)                      
            success = success and pipelineStep(tmp, outFile, ["samtools", "sort", tmp, replaceExtension(outFile, "")], shell=True)
            if(success):
                removeFile(tmp)
        if(index):
            success = success and indexBam(outFile)
            
        if(success and delinFile):
            removeFile(inFile)
    else:
        logging.info("Skipping sam2bam. " + outFile + " already exists.");
    return success

def flagstat( bam ):
    success = True
    flagstat = bam + ".flagstat"
    if not files_exist( flagstat):
        success = success and pipelineStep(None, flagstat, [ "samtools", "flagstat", bam], shell=True, stdout=flagstat)
    return success

def mapqHist( bam ):
    success = True
    mapqhist = bam + ".mapqhist"
    if not files_exist( mapqhist ):
        success = success and pipelineStep(bam, mapqhist, ["ngs-tools-java", "QualityDistribution", "printMAPQHistogram", "-r", bam, "-o", mapqhist ])
    return success

