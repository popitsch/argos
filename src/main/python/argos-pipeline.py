'''
Created on Mar 20, 2014

@author: niko
'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter, \
    ArgumentTypeError
from subprocess import *

import sys, os
from pytz import reference
# Necessary for including python modules from a parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import *
import string

import datetime, time
import logging
import tempfile


usage = '''



    _   __  __  __  __ 
   /_| /__)/_  /  )(   
  /  |/ | (__)(__/__)  pipeline

  Copyright 2014 CIBIV. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.
  
  Calculates ARGOS score files for the passed fasta sequences.
  This pipeline requires access to a number of external execultables (see below)
  as well as some general linux tools (such as "cat").
  The output is a number of BigWig and corresponding CODOC files 
  containing the ARGOS signals.

USAGE
'''

# Required executables 
EXE_NGM="/software/ngm/ngm-0.4.12-pre/ngm"
EXE_NGSTOOLS="ngs-tools-java-32G"
EXE_CODOC="codoc"
EXE_WTBW="wigToBigWig"

def wigToBW(inWig, outBW, fasta, additionalParameters="-clip", override=True):
    success = True
    if(not files_exist(outBW) or override):
        chrSizes = fasta + ".chrSizes"
        created = False
        if not os.path.exists(chrSizes):
            created = True
            success = success and pipelineStep(fasta, chrSizes, [EXE_NGSTOOLS, "mfasta", "stats", "-f", fasta], shell=True, stdout=chrSizes)        
        success = success and pipelineStep(inWig, outBW, [EXE_WTBW + " " + additionalParameters + " " + inWig + " " + chrSizes + " " + outBW], shell=True)
        if created:
            removeFile(chrSizes)      
    return success


def calculateWigs(fastas, SCORES, ISS, AMB, MSD, DATA, dontclean):
 
    for f in fastas:
        tf = tempfile.NamedTemporaryFile(dir=args.tmpdir, delete=False).name
        cix = f + ".cix"
        ngmtmp1 = f + "-enc.ngm"
        ngmtmp2 = f + "-ht-13-2.ngm"
        bam = tf + ".bam"
        params = ""  # additional NGM parameters. Use --end-to-end to do a global alignment, use -g to use the GPU   
        if args.useGpu:
            params="-g"
        os.unlink(tf)
        
        try:
            print("Processing " + f);
            logging.info("Calculating ISS (%s), AMB (%s), MSD (%s), local signals with context-size (%s), and DATA (%s).", ISS, AMB, MSD, str(args.step), DATA);
            
            if not os.path.exists(SCORES) or args.recalcScores:
                if not pipelineStep(f, bam, [EXE_NGSTOOLS, "NGMTools", "mfasta2bam",
                                         "-f", f, "-b", bam, "-c", "all", "-rl", str(args.rl), "-step", str(args.step)]):
                    return False
                
                if not pipelineStep([f, bam], SCORES, 
                                    [EXE_NGM,
                                     "-r", f,
                                     "-q", bam,
                                     "-o", SCORES,
                                     "--argos",
                                     "--argos-min-score", "0.5",
                                     "--gz",
                                     params]):
                    return False
            else:
                print("SKIPPED creation of existing SCORES file " + SCORES);
                logging.info("SKIPPED creation of existing SCORES file (%s).", SCORES);
            
            if not pipelineStep([SCORES], [ISS, AMB, MSD, DATA], ["argos", "compress",
                "-s", SCORES,
                "-rl", str(args.rl),
                "-step", str(args.step),
                "-outWigISS", ISS,
                "-outWigAMB", AMB,
                "-outWigMSD", MSD,
                "-outData", DATA,
                "-ctxSize", str(args.ctxSize),
                "-t", args.tmpdir
                ]):
                return False
        finally:
            if not dontclean:
                print "removing temp files..."
                if os.path.exists(bam): 
                    removeFile([bam])
                if os.path.abspath(cix).startswith(os.path.abspath(args.tmpdir)):
                    removeFile([cix, ngmtmp1, ngmtmp2])                            
    return True

def compressCodoc(wig, codoc):
    if not pipelineStep([wig], [codoc], [EXE_CODOC, "compress", "-cov", wig, "-o", codoc, "-qparam", "0", "-scaleCoverage", "100"]):
        return False
    return True

parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-g", "--genome", type=existing_file, required=True, dest="genome", metavar="genome", help="The mfasta file containing the considered genome")
parser.add_argument("-c", "--context", type=existing_file, required=False, dest="context", metavar="context", help="The mfasta files containing the context-genome(s)", nargs='*')
parser.add_argument("-o", "--outdir", type=str, required=True, dest="outdir", metavar="outdir", help="output directory")
parser.add_argument("-t", "--tmpdir", type=str, required=True, dest="tmpdir", metavar="tmpdir", help="temp directory")
parser.add_argument("-rl", type=int, required=False, dest="rl", metavar="rl", help="read length", default=100)
parser.add_argument("-step", type=int, required=False, dest="step", metavar="step", help="step size", default=5)
parser.add_argument("-ctxSize", type=int, required=False, dest="ctxSize", metavar="ctxSize", help="context size for local signals", default=10000)
parser.add_argument("-dontclean", action="store_true", dest="dontclean", help="if set to true, the temp files will not be removed (for debugging purposes only!)", default=False)
parser.add_argument("-recalcScores", type=bool, required=False, dest="recalcScores", metavar="recalcScores", help="force recalculation of mapping scores", default=False)
parser.add_argument("-calcChrom", type=bool, required=False, dest="calcChrom", metavar="calcChrom", help="Calculation CHR ctx scores", default=False)
parser.add_argument("-gpu", action="store_true", dest="useGpu", help="Use GPU for alignment", default=False)


args = parser.parse_args() 
    
# create directories
if not os.path.exists(os.path.dirname(args.outdir)):
    print("Creating dir " + os.path.dirname(args.outdir))
    os.makedirs(os.path.dirname(args.outdir))
if not os.path.exists(os.path.dirname(args.tmpdir)):
    print("Creating dir " + os.path.dirname(args.tmpdir))
    os.makedirs(os.path.dirname(args.tmpdir))

# start 
log = os.path.join(args.outdir, "genamb-pipeline.log")
logging.basicConfig(filename=log, level=logging.DEBUG)                
logging.info("==========================================================")
logging.info("    _   __  __  __  __             ")
logging.info("   /_| /__)/_  /  )(               ") 
logging.info("  /  |/ | (__)(__/__)  pipeline    ")    
logging.info("                                   ")
logging.info("==========================================================")
logging.info("Started script at %s.", datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
print "genome " + os.path.abspath(args.genome)

contextFile = None
chrFiles = []
try:
    # create context genome
    if args.context is not None:
        contextFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, delete=False).name
        logging.info("Create joined contextFile " + contextFile)
        genomeInCtx = False
        totalsize = 0
        for ctx in args.context:
            print os.path.abspath(ctx)
            if os.path.abspath(ctx) == os.path.abspath(args.genome):
                genomeInCtx = True
            if not pipelineStep(ctx, contextFile, ["cat", ctx], stdout=contextFile, append=True, shell=True):
                raise IOError("Could not concatenate files")
            totalsize += os.stat(ctx).st_size
        if not genomeInCtx:
            if not pipelineStep(args.genome, contextFile, ["cat", args.genome], stdout=contextFile, append=True, shell=True):
                raise IOError("Could not concatenate contextFile files")
            totalsize += os.stat(args.genome).st_size
        if not totalsize == os.stat(contextFile).st_size:
             raise StandardError("Lens do not match!")
         
    # split genome (if there are multiple chroms)
    if args.calcChrom:
        tf = tempfile.NamedTemporaryFile(dir=args.tmpdir, delete=False).name
        if not pipelineStep(args.genome, None, [EXE_NGSTOOLS, "fasta", "split", "-f", args.genome, "-o", tf]):
            raise IOError("Could not split genome")
        num = 1
        while True:
            cf = tf + "." + str(num)
            if not os.path.exists(cf):
                break
            chrFiles += [cf] 
            num += 1
        os.unlink(tf)
        if num == 2:
            print "found only one chrom. do not split..."
            os.unlink(chrFiles[0])
            chrFiles = []       
    
    # calculate WIG files
    prefix = args.outdir + os.path.splitext(os.path.basename(args.genome))[0]
    SCORES = prefix + ".scores.gz"
    ctxDATA = prefix + "-CTX.genAmbData"
    genomeDATA = prefix + "-GENOME.genAmbData"
    chrDATA = prefix + "-CHR.genAmbData"
    
    ctxISSWig = prefix + "-CTX.ISS.wig"  
    ctxAMBWig = prefix + "-CTX.AMB.wig"
    ctxMSDWig = prefix + "-CTX.MSD.wig"
    
    genomeISSWig = prefix + "-GENOME.ISS.wig"
    genomeAMBWig = prefix + "-GENOME.AMB.wig"
    genomeMSDWig = prefix + "-GENOME.MSD.wig"
     
    chrISSWig = prefix + "-CHR.ISS.wig"
    chrAMBWig = prefix + "-CHR.AMB.wig"
    chrMSDWig = prefix + "-CHR.MSD.wig"
    
    CTX = ".ctx" + str(args.ctxSize) + ".wig"
    
    # GENOMIC CONTEXT
    if args.context is not None:
        logging.info("--------------------------------------------------")
        logging.info("CALCULATE GENOMIC CTX SCORES")
        logging.info("--------------------------------------------------")
        if not calculateWigs([contextFile], SCORES, ctxISSWig, ctxAMBWig, ctxMSDWig, ctxDATA, args.dontclean):
            raise IOError("Could not calculate WIGs")    
        # now filter the WIG files and remove all entries not found in the calculated genome.
        pipelineStep(args.genome, None, [EXE_NGSTOOLS, "WigTools", "filterByFasta", "-w", ctxISSWig, "-f", args.genome, "-o", ctxISSWig + ".filtered"])
        pipelineStep(args.genome, None, [EXE_NGSTOOLS, "WigTools", "filterByFasta", "-w", ctxAMBWig, "-f", args.genome, "-o", ctxAMBWig + ".filtered"])
        pipelineStep(args.genome, None, [EXE_NGSTOOLS, "WigTools", "filterByFasta", "-w", ctxMSDWig, "-f", args.genome, "-o", ctxMSDWig + ".filtered"])
        pipelineStep(args.genome, None, [EXE_NGSTOOLS, "WigTools", "filterByFasta", "-w", ctxISSWig + CTX, "-f", args.genome, "-o", ctxISSWig + CTX + ".filtered"])
        pipelineStep(args.genome, None, [EXE_NGSTOOLS, "WigTools", "filterByFasta", "-w", ctxAMBWig + CTX, "-f", args.genome, "-o", ctxAMBWig + CTX + ".filtered"])
        pipelineStep(args.genome, None, [EXE_NGSTOOLS, "WigTools", "filterByFasta", "-w", ctxMSDWig + CTX, "-f", args.genome, "-o", ctxMSDWig + CTX + ".filtered"])

        wigToBW(ctxISSWig + ".filtered", ctxISSWig + ".bw", args.genome)
        wigToBW(ctxAMBWig + ".filtered", ctxAMBWig + ".bw", args.genome)
        wigToBW(ctxMSDWig + ".filtered", ctxMSDWig + ".bw", args.genome)        
        wigToBW(ctxISSWig + CTX + ".filtered", ctxISSWig + CTX + ".bw", args.genome)
        wigToBW(ctxAMBWig + CTX + ".filtered", ctxAMBWig + CTX + ".bw", args.genome)
        wigToBW(ctxMSDWig + CTX + ".filtered", ctxMSDWig + CTX + ".bw", args.genome)

        compressCodoc(ctxISSWig + ".filtered", ctxISSWig + ".codoc")
        compressCodoc(ctxAMBWig + ".filtered", ctxAMBWig + ".codoc")
        compressCodoc(ctxMSDWig + ".filtered", ctxMSDWig + ".codoc")
        compressCodoc(ctxISSWig + CTX + ".filtered", ctxISSWig + CTX + ".codoc")
        compressCodoc(ctxAMBWig + CTX + ".filtered", ctxAMBWig + CTX + ".codoc")
        compressCodoc(ctxMSDWig + CTX + ".filtered", ctxMSDWig + CTX + ".codoc")

        removeFile([ctxISSWig, ctxISSWig + ".filtered", ctxAMBWig, ctxAMBWig + ".filtered", ctxMSDWig, ctxMSDWig + ".filtered"])
        removeFile([ctxISSWig + CTX, ctxISSWig + CTX + ".filtered", ctxAMBWig + CTX, ctxAMBWig + CTX + ".filtered", ctxMSDWig + CTX, ctxMSDWig + CTX + ".filtered"])

    # GENOME
    logging.info("--------------------------------------------------")
    logging.info("CALCULATE GENOMIC SCORES")
    logging.info("--------------------------------------------------")
    if not calculateWigs([args.genome], SCORES, genomeISSWig, genomeAMBWig, genomeMSDWig, genomeDATA, args.dontclean):
        raise IOError("Could not calculate WIGs")    
    wigToBW(genomeISSWig, genomeISSWig + ".bw", args.genome)
    wigToBW(genomeAMBWig, genomeAMBWig + ".bw", args.genome)
    wigToBW(genomeMSDWig, genomeMSDWig + ".bw", args.genome)
    wigToBW(genomeISSWig + CTX, genomeISSWig + CTX + ".bw", args.genome)
    wigToBW(genomeAMBWig + CTX, genomeAMBWig + CTX + ".bw", args.genome)
    wigToBW(genomeMSDWig + CTX, genomeMSDWig + CTX + ".bw", args.genome)
    compressCodoc(genomeISSWig, genomeISSWig + ".codoc")
    compressCodoc(genomeAMBWig, genomeAMBWig + ".codoc")
    compressCodoc(genomeMSDWig, genomeMSDWig + ".codoc")
    compressCodoc(genomeISSWig + CTX, genomeISSWig + CTX + ".codoc")
    compressCodoc(genomeAMBWig + CTX, genomeAMBWig + CTX + ".codoc")
    compressCodoc(genomeMSDWig + CTX, genomeMSDWig + CTX + ".codoc")

    removeFile([genomeISSWig, genomeAMBWig, genomeMSDWig])
    removeFile([genomeISSWig + CTX, genomeAMBWig + CTX, genomeMSDWig + CTX])

    # CHROMS
    if args.calcChrom and num > 2:
        logging.info("--------------------------------------------------")
        logging.info("CALCULATE CHR SCORES")
        logging.info("--------------------------------------------------")
        if not calculateWigs(chrFiles, SCORES, chrISSWig, chrAMBWig, chrMSDWig, chrDATA, args.dontclean):
            raise IOError("Could not calculate WIGs")       
        wigToBW(chrISSWig, chrISSWig + ".bw", args.genome)
        wigToBW(chrAMBWig, chrAMBWig + ".bw", args.genome)
        wigToBW(chrMSDWig, chrMSDWig + ".bw", args.genome)
        wigToBW(chrISSWig + CTX, chrISSWig + CTX + ".bw", args.genome)
        wigToBW(chrAMBWig + CTX, chrAMBWig + CTX + ".bw", args.genome)
        wigToBW(chrMSDWig + CTX, chrMSDWig + CTX + ".bw", args.genome)
        compressCodoc(chrISSWig, chrISSWig + ".codoc")
        compressCodoc(chrAMBWig, chrAMBWig + ".codoc")
        compressCodoc(chrMSDWig, chrMSDWig + ".codoc")
        compressCodoc(chrISSWig + CTX, chrISSWig + CTX + ".codoc")
        compressCodoc(chrAMBWig + CTX, chrAMBWig + CTX + ".codoc")
        compressCodoc(chrMSDWig + CTX, chrMSDWig + CTX + ".codoc")
        removeFile([chrISSWig, chrAMBWig, chrMSDWig])
        removeFile([chrISSWig + CTX, chrAMBWig + CTX, chrMSDWig + CTX])
    
    
    print "Finished." 
    logging.info("==========================================================")
    logging.info("FINISHED script at %s.", datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
    logging.info("==========================================================")
    
except Exception as e:
    print "Error!"
    raise 
finally:
    if not args.dontclean:
        print "CLEANING " + str(args.dontclean)
        if contextFile is not None:
            print "remove " + contextFile
            os.unlink(contextFile)
        for f in chrFiles:
            print "remove " + f
            os.unlink(f)
        
logging.info("Finished.")
                
