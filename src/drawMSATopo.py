#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Filename:  drawMSATopo.py
#
# Description:
#   Visualize aligned membrane proteins by highlighting features of membrane
#   topology
#
# Author:
#   Nanjiang Shu  nanjiang.shu@scilifelab.se

import string
import sys
import re
import os
import myfunc
import math
import libtopologycmp as lcmp
import numpy as np
import Bio.SubsMat.MatrixInfo
import subprocess
from matplotlib.lines import *
from colour import Color

alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

rundir = os.path.dirname(os.path.realpath(__file__))
# pySVG
# import pysvg.structure
# import pysvg.builders
# import pysvg.text

# matplotlib
# from matplotlib.font_manager import FontProperties
# from pylab import *
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

import tempfile

import logging
import logging.config
import yaml

logger = logging.getLogger(__name__)

# PyX
#import pyx


nodename = os.uname()[1]
colorList = ["red", "blue", "green", "cyan","pink"]
colorList_DG_profile = ["red", "black", "green", "cyan", "yellow"]

ylabel_DGprofile =  "\u0394G (kcal/mol)"

PIL_user_path = os.environ['HOME'] + "/usr/lib64/python2.6/site-packages/PIL"
if nodename.find("uppmax") != -1:
    sys.path.append(PIL_user_path)

from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont

GAP = myfunc.GAP

#import cProfile

# ChangeLog 2011-11-21 #{{{
#    1. Bug solved for height overflow. When sepration line is added, canvas
#    height should be extended as well.
#    2. AA sequences are obtained from the original sequence source, not the
#    alignment. In this case, un-aligned sequences will be matched to aligned
#    topology accordingly.
# ChangeLog 2011-11-23
#    1. columns with mostly gaps and not TM regions are shrinked according to
#    its reciprocal gap percentage.
# ChangeLog 2011-11-24
#    1. Bug solved for ShrinkGapInMSA. The function ShrinkGapInMSA_0 is used by
#    default now.
# ChangeLog 2011-11-25
#    1. Bug solved in IsSafetoDeleteTheRegion so that the beginning and end
#    topology state will not be deleted.
#    2. For blocks with no 'M' state but with both 'i' and 'o', if for each
#    subsequence in the block there are at most one state, it is OK to shrink
#    the block to one column
# ChangeLog 2011-12-01 
#    1. Bug in ShrinkGapInMSA_0 fixed
# ChangeLog 2012-02-22
#    1. For pairwise topology alignment, annotation is the sequence id
# ChangeLog 2012-02-23
#    1. Draw DG profile
# ChangeLog 2012-03-21 
#    when aaseq file is added, only use this file 
# ChangeLog 2012-09-28
#    Location of dgscanProg and font_dir set to user folder. Using environ
#    DATADIR3
# ChangeLog 2013-04-17
#    Add krbias to outfile name when -krbias is enabled
# ChangeLog 2013-11-20
#    Add option -colorTMbox, so that the whole TM helix will be colored red,
#    including the gaps within the TM helix
#    Add option -showTMidx, so that the index of TM helices will be showin
#    instead of the sequence.
# ChangeLog 2015-03-17
#    1. for matplotlib, set the font by fontpath, avoid "font not found
#       problem" if not in the system
#    2. the shrinkrate is set so that the length of the normalized sequence
#       length is 100, unless set globally
#    3. set also MAX_SIZE_ANNOTATION
# ChangeLog 2016-09-13 
#    1. draw signal peptide as green, now only in PIL 
# ChangeLog 2017-06-29
#    Added a new function to draw core-rainbow
# ChangeLog 2017-07-03
#    Moved the fond_dir fonts within the script folder, removed the DATADIR3
#    env variable
# ChangeLog 2017-08-07 
#   1. add the flag "-ptag"
#}}}

# global constant


usage="""
Usage:   drawMSATopo.py [-i] topomsa-in-fasta-format
Options:
  -method    STR   Modules to use for plotting, (default: pil)
                   Can be pyx, svg, pil, mat, core-rainbow
  -of        STR   Output format, can be png
  -l        FILE   Set input file list
  -mode      STR   Image mode, P or RGB, (default: P)
  -fontsize  INT   Set the font size, (default: 9)
  -text  y|n       Wether draw text i, o or M in the alignment, (default:  yes)
                   if no, then only the background color is shown
                   red for M, faded-yellow for i and faded-blue for o
  -sep  y|n        Whether draw a seprator line in grey between each group
                   (default: yes)
  -pfm  y|n        Whether draw profile for 'M' state, (default: yes)
  -pdg  y|n        Whether draw DG profile, (default: no)
  -pmsa y|n        Whether draw MSA region, (default: yes)
  -pscale y|n      Whether draw the residue position scale bar, (default: yes)
                   If MSA region is not drawn, this will also be disabled.
  -ptag y|n        Whether draw a vertical bar for proteins in different groups, (default: no)
  -dgpfile FILE    DG profile file produced by myscanDG.pl
  -aapath  DIR     Set path for amino acid sequence file, if set, sequence file
                   will be searched as $ID.fa
  -outpath DIR     Set outpath, (default: $dirname(infile))
  -autosize y|n    Whether autosize font, (default: yes)
  -shrink   y|n    Whether shrink gap regions, (default: yes)
  -m-shrink INT    method of shrinking, (default: 1)
                   0: shrink both non-TM region and TM region
                   1: just shrink non-TM region
  -aaseq   FILE    Set aaseq file for all 
  -krbias          Draw krbias
  -maxdistkr INT   Maximal distance to TM for KR residue, (default: 100)
  -win       INT   Window size for text and html format alignment output. (default: 70)
  -h, --help       Print this help message and exit
  -htmlheader STR  Set header text for HTML output
  -colorhtml       Use colorful output for HTML alignment
  -colorTMbox      Color the whole TM helix as a red rectangle, even gaps
  -colorkingdom    Color the TM regions by kingdom (Archaea green, Bacteria Red, Eukaryota blue) as present in annotations.
  -advtopo         Show advanced topology at top (reentrant regions and in/out helices)
  -showTMidx       Display index of TM helix as the text for the sequence, e.g. TM1 TM2
  -shrinkrate FLOAT     Proportional shrink rate, (default: 1.0)
  -shrinkrateTM FLOAT   Proportional shrink rate for TM regions, (default: 2.0)
  -max-hold-loop INT    Maximal positions to keep for loop regions (default: 12)
  -imagescale FLOAT     Overal scale of the image (default: None). If not set, it will be calculated automatically
  -h2wratio FLOAT       Set the height to width ratio (default: None). If not set, it use the original ratio
  -cleanplot       Make clean plot, works only for the PIL mode
  -showgap         Show gap as blank region, works only for the PIL mode
  -debug                Print debug information, (default: no)

Created 2011-09-05, updated 2020-06-26, Nanjiang Shu

Examples:
# Draw topology alignment with aligned and unaligned regions, aligned regions
# shown in rainbow color

%s -pdg n -shrink yes -method core-rainbow -text no -shrink yes topoalignFile


"""%(sys.argv[0])
def PrintHelp():
    print(usage)

def sig2(x, scale=1):
    """ Calculate sigmoid value
    """
    return 1/(1+math.exp(-x/scale))


def WriteTXTAlignment(idList, annoList, alignedTopoSeqList,#{{{
        originalAlignedTopoSeqList, aaSeqList, final2seq_idxMapList,
        outfile):
    WIDTH = g_params['window_size']
    maxSizeAnno = max([len(x) for x in annoList])
    lengthAlignment = len(alignedTopoSeqList[0])
    numSeq = len(idList)
    posTMList = [myfunc.GetTMPosition(x) for x in alignedTopoSeqList]

    fpout = open(outfile, "w")

    strs = [""]*numSeq
    j = 0 # iterator for the alignment position
    isStart = True
    cnt = 0
    while j < lengthAlignment:
        if isStart:
            strs = [""]*numSeq
            for i in range(numSeq):
                try:
                    strs[i] += "%-*s %4d "%(maxSizeAnno, annoList[i],
                            final2seq_idxMapList[i][j])
                except KeyError:
                    print("final2seq_idxMapList error  i=%d, j=%d"%(i,j))
                    pass
            isStart = False
        isWithinTMregion = False
        for i in range(numSeq):
            if lcmp.IsWithinTMRegion(j, posTMList[i]):
                aa = aaSeqList[i][j].upper()
                isWithinTMregion = True # if hit TM region of any sequence, set as TRUE
            else:
                aa = aaSeqList[i][j].lower()
            strs[i] += aa
        if (cnt >= WIDTH and isWithinTMregion == False) or (j >= lengthAlignment-1):
            for i in range(numSeq):
                strs[i] += " %4d"%(final2seq_idxMapList[i][j])
            for i in range(numSeq):
                fpout.write("%s\n"%(strs[i]))
            fpout.write("\n")
            fpout.write("\n")
            strs = [""]*numSeq
            isStart = True
            cnt = 0
        j += 1
        cnt += 1
    fpout.close()
    #}}}
def WriteHTMLAlignment(idList, annoList, alignedTopoSeqList,#{{{
        originalAlignedTopoSeqList, aaSeqList, final2seq_idxMapList,
        outfile):
    # but do not break in the middle of a helix
    WIDTH = g_params['window_size']
    maxSizeAnno = max([len(x) for x in annoList])
    lengthAlignment = len(alignedTopoSeqList[0])
    numSeq = len(idList)
    posTMList = [myfunc.GetTMPosition(x) for x in alignedTopoSeqList]

    fpout = open(outfile, "w")
    header = """
<!DOCTYPE html>
<html>
<body>
<table border="0" cellspacing="0" cellpadding="0">
"""
    tail = """
</table>
</body>
</html>
    
"""
    print(header, file=fpout)


    strs = [""]*numSeq
    j = 0 # iterator for the alignment position
    isStart = True
    cnt = 0
    while j < lengthAlignment:
        if isStart:
            strs = [""]*numSeq
            for i in range(numSeq):
                strs[i] += "<tr><td>%s</td><td>%d</td>"%(annoList[i],
                        final2seq_idxMapList[i][j])
            isStart = False
        isWithinTMregion = False
        for i in range(numSeq):
            if lcmp.IsWithinTMRegion(j, posTMList[i]):
                aa = aaSeqList[i][j].upper()
                isWithinTMregion = True # if hit TM region of any sequence, set as TRUE
                strs[i] += "<td><b><font color=\"black\">%s</font></b></td>"%(aa)
            else:
                aa = aaSeqList[i][j].lower()
                strs[i] += "<td><font color=\"grey\">%s</font></td>"%(aa)
        if (cnt >= WIDTH and isWithinTMregion == False) or (j >= lengthAlignment-1):
            for i in range(numSeq):
                strs[i] += "<td>%d</td></tr>"%(final2seq_idxMapList[i][j])
            for i in range(numSeq):
                fpout.write("%s\n"%(strs[i]))
            fpout.write("\n")
            strs = [""]*numSeq
            isStart = True
            cnt = 0
        j += 1
        cnt += 1

    print(tail, file=fpout)

    fpout.close()
    #}}}
def WriteHTMLAlignment2(idList, annoList, alignedTopoSeqList,#{{{
        originalAlignedTopoSeqList, aaSeqList, final2seq_idxMapList,
        outfile):
    logger = logging.getLogger(__name__)
# for two sequence pairwise alignment
# assign also the identical and similarity by using BLOSUM62
    annoList = idList
    #WIDTH = 90 # but do not break in the middle of a helix, adjust 1
    #WIDTH = 60 # but do not break in the middle of a helix
    WIDTH = g_params['window_size']

    maxSizeAnno = max([len(x) for x in annoList])
    lengthAlignment = len(alignedTopoSeqList[0])
    numSeq = len(idList)
    posTMList = [myfunc.GetTMPosition(x) for x in alignedTopoSeqList]

    blosum62 = Bio.SubsMat.MatrixInfo.blosum62

    if g_params['colorhtml']:
        color_TM = 'red'
        color_nonTM = 'grey'
    else:
        color_TM = 'black'
        color_nonTM = 'grey'

    fpout = open(outfile, "w")
    header = """
<!DOCTYPE html>
<html>
<body>
<h3>%s</h3>
<pre>
"""%(g_params['htmlheader'])
    tail = """
</pre>
</body>
</html>
    
"""
    print(header, file=fpout)



    strs = [""]*numSeq
    j = 0 # iterator for the alignment position
    isStart = True
    cnt = 0
    while j < lengthAlignment:
        if isStart:
            strs = [""]*(numSeq+1)
            for i in range(numSeq):
                try:
                    strs[i] += "%-*s %4d "%(maxSizeAnno, annoList[i],
                            final2seq_idxMapList[i][j])
                except KeyError:
                    logger.debug("final2seq_idxMapList error  i=%d, j=%d"%(i,j))
                    pass
            strs[2] += "%-*s %4s "%(maxSizeAnno, "", "")
            isStart = False
        isWithinTMregion = False

        aa1 = aaSeqList[0][j].upper()
        aa2 = aaSeqList[1][j].upper()
        if aa1 == GAP or aa2 == GAP:
            char_rel = " "
        else:
            if (aa1,aa2) in blosum62:
                blosum_score = blosum62[(aa1,aa2)]
            elif (aa2,aa1) in blosum62:
                blosum_score = blosum62[(aa2,aa1)]
            else:
                blosum_score = -1

            if aa1 == aa2:
                char_rel =  "|"
            elif blosum_score > 0:
                char_rel = "."
            else:
                char_rel = " "
        strs[2] += char_rel

        for i in range(numSeq):
            if lcmp.IsWithinTMRegion(j, posTMList[i]):
                aa = aaSeqList[i][j].upper()
                isWithinTMregion = True # if hit TM region of any sequence, set as TRUE
                strs[i] += "<b><font color=\"%s\">%s</font></b>"%(color_TM, aa)
            else:
                aa = aaSeqList[i][j].lower()
                strs[i] += "<font color=\"%s\">%s</font>"%(color_nonTM, aa)
        if ((cnt >= WIDTH and isWithinTMregion == False) 
                or (j >= lengthAlignment-1)
                or j == 190):
            for i in range(numSeq):
                strs[i] += " %4d"%(final2seq_idxMapList[i][j])

            fpout.write("%s\n"%(strs[0]))
            fpout.write("%s\n"%(strs[2])) #relationship
            fpout.write("%s\n"%(strs[1]))
            fpout.write("\n\n")

            strs = [""]*(numSeq+1)
            isStart = True
            cnt = 0
        j += 1
        cnt += 1

    print(tail, file=fpout)

    fpout.close()
    #}}}
def WriteHTMLAlignment3(idList, annoList, alignedTopoSeqList,#{{{
        originalAlignedTopoSeqList, aaSeqList, final2seq_idxMapList,
        outfile):
    logger = logging.getLogger(__name__)
    annoList = idList
    WIDTH = g_params['window_size']
    maxSizeAnno = max([len(x) for x in annoList])
    lengthAlignment = len(alignedTopoSeqList[0])
    numSeq = len(idList)
    posTMList = [myfunc.GetTMPosition(x) for x in alignedTopoSeqList]

    fpout = open(outfile, "w")
    header = """
<!DOCTYPE html>
<html>
<body>
<pre>
"""
    tail = """
</pre>
</body>
</html>
    
"""
    print(header, file=fpout)


    strs = [""]*numSeq
    j = 0 # iterator for the alignment position
    isStart = True
    cnt = 0
    while j < lengthAlignment:
        if isStart:
            strs = [""]*numSeq
            for i in range(numSeq):
                try:
                    strs[i] += "%-*s %4d "%(maxSizeAnno, annoList[i],
                            final2seq_idxMapList[i][j])
                except KeyError:
                    logger.debug( "final2seq_idxMapList error  i=%d, j=%d"%(i,j))
                    pass
            isStart = False
        isWithinTMregion = False
        for i in range(numSeq):
            if lcmp.IsWithinTMRegion(j, posTMList[i]):
                aa = aaSeqList[i][j].upper()
                isWithinTMregion = True # if hit TM region of any sequence, set as TRUE
                strs[i] += "<b><font color=\"black\">%s</font></b>"%(aa)
            else:
                aa = aaSeqList[i][j].lower()
                strs[i] += "<font color=\"grey\">%s</font>"%(aa)
        #print "isWithinTMregion=", isWithinTMregion
        if ((cnt >= WIDTH and isWithinTMregion == False) 
                or (j >= lengthAlignment-1)
                or j == 190):
            for i in range(numSeq):
                strs[i] += " %4d"%(final2seq_idxMapList[i][j])
            for i in range(numSeq):
                fpout.write("%s\n"%(strs[i]))
            fpout.write("\n")
            strs = [""]*numSeq
            isStart = True
            cnt = 0
        j += 1
        cnt += 1

    print(tail, file=fpout)

    fpout.close()
    #}}}
def GetAAPath(topomsafile):#{{{
    """
    Get the path for amino acid sequences
    """
    if g_params['aapath'] != "":
        return g_params['aapath']
    else:
        return myfunc.my_dirname(topomsafile)
#}}}
def GetDGProfileFileName(inFile, seqID):# {{{
    """
    Auto determine the DG profile file name based on the inFile
    """
    DGProfileFile = ""
    dirname_infile = myfunc.my_dirname(inFile)
    rootname_infile = os.path.basename(os.path.splitext(inFile)[0])
    if not os.path.exists(g_params['DGProfileFile']):
        DGProfileFile = dirname_infile + os.sep + rootname_infile + "_dg.txt"
        if not os.path.exists(DGProfileFile):
            DGProfileFile = dirname_infile + os.sep + rootname_infile + "-%s"%(seqID) + "_dg.txt"
            if not os.path.exists(DGProfileFile):
                logger.debug("DGProfileFile not found for inFile %s"%(inFile))
                DGProfileFile = ""
    else:
        DGProfileFile = g_params['DGProfileFile']

    logger.debug("In function GetDGProfileFileName: seqID=%s, DGProfileFile=%s"%(seqID, DGProfileFile))
    return DGProfileFile
# }}}

def GetAASeqDict(topomsafile):#{{{
    """
    Get the amino acid sequence dictionary, keys are seqids
    """
    #if (not g_params['isDrawText']) or g_params['isShrink']:
    #if (not g_params['isDrawText']):
    #    return {}
    #else:
    logger = logging.getLogger(__name__)
    if g_params['aaSeqDict'] != {}:
        return g_params['aaSeqDict']
    else:
        aaSeqDict = {}
        aapath = GetAAPath(topomsafile)
        fastaID = os.path.basename(topomsafile).split('.')[0]
        if topomsafile.find('homology') >= 0:
            fastaAASeqFile = aapath + os.sep + fastaID + '.homology.fa'
        else:
            fastaAASeqFile = aapath + os.sep + fastaID + '.fa'

        if os.path.exists(fastaAASeqFile):
            logger.info("Seqfile %s found"%fastaAASeqFile)
            (aaSeqIDList, aaSeqList) = myfunc.ReadFasta_without_annotation(
                    fastaAASeqFile)
            if len(aaSeqList) <= 0:
                msg = "Failed to read aaSeqFile %s"
                logger.debug(msg)
            else:
                for i in range (len(aaSeqIDList)):
                    aaSeqDict[aaSeqIDList[i]] = aaSeqList[i]
        else:
            logger.debug("aaSeqFile %s does not exist."%(fastaAASeqFile))
        return aaSeqDict
#}}}
def HideNonKRResidue(aaseq):#{{{
    newseq = ""
    for aa in aaseq:
        if  aa in ["K","R"]:
            newseq += aa
        else:
            newseq += " "
    return newseq
#}}}
def GetKRStateFraction(alignedSeqList):#{{{
    """return (cnt_K, cnt_R,  per_K, per_R)"""
    lengthAlignment=len(alignedSeqList[0])
    numSeq = len(alignedSeqList)
    cnt_K = [0]*lengthAlignment
    cnt_R = [0]*lengthAlignment
    for i in range(numSeq):
        alignedSeq = alignedSeqList[i]
        for j in range(lengthAlignment):
            s = alignedSeq[j]
            if s == 'K':
                cnt_K[j] += 1
            elif s == 'R':
                cnt_R[j] += 1

    per_K = [0.0]*lengthAlignment
    per_R = [0.0]*lengthAlignment
    numSeq_float = float(numSeq)
    for j in range(lengthAlignment):
        per_K[j] = cnt_K[j]/(numSeq_float)
        per_R[j] = cnt_R[j]/(numSeq_float)
    return (cnt_K, cnt_R, per_K, per_R)

#}}}
def MatchToAlignedSeq(unalignedseq, alignedseq, seqID): #{{{
    """match the unaligned seq to the aligned seq, gaps are added
    return alignedseq at failure"""
    logger = logging.getLogger(__name__)
    newseq = ""
    j = 0
    GAP = g_params['GAP']
    for i in range(len(alignedseq)):
        if alignedseq[i] != GAP:
            newseq += unalignedseq[j]
            j += 1
        else:
            newseq += GAP
    if len(newseq) != len(alignedseq):
        logger.debug("failed to match sequence for ID %s" %seqID)
        return alignedseq
    else:
        return newseq
#}}}
def ReadInDGProfile(infile):#{{{
    """Read in DG profile output by myscanDG.pl"""
    dgpDict = {}
    try:
        fpin = open(infile, 'r')
        buff = fpin.read()
        fpin.close()
        lines = buff.split('\n')
        numLine = len(lines)
        i = 0
        seqid = ''
        while i < numLine:
            line = lines[i]
            if line.find("#SeqID") == 0:
                seqid = line.split()[1]
            elif line.find("#Number of sliding windows") == 0:
                numWin = int(line.split(':')[1])
                dgp = []
                for j in range(numWin):
                    strs = lines[i+j+1].split()
                    if len(strs) != 2:
                        logger.debug("dgscan file error. strs=%s"%(str(strs)))
                        sys.exit(1)
                    dgp.append((int(strs[0]), float(strs[1])))
                i += numWin
                dgpDict[seqid] = dgp
            i += 1
        #print (dgpDict)
        return dgpDict
    except IOError:
        print("Failed to read dgprofile", infile, file=sys.stderr)
#}}}

def MatchAlignedDGP(dgp, idxmap_aligne2seq, posindexmap, aligned_toposeq):#{{{
    """
    match dgp (a list of tuples) to the aligned toposeq
    posindexmap is the index map from the shrinked seq to the original seq
    idxmap_aligne2seq is a dictionary of index map from the original (non-shrinked) MSA to the gapless seq
    """
    aligned_dgp = []
    lenAlignedSeq = len(aligned_toposeq)
    resMap = {}
    inew = 0
    if len(posindexmap) == 0:
        isShrink = False
    else:
        isShrink = True

    # convert dgp in to dictionary
    dgp_dt = {}
    for (idx, dg) in dgp:
        dgp_dt[idx] = dg

    for j in range(lenAlignedSeq):
        if aligned_toposeq[j] != '-':
            if isShrink:
                j_origseq = idxmap_aligne2seq[posindexmap[j]]
            else:
                j_origseq = idxmap_aligne2seq[j]
            try:
                dg = dgp_dt[j_origseq]
                aligned_dgp.append((j, dg))
            except KeyError:
                pass
    return aligned_dgp
#}}}


def GetFontDimension(font_size):#{{{
    if font_size == 3:
        return (2,4)
    elif font_size == 4:
        return (3,5)
    elif font_size == 5:
        return (4,6)
    elif font_size == 6:
        return (5,7)
    elif font_size == 7:
        return (5,8)
    elif font_size == 8:
        return (6,9)
    elif font_size == 9:
        return (7,10)
    elif font_size == 10:
        return (7,11)
    elif font_size == 11:
        return (8,12)
    elif font_size == 12:
        return (8,13)
    elif font_size == 13:
        return (9,14)
    elif font_size == 14:
        return (10,15)
    else :
         return (8,13)
#}}}

def AutoSizeFontHistogram(ylabel, yticList, widthBox, heigthBox, #{{{
        spaceToLeftBorder):
    maxSizeYtic = max([len(str(x)) for x in yticList])
    numCharSpaceRequied = len(ylabel) + maxSizeYtic*2  + 1
    maxAllowedFontWidth = spaceToLeftBorder / numCharSpaceRequied
    maxAllowdFontHeight = heigthBox / (len(yticList)-1)

#     print "spaceToLeftBorder=",spaceToLeftBorder
#     print "heightBox=",heigthBox
#     print "(mw, mh)=",(maxAllowedFontWidth,maxAllowdFontHeight)
    fs = 100
    while 1:
        if fs < 9:
            break
        fnt = ImageFont.truetype(g_params['font_dir'] + g_params['font'], fs)
        (fw, fh) = fnt.getsize("a")
        if fw <= maxAllowedFontWidth  and fh <= maxAllowdFontHeight:
            break
        else:
            fs -= 1
#     print "fs=",fs
    return fs
    #}}}

def AutoSizeFontTMBox(fontWidthAlign, fontHeightAlign, numSeq, specialProIdxDict, posTMList, TMnameList ): #{{{
    """Autosize the font for text written in TM box so that it fits the
       narrowest box """
# Get the maximum allowd fontWidth for each box
    maxAllowedFontWidthList = []
    margin = 1; #pixels
    # scale is roughly 50 seqs -> 0.5, 1500 seqs -> 1.5
    #scaleTMBox = myfunc.FloatDivision(numSeq, 1450)+ 27.0/58.0
    scaleTMBox = 1
    specialProIdxList = specialProIdxDict['reppro'] + specialProIdxDict['pdb'] + specialProIdxDict['final']

    fs = 50
    itr = 0
    MAX_ITR = 300
    while 1:
        if fs < 2:
            break
        fnt = ImageFont.truetype(g_params['font_dir'] + g_params['font'], fs)
        maxMargin = -9999999
        minMargin = 9999999
        for idx in specialProIdxList: # check all topologies with TMbox
            posTM = posTMList[idx]
            TMname = TMnameList[idx]
            for j in range(len(posTM)):
                (b, e) = posTM[j]
                try:
                    ss = TMname[j]
                except IndexError:
                    ss = "TM %d"%(j+1)
                ss = "%d"%(j+1)
                boxWidth = fontWidthAlign * (e-b)
                textWidth, textHeight = fnt.getsize(ss)
                margin = boxWidth - textWidth
                #print ("margin, boxwidth, textwidth)=", (margin, boxWidth, textWidth))
                if margin > maxMargin:
                    maxMargin = margin
                if margin < minMargin:
                    minMargin = margin
        #logger.debug("itr=%s"%str(itr) + " fs=%s"%str(fs)+ " (minMargin,maxMargin, fontHeightAlign) = "%(minMargin, maxMargin, fontHeightAlign))

        if minMargin < 5:
            fs -= 1
        elif minMargin >= fontHeightAlign*2:
            fs += 1
        else:
            break
        itr += 1
        if textHeight < fontHeightAlign*numSeq*0.1 and textHeight > fontHeightAlign*numSeq*0.05:
            break
        if itr > MAX_ITR:
            break

    g_params['font_size_TMbox'] = fs
    g_params['font_size_scalebar'] = int(fs*0.9+0.5)
    g_params['fntScaleBar'] = ImageFont.truetype(g_params['font_dir'] +
            g_params['font'], g_params['font_size_scalebar'])
    g_params['fntTMbox'] = ImageFont.truetype(g_params['font_dir'] +
            g_params['font'], g_params['font_size_TMbox'])
    g_params['fntTMbox_label'] = ImageFont.truetype(g_params['font_dir'] +
            "DejaVuSerif.ttf", g_params['font_size_TMbox']+1)
    fnt = ImageFont.truetype(g_params['font_dir'] + g_params['font'], fs)

    logger.debug("font_size_TMbox=%d", g_params['font_size_TMbox'])
    #print "fs=",fs
    return fnt.getsize("M")
#}}}
def AutoSizeFontDGProfileLabel(dgprofileRegionHeight):# {{{
    """
    Resize the font for the label of dg profile
    """
    margin = int(dgprofileRegionHeight*0.1+0.5)
    fs = 50
    itr = 0
    MAX_ITR = 300

    while 1:
        if fs < 2:
            break
        fnt = ImageFont.truetype(g_params['font_dir']+"DejaVuSerif-Bold.ttf", fs)
        textWidth = fnt.getsize(ylabel_DGprofile)[0]
        diff = dgprofileRegionHeight - textWidth

        if  diff < margin:
            fs -= 1
        elif diff >= margin*2:
            fs += 1
        else:
            break
        itr += 1
        #print ("itr=", itr, "fs=", fs, "diff=", diff, "margin=", margin)
        if itr > MAX_ITR:
            break
    g_params['fntDGprofileLable'] = fnt
    g_params['fntDGprofileTic'] = ImageFont.truetype(g_params['font_dir']+"DejaVuSerif.ttf", max(2, int(fs*0.65)))
    g_params['fntDGprofileLegend'] = ImageFont.truetype(g_params['font_dir']+"DejaVuSerif.ttf", max(2, int(fs*0.7)))

# }}}

def GetPositionIdenticalAdjacentNumber(lst, start, minLength): #{{{
# given a list of numbers 
# e.g.
# [ 1, 1, 0, 4,3, 4,2,3,3,3,54, 4,3, 44,44,44,44,3,3,3,3]
# get the position of identical adjacent numbers with at least minLength
    posList = []
    N = len(lst)
    if N  <= 0 :
        return posList
    i = 0
#     print 'N=', N
    while i < N:
        j = 0
        while i+j < N and lst[i+j] == lst[i]:
            j += 1
        if j > 0:
            if j >= minLength:
                posList.append((i+start, i+j+start))
            i += j
        else:
            i += 1
            
    return posList

#}}}

def GetRemainingSegmentList(start, end, posListToRemove):#{{{
    """get the remaining segments from start to end by removing segment defined
    in posListToRemove"""
    numPosToRemove = len(posListToRemove)
    if numPosToRemove < 1:
        return [(start, end)]
    length = end-start
    if length <= 0:
        return  [(start, end)]
    lst = [1]*length
#     print "====================="
#     print lst
#     print "====================="
    for (b, e) in posListToRemove:
        b1 = max(0, b-start)
        e1 = max(0, e-start)
        for i in range(b1, e1):
            lst[i] = 0
#     print lst
#     print "====================="
    posRemainList = []
    i = 0
    while i < length:
        j = 0
        while i+j < length and lst[i+j] == 1:
            j += 1
        if j > 0:
            posRemainList.append((i+start, i+j+start))
            i += j
        else:
            i += 1
    return posRemainList
#}}}

def CalDistPointToFragment(x, fragment):#{{{
    """
    MMMMMMM
   K            dist = 1
    K           dist = 0
          R     dist = 0
    """
    if x <= fragment[0]:
        dist = fragment[0]-x
    elif x >= fragment[1]:
        dist = x-fragment[1]+1
    else:
        d1 = x - fragment[0]
        d2 = (fragment[1]-1) - x
        dist = min(d1,d2)
    return dist
#}}}
def IsOutofMaxDistKR(posTM, x, maxDistKR):#{{{
    numTM = len(posTM)
    for i in range(numTM):
        d = CalDistPointToFragment(x, posTM[i])
        if d > 0 and d <= maxDistKR:
            if not lcmp.IsWithinTMRegion(x, posTM):
                return False
    return True
#}}}
def IsSafetoDeleteTheRegion(topoSeqList, start, end):#{{{
    """Check whether the deletion of one block in the topology MSA will affect
    the topology of any of the sequence"""
    GAP = g_params['GAP']
    lengthAlignment = len(topoSeqList[0])
    numSeq = len(topoSeqList)
    stateList = 'ioM'
    numState = len(stateList)
    if start < 1 or end >= lengthAlignment -1:
        return False
    for i in range(numSeq):
        topo = topoSeqList[i]
        cntFoundState = 0
        pList = [-1]*numState
        for j in range(numState):
            pList[j] = topo[start:end].find(stateList[j])
            if pList[j] >= 0:
                cntFoundState += 1
        if cntFoundState >= 3:
            return False
        elif cntFoundState >= 1:
            gaplesssubstr = topo[start:end].replace(GAP, '')
            if len(gaplesssubstr) > 0:
# get the first char and last char of the gapless substr
                firstT = gaplesssubstr[0]
                lastT = gaplesssubstr[len(gaplesssubstr)-1]
# check the first non GAP state on the right side and left side, if both are
# 'M', it is not safe to delete, otherwise safe 
# scan the left side
                j = start -1
                while j >= 1:
                    if topo[j] != GAP:
                        break
                    j -= 1
                firstLeftSideState = topo[j]
                p1 = j
# go the right side
                j = end
                while j < lengthAlignment -1:
                    if topo[j] != GAP:
                        break
                    j += 1
                firstRightSideState = topo[j]
                p2 = j
        # 1. leaving the beginning and end topology state unchanged
        # 2. do not remove the region of both sides are TM helices, otherwise,
        # two TM helices will be merged into one
                if (p1 == 0 or p2 == lengthAlignment-1 ):
                    return False
                else:
                    if cntFoundState == 2:
                        if not (lastT == firstRightSideState and
                                firstT == firstLeftSideState):
                            return False
                    elif cntFoundState == 1:
                        if not (lastT == firstRightSideState or
                                firstT == firstLeftSideState):
                            return False
    return True
#}}}
def IsSafetoDeleteTheRegionNew(origTopoSeqList, startOrig, endOrig,#{{{
        newTopoSeqList, startNew, per_K, per_R):
    """Check whether the deletion of one block in the topology MSA will affect
    the topology of any of the sequence"""
# The left side should be checked with the newTopoSeqList
# subsequence is obtained from origTopoSeqList[i][startOrig:endOrig]
# startNew is the position in the newTopoSeqList
    try:
        GAP = g_params['GAP']
        lengthAlignment = len(origTopoSeqList[0])
        numSeq = len(origTopoSeqList)
        stateList = 'ioM'
        numState = len(stateList)
        if g_params['isDrawKRBias'] and (sum(per_K[startOrig:endOrig]) +
            sum(per_R[startOrig:endOrig])) > 0.0:
            return False

        if startOrig < 1 or endOrig >= lengthAlignment -1:
            return False
        for i in range(numSeq):
            topoOrig = origTopoSeqList[i]
            topoNew = newTopoSeqList[i]
            cntFoundState = 0
            pList = [-1]*numState
            for j in range(numState):
                pList[j] = topoOrig[startOrig:endOrig].find(stateList[j])
                if pList[j] >= 0:
                    cntFoundState += 1
            if cntFoundState >= 3:
                return False
            elif cntFoundState >= 1:
                gaplesssubstr = topoOrig[startOrig:endOrig].replace(GAP, '')
                if len(gaplesssubstr) > 0:
# get the first char and last char of the gapless substr
                    firstT = gaplesssubstr[0]
                    lastT = gaplesssubstr[len(gaplesssubstr)-1]
# check the first non GAP state on the right side and left side, if both are
# 'M', it is not safe to delete, otherwise safe 
# scan the left side
                    j = startNew -1
                    while j >= 1:
                        if topoNew[j] != GAP:
                            break
                        j -= 1
                    if j >= 0:
                        firstLeftSideState = topoNew[j]
                    else:
                        firstLeftSideState = 'X'
                    p1 = j
# go the right side
                    j = endOrig
                    while j < lengthAlignment -1:
                        if topoOrig[j] != GAP:
                            break
                        j += 1
                    firstRightSideState = topoOrig[j]
                    p2 = j
            # 1. leaving the beginning and end topology state unchanged
            # 2. do not remove the region of both sides are TM helices, otherwise,
            # two TM helices will be merged into one
                    if (p1 < 0 or p2 == lengthAlignment-1 ):
                        return False
                    else:
                        if cntFoundState == 2:
                            if not (lastT == firstRightSideState and
                                    firstT == firstLeftSideState):
                                return False
                        elif cntFoundState == 1:
                            if not (lastT == firstRightSideState or
                                    firstT == firstLeftSideState):
                                return False
    except IndexError:
        return False
    return True
#}}}
def IsAtTMregionOfSpecialPro(i, topoSeqList, specialProIdxList):# {{{
    """
    Check if the residue position is at TM region of the speical proteins
    """
    for idx in specialProIdxList:
        if topoSeqList[idx][i] == "M":
            return True
    else:
        return False
# }}}
def GetPosTM_MSA(posTMList, specialProIdxList):# {{{
    """
    Get the beginning and end position of the MSA which has TM helices
    """
    beg = 9999999999
    end = -1
    for i in range(len(posTMList)):
        if not i in specialProIdxList:
            posTM = posTMList[i]
            if posTM[0][0] < beg:
                beg = posTM[0][0]
            if posTM[-1][1] > end:
                end = posTM[-1][1]
    return (beg, end)

# }}}
def ShrinkSeq(seq, shrinkedwidth):#{{{
    """Shrink the seq to shrinkedwidth"""
    N = len(seq)
    if N <= shrinkedwidth:
        return seq
    newseq = ""
    for i in range(shrinkedwidth):
        idx = int(round(i/float(shrinkedwidth-1)*(N-1)))
        newseq += seq[idx]
    return newseq
#}}}
def ShrinkGapInMSA_obsolete(topoSeqList, specialProIdxList=[]):#{{{
    """Shrink the gap regions
    topoSeqList will be updated and return the maparray"""
# For columns without 'M', shrink the region of each sequencs in the block to 
#     1. '' if there are no 'i' or 'o' in the block
#     2. 'i' or ' ' if there is no 'o' in the block
#     3. 'o' or ' ' if there is no 'i' in the block
#     4. 'io' or 'oi' or '  ' if there exist both 'i' and 'o' in the block
#     To be 'io' or 'i' or ' ' is depending on the subsequence in the region
#     for each topology
# further if there exist i or o but if the removal of this colum will not
    # change the topology of any sequence, this one can be removed.
# For columns with 'M', shrink the region with continous 'M' depending on the
# number of M in the  column
#     For the continous 'M' region, make a profile of 'M' percentage
#     For flat regions with length > 5, shrink them to  min(L, L/5*N/2)
#     For smooth profile with a peak, take the region above 50% 
#
# For TM regions in specialProIdxList, do not shrink it if it is at the
# beginning and end
    (cnt_i, cnt_o, cnt_M, cnt_GAP, 
            per_i, per_o, per_M, per_GAP) = lcmp.GetTopoStateFraction(
                    topoSeqList)
    lengthAlignment = len(topoSeqList[0])
    i = 0
    numSeq = len(topoSeqList)
    newList = [""]*numSeq
    posindexmap = {}
    num_specialpro = len(specialProIdxList)
    posTMList = [myfunc.GetTMPosition(x) for x in topoSeqList]
    (begTM_MSA, endTM_MSA) = GetPosTM_MSA(posTMList, specialProIdxList)

    cnt = 0
    while i < lengthAlignment:
        j = 0
        sumPer_i = 0.0
        sumPer_o = 0.0
        while i+j < lengthAlignment and per_M[i+j] == 0.0:
            sumPer_i += per_i[i+j]
            sumPer_o += per_o[i+j]
            j += 1
        if j >= 1:  #{{{
            if sumPer_i > 0.0 or sumPer_o > 0.0:
                if sumPer_i == 0.0:
                    for iseq in range(numSeq):
                        if topoSeqList[iseq][i:i+j].find("o") >= 0:
                            newList[iseq] += 'o'
                        else:
                            newList[iseq] += ' '
                    posindexmap[cnt] = i+j-1
                    cnt += 1
                elif sumPer_o == 0.0:
                    for iseq in range(numSeq):
                        if topoSeqList[iseq][i:i+j].find("i") >= 0:
                            newList[iseq] += 'i'
                        else:
                            newList[iseq] += ' '
                    posindexmap[cnt] = i+j-1
                    cnt += 1
                else:
                    for iseq in range(numSeq):
                        ss = topoSeqList[iseq][i:i+j]
                        p1 = ss.find('i')
                        p2 = ss.find('o')

                        if p1 >= 0 and p2 >= 0:
                            if p1 < p2:
                                newList[iseq] += 'io'
                            else:
                                newList[iseq] += 'oi'
                        else:
                            if p1 >= 0:
                                newList[iseq]+='ii'
                            elif p2 >= 0:
                                newList[iseq] += 'oo'
                            else:
                                newList[iseq] += '  '
                    posindexmap[cnt] = i
                    posindexmap[cnt+1] = i+j-1
                    cnt += 2
            i += j;#}}}
        else: # starts a region with M#{{{
            if num_specialpro > 0:
                if (IsAtTMregionOfSpecialPro(i, topoSeqList, specialProIdxList) and
                    (i < begTM_MSA and i>=endTM_MSA)):
                    for iseq in range(numSeq):
                        newList[iseq] =  topoSeqList[iseq][i]
                i += 1
            else:
                sumPer_i = 0.0
                sumPer_o + 0.0
                while i+j < lengthAlignment and per_M[i+j] > 0.0:
                    sumPer_i += per_i[i+j]
                    sumPer_o += per_i[i+j]
                    j += 1
#find all flat regions with >=5 residues
#             print cnt_M[i:i+j]
                posFlatRegionList = GetPositionIdenticalAdjacentNumber(
                        cnt_M[i:i+j], i, 5)
# get the rest regions
                posNonFlatRegionList = GetRemainingSegmentList(i, i+j,
                        posFlatRegionList)

                mergedRegionList = []
                for (b,e)in posFlatRegionList:
                    mergedRegionList.append(('flat', b, e))
                for (b,e)in posNonFlatRegionList:
                    mergedRegionList.append(('nonflat', b, e))
                mergedRegionList = sorted(mergedRegionList, key=lambda tup:tup[1])

                for (state, b, e) in mergedRegionList:
                    if state == 'flat':
                        shrinkedwidth = max(2, int(round((e-b)* min(1.0,
                            per_M[b]*10))))
                        for iseq in range(numSeq):
                            newList[iseq] += ShrinkSeq(topoSeqList[iseq][b:e],
                                    shrinkedwidth)
                        for k in range(cnt, cnt+shrinkedwidth):
                            posindexmap[k] = (b +
                                int(round((k-cnt)*(e-b)/float(shrinkedwidth-1))))
                        cnt += shrinkedwidth
                    else:
                        selectedPosList = []
                        minPerM = min(per_M[b:e])
                        maxPerM = max(per_M[b:e])
                        middlePerM = minPerM + (maxPerM - minPerM)*0.6 
                        for k in range(b,e):
                            if per_M[k] >= middlePerM:
                                selectedPosList.append(k)
                        selectedPosListSet = set(selectedPosList)
                        for k in range(b, e):
                            if (k in selectedPosListSet or not
                                    IsSafetoDeleteTheRegion(topoSeqList, k, k+1)):
                                for iseq in range(numSeq):
                                        newList[iseq] += topoSeqList[iseq][k]
                                posindexmap[cnt] = k
                                cnt += 1
                i += j
#}}}
    for iseq in range(numSeq):
        topoSeqList[iseq] = newList[iseq]
    return posindexmap

#}}}
def ShrinkGapInMSA_0(idList, topoSeqList, specialProIdxList=[]): #{{{
    """Shrink the gap regions
    topoSeqList will be updated and return the maparray
    by default ShrinkGapInMSA_0 is used"""
# For columns without 'M', shrink the region of each sequencs in the block to 
#     1. '' if there are no 'i' or 'o' in the block
#     2. 'i' or ' ' if there is no 'o' in the block
#     3. 'o' or ' ' if there is no 'i' in the block
#     4. 'io' or 'oi' or '  ' if there exist both 'i' and 'o' in the block
#     To be 'io' or 'i' or ' ' is depending on the subsequence in the region
#     for each topology
# further if there exist i or o but if the removal of this colum will not
    # change the topology of any sequence, this one can be removed.
# For columns with 'M', shrink the region with continous 'M' depending on the
# number of M in the  column
#     For the continous 'M' region, make a profile of 'M' percentage
#     For flat regions with length > 5, shrink them to  min(L, L/5*N/2)
#     For smooth profile with a peak, take the region above 50% 
#
    (cnt_i, cnt_o, cnt_M, cnt_GAP, 
            per_i, per_o, per_M, per_GAP) = lcmp.GetTopoStateFraction(
                    topoSeqList)
    isDrawKRBias = g_params['isDrawKRBias']

    num_specialpro = len(specialProIdxList)
    #print ("num_specialpro=%d"%(num_specialpro))
    posTMList = [myfunc.GetTMPosition(x) for x in topoSeqList]
    (begTM_MSA, endTM_MSA) = GetPosTM_MSA(posTMList, specialProIdxList)

    if isDrawKRBias:
        aaSeqDict = g_params['aaSeqDict']
        alignedSeqList = []
        for i  in range(len(idList)):
            toposeq = topoSeqList[i]
            seqid = idList[i]
            try:
                aaseq = aaSeqDict[seqid]
                aaseq = MatchToAlignedSeq(aaseq, toposeq, seqid)
                alignedSeqList.append(aaseq)
            except KeyError: 
                pass
        (cnt_K, cnt_R, per_K, per_R) = GetKRStateFraction(alignedSeqList)
    else:
        (cnt_K, cnt_R, per_K, per_R) = ([],[],[],[])


    lengthAlignment = len(topoSeqList[0])
    i = 0
    numSeq = len(topoSeqList)
    newList = [""]*numSeq
    posindexmap = {}

    cnt = 0
    while i < lengthAlignment:
        j = 0
        sumPer_i = 0.0
        sumPer_o = 0.0
        while i+j < lengthAlignment and per_M[i+j] == 0.0:
            sumPer_i += per_i[i+j]
            sumPer_o += per_o[i+j]
            j += 1
        if j >= 1:  #{{{ # non TM region
#             print "per_i:", per_i[i:i+j]
#             print "per_o:", per_o[i:i+j]
#             print "sumPer_i:", sumPer_i, "sumPer_o:", sumPer_o
#             print "Non M region: (%d, %d)"%(i,i+j)
            if sumPer_i > 0.0 or sumPer_o > 0.0: 
#                 print "With i or o: region: (%d, %d)"%(i,i+j)
                if not IsSafetoDeleteTheRegionNew(topoSeqList, i, i+j, newList,
                        cnt, per_K, per_R):
                    # otherwise, just delete this region
                    if isDrawKRBias:
                        repStatList = [] # state to be replaced
                        for iseq in range(numSeq):
                            subseq = topoSeqList[iseq][i:i+j].replace('-','')
                            if len(subseq) == 0:
                                repStatList.append(' ')
                            else:
                                repStatList.append(subseq[0])

                        tmpcnt = 0
                        for pp in range(i, i+j):
                            if per_K[pp] > 0.0 or per_R[pp] > 0.0:
                                for iseq in range(numSeq):
                                    newList[iseq] += repStatList[iseq]
                                posindexmap[cnt] = pp
                                cnt += 1
                                tmpcnt += 1
                        if tmpcnt == 0:
                            pp = i
                            for iseq in range(numSeq):
                                newList[iseq] += repStatList[iseq]
                            posindexmap[cnt] = pp
                            cnt += 1
                    else:
                        if sumPer_i == 0.0 or sumPer_o == 0.0:#{{{
                            for iseq in range(numSeq):
                                segment = topoSeqList[iseq][i:i+j]
                                if segment.find('o') >= 0:
                                    newList[iseq] += 'o'
                                elif segment.find('i') >= 0:
                                    newList[iseq] += 'i'
                                else:
                                    newList[iseq] += ' '
                            posindexmap[cnt] = i+j-1
                            cnt += 1
                        else:
                            stateList = 'io'
                            maxCntFoundState = 0
                            for iseq in range(numSeq):
                                cntFoundState = 0
                                segment = topoSeqList[iseq][i:i+j]
                                for state in stateList:
                                    if segment.find(state) >= 0:
                                        cntFoundState += 1
                                if cntFoundState > maxCntFoundState:
                                    maxCntFoundState = cntFoundState
                                if maxCntFoundState >= 2:
                                    break
                            if maxCntFoundState == 2:
                                for iseq in range(numSeq):
                                    ss = topoSeqList[iseq][i:i+j]
                                    p1 = ss.find('i')
                                    p2 = ss.find('o')

                                    if p1 >= 0 and p2 >= 0:
                                        if p1 < p2:
                                            newList[iseq] += 'io'
                                        else:
                                            newList[iseq] += 'oi'
                                    else:
                                        if p1 >= 0:
                                            newList[iseq]+='ii'
                                        elif p2 >= 0:
                                            newList[iseq] += 'oo'
                                        else:
                                            newList[iseq] += '  '
                                posindexmap[cnt] = i
                                posindexmap[cnt+1] = i+j-1
                                cnt += 2
                            else:
                                for iseq in range(numSeq):
                                    segment = topoSeqList[iseq][i:i+j]
                                    if segment.find('o') >= 0:
                                        newList[iseq] += 'o'
                                    elif segment.find('i') >= 0:
                                        newList[iseq] += 'i'
                                    else:
                                        newList[iseq] += ' '
                                posindexmap[cnt] = i+j-1
                                cnt += 1#}}}

            i += j;#}}}
        else: # starts a region with M#{{{
            if num_specialpro > 0:
                if (IsAtTMregionOfSpecialPro(i, topoSeqList, specialProIdxList) 
                    #and (i < begTM_MSA and i>=endTM_MSA)
                    ):
                    for iseq in range(numSeq):
                        newList[iseq] +=  topoSeqList[iseq][i]
                    posindexmap[cnt] = i
                    cnt += 1
                    print((i, lengthAlignment))
                i += 1
            else:
                sumPer_i = 0.0
                sumPer_o + 0.0
                while i+j < lengthAlignment and per_M[i+j] > 0.0:
                    sumPer_i += per_i[i+j]
                    sumPer_o += per_i[i+j]
                    j += 1
                if j > 0:
#                print "M region: (%d, %d)"%(i,i+j)
#find all flat regions with >=5 residues
                    posFlatRegionList = GetPositionIdenticalAdjacentNumber(
                            cnt_M[i:i+j], i, 5)
# get remaining regions
                    posNonFlatRegionList = GetRemainingSegmentList(i, i+j,
                            posFlatRegionList)
                    mergedRegionList = []
                    for (b,e)in posFlatRegionList:
                        mergedRegionList.append(('flat', b, e))
                    for (b,e)in posNonFlatRegionList:
                        mergedRegionList.append(('nonflat', b, e))
                    mergedRegionList = sorted(mergedRegionList, key=lambda
                            tup:tup[1])

#             if i >= 1320 and i+j <= 1460:
#                 print "region (%d, %d)"%(i, i+j)
#                 print cnt_M[i:i+j]
#                 print "posFlatRegionList:", posFlatRegionList
                    for (state, b, e) in mergedRegionList:
                        if state == 'flat':
                            if (per_GAP[b] > 0.95 and
                                    IsSafetoDeleteTheRegionNew(topoSeqList, b, e,
                                        newList, cnt, per_K, per_R)):
                                shrinkedwidth = 0
                            else:
                                shrinkedwidth = max(5, int(round((e-b)* min(1.0,
                                    per_M[b]*1.5))))
#                     if b >= 1320 and e <= 1460:
#                         print ("per_M[b]:",per_M[b], "len(%d, %d)="%(b, e),
#                         e-b, "shrinkedwidth=",shrinkedwidth)
                                selectedIndexList = [b +
                                        int(round(k*(e-b-1)/float(shrinkedwidth-1)))
                                        for k in range(shrinkedwidth)]
                                if isDrawKRBias:
                                    for pp in range(b, e):
                                        if (per_K[pp] + per_R[pp] > 0.0):
                                            selectedIndexList.append(pp)
                                    selectedIndexList = sorted(
                                            list(set(selectedIndexList)))
                                for k in range(b, e):
                                    if (k in selectedIndexList or
                                            not IsSafetoDeleteTheRegionNew(
                                                topoSeqList, k, k+1, newList, 
                                                cnt, per_K, per_R)):
                                        for iseq in range(numSeq):
                                            newList[iseq] += topoSeqList[iseq][k]
                                            posindexmap[cnt] = k
                                        cnt += 1
                        else: #'nonflat'
                            minPerM = min(per_M[b:e])
                            maxPerM = max(per_M[b:e])
                            middlePerM = minPerM + (maxPerM - minPerM)*0.5 
                            selectedIndexList = []
                            for k in range(b,e):
                                if ((per_GAP[k] < 0.95 and per_M[k] > middlePerM) or
                                    per_M[k] > 0.65 or
                                    (isDrawKRBias and (per_K[k]+per_R[k])>0.0)):
                                    selectedIndexList.append(k)
                            for k in range(b, e):
                                if (k in selectedIndexList or
                                        not IsSafetoDeleteTheRegionNew(topoSeqList,
                                            k, k+1, newList, cnt, per_K, per_R)):
                                    for iseq in range(numSeq):
                                        newList[iseq] += topoSeqList[iseq][k]
                                        posindexmap[cnt] = k
                                    cnt += 1
#                 if b >= 1320 and e <= 1460:
#                     print ("numSelectedColumn=", numSelectedColumn, maxPerM,
#                     "len(%d, %d)=%d"%(b, e, e-b))
                    i += j
                else:
                    i += 1
#}}}
    for iseq in range(numSeq):
        topoSeqList[iseq] = newList[iseq].replace(" ", "-")
        #print ("%10s: %s"%(idList[iseq], topoSeqList[iseq]))
    return posindexmap
#}}}
def ShrinkGapInMSA_exclude_TMregion(idList, topoSeqList): #{{{
    """
    Shrink non TM region and gap region
    topoSeqList will be updated and return the index map 
    Return posindexmap
    posindexmap           shrink -> non-shrink
    """
# updated 2016-09-13, shrink exclude signal peptide
# For columns without 'M', shrink the region of each sequencs in the block to 
#     1. '' if there are no 'i' or 'o' in the block
#     2. 'i' or ' ' if there is no 'o' in the block
#     3. 'o' or ' ' if there is no 'i' in the block
#     4. 'io' or 'oi' or '  ' if there exist both 'i' and 'o' in the block
#     To be 'io' or 'i' or ' ' is depending on the subsequence in the region
#     for each topology
# further if there exist i or o but if the removal of this colum will not
    # change the topology of any sequence, this one can be removed.
    (cnt_i, cnt_o, cnt_M, cnt_SP, cnt_GAP, 
            per_i, per_o, per_M, per_SP, per_GAP) = lcmp.GetTopoStateFraction_withSP(
                    topoSeqList)
    NMargin = 4 # keep <= NMargin/2 residues at two sides of aligned TM helices
    halfMargin = NMargin/2

    lengthAlignment = len(topoSeqList[0])
    i = 0
    numSeq = len(topoSeqList)
    newList = [""]*numSeq
    posindexmap = {}

    cnt = 0
    while i < lengthAlignment:
        j = 0
        sumPer_i = 0.0
        sumPer_o = 0.0
        while i+j < lengthAlignment and (per_M[i+j] == 0.0 and per_SP[i+j] == 0.0):
            sumPer_i += per_i[i+j]
            sumPer_o += per_o[i+j]
            j += 1
        poslist_to_keep = []
        if j >= 1:  # non TM region, non SP region
#             print "per_i:", per_i[i:i+j]
#             print "per_o:", per_o[i:i+j]
#             print "sumPer_i:", sumPer_i, "sumPer_o:", sumPer_o
#             print "Non M region: (%d, %d)"%(i,i+j)
            #if ((sumPer_i > 0.0 and sumPer_o == 0.0) or (sumPer_o > 0.0 and sumPer_i == 0.0)):
            if (sumPer_i > 0.0 or sumPer_o > 0.0):
                if i == 0:
                    poslist_to_keep = list(range(max(i+j-halfMargin,0), i+j))
                else:
                    poslist_to_keep = list(range(i,min(i+halfMargin,
                        lengthAlignment)))+list(range(max(i+j-halfMargin,0),i+j))
            else:
                poslist_to_keep = list(range(i, i+j))
            i += j
        else:
            poslist_to_keep = list(range(i, i+1))
            i += 1

        poslist_to_keep = sorted(set(poslist_to_keep), reverse=False)
        for pp in poslist_to_keep:
            for iseq in range(numSeq):
                try:
                    newList[iseq] += topoSeqList[iseq][pp]
                    posindexmap[cnt] = pp
                except IndexError:
                    print("Error! iseq=%d, pp=%d, lengthAlignment=%d"%(iseq, pp, lengthAlignment))
            cnt += 1

    for iseq in range(numSeq):
        topoSeqList[iseq] = newList[iseq].replace(" ", "-")
    return posindexmap
#}}}
def ShrinkMSA_Method_2(topoSeqList, aaSeqList=[], posTMList=[],#{{{
        shrinkrate_TM=2.0, max_hold_loop=3, isDrawKRBias=False):
    """
    Shrink multiple alignment of topologies
    Input:
        topoSeqList     A list of aligned topologies
        aaSeqList       A list of aligned sequences
        shrinkrate_TM   shrink rate for aligned TM regions, 
        max_hold_loop   maximal positions to hold for the loop region
    Output:
        updated topoSeqList
        idxmap_align2shrink       map original-aligned-position -> shrink-alignment
        idxmap_shrink2align       map shrink-alignemnt          -> original-alignment
    """
    if posTMList == []:
        for topo in topoSeqList:
            posTMList.append(myfunc.GetTMPosition(topo))

    numSeq = len(topoSeqList)

    lengthAlignment = len(topoSeqList[0])
# first get positions to keep
    array = ["l"]*lengthAlignment # l -- loop
                                  # P -- positively charged residues K or R
                                  # M -- common TM region
    for posTM in posTMList:
        for tup in posTM:
            for j in range(tup[0], tup[1]):
                array[j] = "M"
    if isDrawKRBias:
        for i in range(numSeq):
            seq = aaSeqList[i]
            for j in range(len(seq)):
                if (seq[j] in ["K", "R"] and
                        (not lcmp.IsWithinTMRegion(j, posTMList[i]))):
                    array[j] = "P"

    #debug
    if g_params['isPrintDebugInfo']:
        print("array for KR and TM regions, (l: loop, P: K or R, M, TM region)")
        print("%s"%("".join(array)))

    poslist_to_keep  = []
    i = 0
    while i < lengthAlignment:
        if array[i] == "M":
            j = 0
            while i+j<lengthAlignment and array[i+j] == "M":
                j+=1
            length_segment = j
            shrinked_len_seg = max(2, int(round(length_segment/shrinkrate_TM)))
            for k in range(shrinked_len_seg):
                poslist_to_keep.append(i +
                    int(round((length_segment-1)*k/float(shrinked_len_seg-1))))
            i += j
        else:
            j = 0
            while i+j<lengthAlignment and array[i+j] != "M":
                j+=1
            length_segment = j
            if length_segment < max_hold_loop:
                poslist_to_keep += list(range(i,i+j))
            else:
                for k in range(i,i+j):
                    if (k-i < max_hold_loop/2 or 
                            i+j-k < max_hold_loop or
                            (isDrawKRBias and array[k] == "P")):
                        poslist_to_keep.append(k)
            i += j

    idxmap_align2shrink = {}
    idxmap_shrink2align = {}
# begin debug:
#     print "ss=",ss
#     for i in xrange(len(poslist_to_keep)):
#         print "ss[ %d ] = %s"%(i, ss[i])
# end debug
    for i  in range(len(poslist_to_keep)):
        pp = poslist_to_keep[i]
        idxmap_align2shrink[pp] = i
        idxmap_shrink2align[i] = pp

    poslist_to_keep = sorted(set(poslist_to_keep), reverse=False)
    cnt = 0
    newList = [""]*numSeq
    for pp in poslist_to_keep:
        for iseq in range(numSeq):
            try:
                newList[iseq] += topoSeqList[iseq][pp]
            except IndexError:
                print("Error! iseq=%d, pp=%d, lengthAlignment=%d"%(iseq, pp,
                        lengthAlignment))
        cnt += 1

    for iseq in range(numSeq):
        topoSeqList[iseq] = newList[iseq].replace(" ", "-")
    return (idxmap_align2shrink, idxmap_shrink2align)
#}}}

def RunDGScan(aaseq, seqID):# {{{
    """
    Calculate the DG profile by using the dgscanProg
    return dgp
    """
    dgp = None
    tmpaaseqfile = tempfile.mktemp()
    tmpdgpfile = tempfile.mktemp()
    tmpfp = open(tmpaaseqfile, 'w')
    tmpfp.write(">%s\n"%seqID)
    tmpfp.write("%s\n"%aaseq)
    tmpfp.close()
    os.system("%s %s -lmin 21 -lmax 21 -o %s"%(g_params['dgscanProg'],
        tmpaaseqfile, tmpdgpfile))
    dgpDict = ReadInDGProfile(tmpdgpfile)
    os.system("rm -f %s %s" %(tmpaaseqfile, tmpdgpfile))
    if dgpDict and seqID in dgpDict:
        dgp = dgpDict[seqID]
    return dgp
# }}}

def DrawTMOfConsensus(posTM, xy0, fontWidth, fontHeight, draw): #{{{
    """Draw TM box"""
    widthAnnotation = g_params['widthAnnotation']
    annoSeqInterval = g_params['annoSeqInterval']
    font_size_TMbox = g_params['font_size_TMbox']
    fntTMbox = g_params['fntTMbox']
    heightTMbox = g_params['heightTMbox']
    (fontWidthTMbox, fontHeightTMbox) = fntTMbox.getsize("M")

    fntTMbox = g_params['fntTMbox']
    (x0,y0) = xy0
    x0 = x0 + widthAnnotation * fontWidth + annoSeqInterval * fontWidthTMbox
    y0 = y0 - fontHeightTMbox/2

    marginTop = 0
    marginBottom = 0

    cnt = 0
    for (b, e) in posTM:
        x1 = x0 + b*fontWidth
        y1 = y0 - marginTop
        x2 = x0 + e*fontWidth
        y2 = y0 + int(heightTMbox*fontHeightTMbox+0.5) + marginBottom
        box=[x1, y1 , x2, y2]
        draw.rectangle(box, fill="violet", outline="black")
# draw text
        s = "TM %d"%(cnt+1)
        (textwidth, textheight) = fntTMbox.getsize(s)
        textheight+=2
        x3 = int(round((x1+x2-textwidth)/2.0))
        y3 = int(round((y1+y2-textheight)/2.0))
        draw.text((x3, y3), s, font=fntTMbox, fill="black")
        cnt += 1
#}}}


def DrawDGProfile(dgpList, lengthAlignment, maxDG, minDG, xy0, #{{{ 
        dgprofileRegionWidth, dgprofileRegionHeight, spaceToLeftBorder,
        isDrawSeqID, line_mode, draw):
    """Draw DG profile
    support drawing multiple DG profiles in one box
    """
    logger.debug("Draw DG profile")
    (x0, y0) = xy0
    paddingtop = int(dgprofileRegionHeight*0.05+0.5)
    paddingbottom = int(dgprofileRegionHeight*0.05+0.5)

    heightDrawRegion = dgprofileRegionHeight - paddingtop - paddingbottom
    widthDrawRegion = dgprofileRegionWidth

    font_size = g_params['font_size_scalebar']
    line_width = max(1, int(dgprofileRegionHeight*0.02+0.5))
    outline_width = max(1, int(line_width*0.75+0.5))
    ticline_width = max(1,  int(outline_width*0.75+0.5))
    logger.debug("(minDG, maxDG)=(%f,%f)"%(minDG, maxDG))

    num_dgp = len(dgpList)

    blue = Color("blue")
    red = Color("red")
    colorpalette = list(blue.range_to(red,num_dgp))


# draw outline box
    x1 = x0
    y1 = y0 + paddingtop
    x2 = x1 + dgprofileRegionWidth
    y2 = y0 + dgprofileRegionHeight - paddingbottom
    box = [x1,y1,x2,y2]
    draw.rectangle(box, outline='black', width=outline_width)
    yMiddle = int(round((y1 + y2) / 2.0))

# draw x, axis
    x1 = x0
    y1 = y0 + paddingtop + int(round(heightDrawRegion*maxDG/(maxDG-minDG)))
    x2 = x1 + widthDrawRegion
    y2 = y1
    draw.line([x1, y1, x2, y2],fill="grey")

    yZero = y1

# draw ytics and text
    fnt = g_params['fntDGprofileTic']
    step = max(0.5, round((maxDG-minDG)/5))
    lengthtic = max(5, int(widthDrawRegion*0.01+0.5))
    ytic = 0.0
    while ytic <= maxDG:
        x1 = x0 - lengthtic
        y1 = yZero - int(round(heightDrawRegion*ytic/(maxDG-minDG))) 
        x2 = x1 + lengthtic
        y2 = y1
        draw.line([x1, y1, x2, y2],fill="black", width=ticline_width)
        text = "%.1f"%ytic
        (textWidth,textHeight) = fnt.getsize(text)
        draw.text((x1-textWidth-lengthtic,y1-textHeight/2), text, font=fnt,
                fill='black')
        ytic += step

    ytic = -step
    while ytic > minDG:
        x1 = x0 - lengthtic
        y1 = yZero - int(round(heightDrawRegion*ytic/(maxDG-minDG))) 
        x2 = x1 + lengthtic
        y2 = y1
        draw.line([x1, y1, x2, y2],fill="black")
        text = "%.1f"%ytic
        (textWidth,textHeight) = fnt.getsize(text)
        draw.text((x1-textWidth-lengthtic,y1-textHeight/2), text, font=fnt,
                fill='black')
        ytic -= step

    for ii in range(num_dgp):
        seqID, dgp = dgpList[ii]
        line_color = colorpalette[ii].get_hex_l() 

    # draw profile
        sizeSquare = int(g_params['image_scale']* 4+0.5)
        pointList= []
        papd = pointList.append
        for (idx, dg) in dgp:
            #print (idx, dg)
            h = int(round(dg/(maxDG-minDG)*heightDrawRegion))
            x1 = (x0 + int(round(widthDrawRegion*float(idx)/lengthAlignment)) -
                    sizeSquare/2)
            y1 = yZero - h - sizeSquare/2
            x2 = x1+sizeSquare
            y2 = y1+sizeSquare
            box=[x1,y1,x2,y2]
            if line_mode == "dot":
                draw.ellipse(box, outline=line_color)
            papd((x1+sizeSquare/2, y1+sizeSquare/2))
        if line_mode == "line":
            for i in range(0,len(pointList)-1,1):
                draw.line([pointList[i],pointList[i+1]],fill=line_color, width=line_width)

    # draw legend
    if g_params['isDrawDGProfileLegend']:
        fnt = g_params['fntDGprofileLegend']
        textWidth, textHeight = fnt.getsize("Initial Topology 1")
        bar_width = textWidth/3
        gap_item = textWidth/4
        gap_text_bar = bar_width/2
        item_width = textWidth + bar_width + gap_item + gap_text_bar
        x_padding = int(item_width*num_dgp*0.05+0.5)
        if isDrawSeqID and num_dgp > 1:
            # draw outline box
            x1 = x0
            y1 = y0 + dgprofileRegionHeight
            x2 = x1 + item_width*num_dgp+x_padding*2
            y2 = y1 + textHeight*2
            box = [x1,y1,x2,y2]
            draw.rectangle(box, outline='black', width=outline_width)

        for ii in range(num_dgp):
            seqID, dgp = dgpList[ii]
            line_color = colorpalette[ii].get_hex_l() 
            if (num_dgp) == 1:
                text = "Initial Topology"
            else:
                text = "Initial Topology %s"%(seqID.lstrip("rep"))
            (textWidth,textHeight) = fnt.getsize(text)
            x = x0 + ii*item_width + x_padding
            y = y0 + dgprofileRegionHeight + textHeight/4
            draw.text((x,y), text, font=fnt, fill='black')
            # draw bar
            x1 = x + textWidth + gap_text_bar
            y1 = y + textHeight/2
            x2 = x1 + bar_width
            y2 = y1
            draw.line([x1, y1, x2, y2],fill=line_color, width=line_width*2)


#}}}
def DrawTMOfConsensus2(topoSeq, posTM, typeTM, TMname, foldType, xy0, fontWidth, fontHeight, draw,length): #{{{
    """Draw TM box, for loop regions, GAPs are shown in blank
    """
    widthAnnotation = g_params['widthAnnotation']
    annoSeqInterval = g_params['annoSeqInterval']
    font_size_TMbox = g_params['font_size_TMbox']
    fntTMbox = g_params['fntTMbox']
    heightTMbox = g_params['heightTMbox']
    (fontWidthTMbox, fontHeightTMbox) = fntTMbox.getsize("M")

    fntTMbox = g_params['fntTMbox']
    (x0,y0) = xy0
    x0 = x0 + widthAnnotation * fontWidth + annoSeqInterval * fontWidthTMbox
    y0 = y0 - fontHeightTMbox/2

    marginTop = 0
    marginBottom = 0
    incolor="#F2EABD"
    #outcolor="#CCFFFF"
    incolor = g_params['loopcolor_in']
    outcolor = g_params['loopcolor_out']
    base_outline_color = "black"
    base_text_color = "black"
    base_outline_width = int(fontHeightTMbox*g_params['heightTMbox']*0.06+0.5)
    loop_height = fontHeightTMbox*0.4

    cnt = 0
    last=x0
    for (b, e) in posTM:
        x1 = x0 + b*fontWidth
        y1 = y0 - marginTop
        x2 = x0 + e*fontWidth
        y2 = y0 + int(heightTMbox*fontHeightTMbox +0.5) + marginBottom
        box=[x1, y1 , x2, y2]

        (text_TMlabel, text_color, outline_color, outline_width) = lcmp.SetMakeTMplotColor(
                cnt, TMname, foldType, base_outline_width, base_text_color, base_outline_color)

        logger.debug("text_TMlabel=%s, outline_color=%s, outline_width=%d"%(text_TMlabel, outline_color, outline_width))

        if (typeTM[cnt]=="M"): # out to in
            draw.rectangle(box, fill=g_params['memcolor_out_to_in'], outline=outline_color, width=outline_width)
            #draw.line([last,y2,x1,y2],g_params['loopcolor_in'])
            draw.rectangle([last,y2-loop_height,x1,y2],fill=incolor) #draw loop
        elif (typeTM[cnt]=="W"): # in to out
            draw.rectangle(box, fill=g_params['memcolor_in_to_out'], outline=outline_color, width=outline_width)
            #draw.line([last,y1,x1,y1],outcolor)
            draw.rectangle([last,y1,x1,y1+loop_height],fill=outcolor)
        elif (typeTM[cnt]=="R"): # Reeentrant inside
            y1 = y0 - marginTop + int(heightTMbox*fontHeightTMbox/3.0+0.5)
            y2 = y0 + int(heightTMbox*fontHeightTMbox+0.5) + marginBottom
            box=[x1, y1 , x2, y2]
            draw.rectangle(box, fill=incolor, outline=outline_color, width=outline_width)
            #draw.line([last,y2,x1,y2],incolor)
            draw.rectangle([last,y2-loop_height,x1,y2],fill=incolor) # draw loop
        elif (typeTM[cnt]=="r"): # Reentrant outside
            y1 = y0 - marginTop
            y2 = y0 + int(heightTMbox*fontHeightTMbox/3.0*2+0.5) + marginBottom
            box=[x1, y1 , x2, y2]
            draw.rectangle(box, fill=outcolor, outline=outline_color, width=outline_width)
            #draw.line([last,y1,x1,y1],outcolor)
            draw.rectangle([last,y1,x1,y1+loop_height],fill=outcolor) # draw loop
        else:
            draw.rectangle(box, fill="violet", outline=outline_color, width=outline_width)
        last=x2
        # draw text_TMlabel
        (textwidth, textheight) = fntTMbox.getsize(text_TMlabel)
        textheight+=2
        x3 = int(round((x1+x2-textwidth)/2.0))
        y3 = int(round((y1+y2-textheight)/2.0))
        draw.text((x3, y3), text_TMlabel, font=fntTMbox, fill=text_color)
        cnt += 1
    if (typeTM[cnt-1]=="R" or typeTM[cnt-1]=="W"):
        #draw.line([x2,y2,x0 + length*fontWidth,y2],incolor)
        draw.rectangle([x2,y2-loop_height,x0 + length*fontWidth,y2],fill=incolor)  # loop
    elif (typeTM[cnt-1]=="r" or (typeTM[cnt-1]=="M")):
        draw.line([x2,y1,x0 + length*fontWidth,y1],outcolor)
        draw.rectangle([x2,y1,x0 + length*fontWidth,y1+loop_height],fill=outcolor) # loop

    if g_params['isShowGap']:
        # show gap with blank, draw rectangle with fill=(0, 0, 0, 0)
        posGAP = myfunc.GetGapPosition(topoSeq)
        for (b, e) in posGAP:
            if not lcmp.IsWithinTMRegion(b, posTM):
                x1 = x0 + b*fontWidth
                y1 = y0 - marginTop
                x2 = x0 + e*fontWidth
                y2 = y0 + int(heightTMbox*fontHeightTMbox +0.5) + marginBottom

                io_left = lcmp.Get_IOState_upstream(topoSeq, b)
                io_right = lcmp.Get_IOState_downstream(topoSeq, e-1)
                if io_left == 'i' or io_right == 'i':
                    loop_type = 'inside'
                    y1 = y2 - loop_height
                else:
                    loop_type = 'outside'
                    y2 = y1 + loop_height
                box=[x1, y1 , x2, y2]
                #fillcolor = (0,0,0,0)
                fillcolor = (255,255,255,255)
                draw.rectangle(box, fill=fillcolor) # erase the gap region

#}}}
def DrawScale(length, posindexmap, xy0, font_size_alignment, #{{{
        fontWidth, fontHeight, draw):
    """Draw horizontal scale bar. font_size_alignment is the font_size for the
       MSA"""
    widthAnnotation = g_params['widthAnnotation']
    annoSeqInterval = g_params['annoSeqInterval']
    fntTMbox = g_params['fntTMbox']
    (fontWidthTMbox, fontHeightTMbox) = fntTMbox.getsize("a")

    isShrinked = False

    if len(posindexmap) > 0:
        length = len(posindexmap)
        isShrinked = True

    font_size_scalebar = g_params['font_size_scalebar']
    fntScaleBar = g_params['fntScaleBar']
    (fontWidthScaleBar, fontHeightScaleBar) = fntScaleBar.getsize("a")
    (x0,y0) = xy0
    x = x0 + widthAnnotation*fontWidth +annoSeqInterval* fontWidthTMbox
    y = y0
    fg = "black"
    step = 20 * max(1,int(math.ceil(fontWidthScaleBar / float(fontWidth) *
        len("%d"%length) /10 + 0.5)))
#     print "step=", step
    i = step
    while i < length:
        if isShrinked:
            s = "%s"%(posindexmap[i])
        else:
            s = "%s"%(i)
        #str=string.rjust(str,step," ")
        numdigit = len(s)
        shiftx = numdigit
        draw.text((x+step*fontWidth - shiftx*fontWidthScaleBar, y), s,
                font=fntScaleBar, fill=fg)
        i += step
        x += step * fontWidth

    y += int(fontHeightScaleBar*1.5+0.5)
    x = x0 + widthAnnotation*fontWidth + annoSeqInterval*fontWidthTMbox
    i = step
    shiftx = step-1
    while i < length:
        s = "|"
        draw.text((x+shiftx*fontWidth,y), s, font=fntScaleBar, fill=fg)
        i += step
        x += step * fontWidth
#}}}
def DrawHistogram(histoList, xy0, histoRegionWidth, histoRegionHeight, #{{{
        color_fill, color_outline, spaceToLeftBorder, draw):
    """Draw histogram"""
    (x0, y0) = xy0
    marginTop = max(20, int(round(histoRegionHeight * 0.1)))
    marginBottom = int(round(histoRegionHeight * 0.1))
    heightBox = histoRegionHeight - marginTop - marginBottom
    widthBox = histoRegionWidth
#   Draw outline box
    box = [x0, y0 + marginTop, x0 + widthBox, y0 + marginTop + heightBox]
    draw.rectangle(box, outline="black")

    yticList = [0, 50, 100]
    ylabel = "% state M"
# draw tics and text
    font_size = AutoSizeFontHistogram(ylabel, yticList, widthBox, heightBox,
            spaceToLeftBorder)
    #font_size = 12
    fnt = ImageFont.truetype(g_params['font_dir'] + g_params['font'],
            font_size)
    (fw, fh) = fnt.getsize('-')
    lengthtic = fw
    maxTicWidth = 0
    for ytic in yticList:
        x1 = x0 - lengthtic
        heightInPixel = int(round(ytic/100.0 * heightBox ))
        y1 = y0 + marginTop + heightBox - heightInPixel
        x2 = x0
        y2 = y1
        box=[x1, y1 , x2, y2]
        draw.line(box, fill="black")
        text = str(ytic)
        (textWidth, textHeight) = fnt.getsize(text)
        if textWidth > maxTicWidth:
            maxTicWidth = textWidth
        draw.text((x1-textWidth-lengthtic-3, y1-textHeight/2), text, font=fnt,
                fill='black')
    text = ylabel
    (textWidth, textHeight) = fnt.getsize(text)
    x = x0
    y = (y0 + marginTop + y0 + marginTop + heightBox) / 2
    xt = x - textWidth - 2*lengthtic - maxTicWidth - 3 - 10 
    yt = y - textHeight/2
    draw.text((xt , yt), text, font=fnt, fill='black')
    
#   Draw histogram bars
    x = x0
    for (w, h) in histoList:
        widthInPixel = int(round(w * widthBox))
        heightInPixel = int(round(h * heightBox ))
        x1 = x
        y1 = y0 + marginTop + heightBox - heightInPixel
        x2 = x1 + widthInPixel
        y2 = y0 + marginTop + heightBox
        box = [x1, y1, x2, y2]
        draw.rectangle(box, fill=color_fill, outline=color_outline)
        x += widthInPixel
#}}}

def GetSizeAnnotationToDraw(annotationList):#{{{
    """Get the size of annotation field"""
    maxSize = 0
    for anno in annotationList:
        m1 = re.search("nTM\s*=\s*[0-9]*", anno)
        m2 = re.search("group of [0-9]*", anno)
        m3 = re.search("[^\s]+", anno)
        pos1 = 0 
        pos2 = 0
        pos3 = 0
        if m1:
            pos1 = m1.end(0)
        if m2:
            pos2 = m2.end(0)
        if m3:
            pos3 = m3.end(0)
        size = max(pos1,pos2, pos3)
        if size > maxSize:
            maxSize = size
    return maxSize
#}}}
def GetSpecialProIndex(seqIDIndexDict):# {{{
    """
    Get index for special proteins
    return a dict with keys
    'reppro' = []
    'pdb' = []
    'final' = []
    """
    li_reppro = []
    li_pdb = []
    li_final = []

    for seqID in sorted(seqIDIndexDict.keys()):
        if seqID.find("rep") == 0:
            li_reppro.append(seqIDIndexDict[seqID])
        elif seqID.find("pdb_") == 0:
            li_pdb.append(seqIDIndexDict[seqID])
        elif seqID.find("final") == 0:
            li_final.append(seqIDIndexDict[seqID])

    dt = {}
    dt['reppro'] = li_reppro
    dt['pdb'] = li_pdb
    dt['final'] = li_final
    return dt

# }}}
def CalculateImageScale(numSeq):# {{{
    """
    Calculate Image Scale based on the number of sequences
    """
    if g_params['image_scale'] == None:
        if numSeq < 50:
            g_params['image_scale'] = 2.0
        elif numSeq < 100:
            g_params['image_scale'] = 1.5
        elif numSeq < 500:
            g_params['image_scale'] = 1.2
        else:
            g_params['image_scale'] = 1.0
    else:
        g_params['image_scale'] = 1.0
# }}}
def VerifyTerminalStatus(topoSeqList, posTMList):# {{{
    """
    Add i/o status before/after the TM helices if it is missing
    """
    numSeq = len(topoSeqList)
    lengthAlignment = len(topoSeqList[0])
    IOSTATE_LIST = ["i","o"]
    for i in range(numSeq):
        topo = topoSeqList[i]
        l_topo = [x for x in topo]
        posTM_Nterm = posTMList[i][0]
        posTM_Cterm = posTMList[i][-1]
        if lcmp.Get_IOState_upstream(topo, posTM_Nterm[0]) == '':
            io_after = lcmp.Get_IOState_downstream(topo, posTM_Nterm[1])
            io_before = IOSTATE_LIST[(IOSTATE_LIST.index(io_after)+1)%2]
            if posTM_Nterm[0] == 0:
                l_topo[posTM_Nterm[0]] = io_before
            else:
                l_topo[posTM_Nterm[0]-1] = io_before
        if lcmp.Get_IOState_downstream(topo, posTM_Cterm[1]) == '':
            io_before = lcmp.Get_IOState_upstream(topo, posTM_Cterm[0])
            io_after = IOSTATE_LIST[(IOSTATE_LIST.index(io_before)+1)%2]
            if posTM_Cterm[1] == lengthAlignment:
                l_topo[posTM_Cterm[1]-1] = io_after
            else:
                l_topo[posTM_Cterm[1]] = io_after
        topoSeqList[i] = ''.join(l_topo)
# }}}

def GetSeqTag(anno):#{{{
    "get sequence tag from annotation"
    tag = ""
    if anno.find("ClusterNo=") != -1:
        tag = re.search("ClusterNo=[0-9]+", anno).group(0)
    elif anno.find(" IDT ") != -1:
        tag = "IDT"
    elif anno.find("Consensus ") != -1:
        tag = "Consensus"
    elif anno.find(" OK ") != -1:
        tag = "OK"
    elif anno.find(" SHIFT ") != -1:
        tag = "SHIFT"
    elif anno.find(" INV ") != -1:
        tag = "INV"
    elif anno.find(" INV_SHIFT ") != -1:
        tag = "INV_SHIFT"
    elif anno.find(" DIFF ") != -1:
        tag = "DIFF"
    elif anno.find(" TM2GAP ") != -1:
        tag = "TM2GAP"
    elif anno.find(" TM2SEQ ") != -1:
        tag = "TM2SEQ"
    elif anno.find(" TM2GAP_AND_TM2SEQ ") != -1:
        tag = "TM2GAP_AND_TM2SEQ"
    elif anno.find(" Mixed ") != -1:
        tag = "Mixed"
    elif anno.find("Eukaryota") != -1:
        tag = "Eukaryota"
    elif anno.find("Archaea") != -1:
        tag = "Archaea"
    elif anno.find("Bacteria") != -1:
        tag = "Bacteria"
    else:
        tag = ""
    return tag
#}}}
def DrawTopology(anno, tag, toposeq, aaseq, xy0, fnt, fontWidth, #{{{
        fontHeight, isDrawText, draw):
    """Draw the topology MSA region with the PIL library"""
# draw background
    annoSeqInterval = g_params['annoSeqInterval']
    widthAnnotation = g_params['widthAnnotation']
    fntTMbox = g_params['fntTMbox']
    (fontWidthTMbox, fontHeightTMbox) = fntTMbox.getsize("a")

    (x0,y0) = xy0
    x = x0
    y = y0

    if isDrawText: # leave the width for the annotation text
        x += widthAnnotation * fontWidth
    else: # We actually need this shift anyhow to not get things shifted
        x += widthAnnotation * fontWidth
        
    #Draw a vertical bar for proteins in different groups
    if g_params['isDrawTagColumn']:
        if tag.find("ClusterNo") != -1:#{{{
            numCluster = int(tag.split("=")[1])
            numTM = int(re.search("nTM=[0-9]+", anno).group(0).split("=")[1])
            cntColor = 0
            fillColor = "white"
            if cntColor > 10 or numTM == 1:
                fillColor = "black"
            else:
                if numCluster == 1:
                    fillColor = "#008000"
                elif numCluster == 2:
                    fillColor = "#239C23"
                elif numCluster == 3:
                    fillColor = "#35A835"
                elif numCluster == 4:
                    fillColor = "#53B953"
                elif numCluster == 5:
                    fillColor = "#6CC66C"
                elif numCluster == 6:
                    fillColor = "#84D084"
                elif numCluster == 7:
                    fillColor = "#A5DEA5"
                elif numCluster == 8:
                    fillColor = "#CEEECE"
                elif numCluster == 9:
                    fillColor = "#E2F5E2"
                elif numCluster == 10:
                    fillColor = "#F5FCF5"
                else:
                    numCluster = "black"
            box=[x+fontWidthTMbox*1,y,x+fontWidthTMbox*3,y+fontHeight]
            draw.rectangle(box, fill=fillColor)
        else:
            fill_color = "white"
            if tag == "IDT":
                fill_color = "red"
            elif tag == "INV":
                fill_color ="blue"
            elif tag == "TM2GAP":
                fill_color ="green"
            elif tag == "TM2SEQ":
                fill_color ="violet"
            elif tag == "TM2GAP_AND_TM2SEQ":
                fill_color ="cyan"
            elif tag == "Consensus":
                fill_color = "purple"
            elif tag == "OK":
                fill_color = "red"
            elif tag == "SHIFT":
                fill_color = "pink"
            elif tag == "INV_SHIFT":
                fill_color ="lightgreen"
            elif tag == "DIFF":
                fill_color ="black"
            elif tag == "Archaea":
                fill_color ="blue"
            elif tag == "Bacteria":
                fill_color ="purple"
            elif tag == "Eukaryota":
                fill_color ="green"
            #print ("TEST",x,y,tag,fill_color)
            #fill_color='green'
            box=[x+fontWidthTMbox*1,y,x+fontWidthTMbox*3,y+fontHeight]
            draw.rectangle(box, fill=fill_color)

        x += annoSeqInterval * fontWidthTMbox
#}}}


    bg="#FFFFFF"; #white
    memcolor = "#FF0000"

    lengthSeq = len(toposeq)
    (posTM, typeTM) = lcmp.GetTMType(toposeq)
    seq_typeTM = [] # seq_typeTM is a list to name the type at every TM residue
                    # position, it defines only the types at the TM position, not loops
    seq_typeTM = [' '] * lengthSeq
    for jj in range(len(posTM)):
        tp = typeTM[jj]
        b,e = posTM[jj]
        for kk in range(b,e):
            seq_typeTM[kk] = tp


# it is much faster to draw a block of text than drawing characters one by one
    i=0
    if g_params['isColorByKingdom']:
        if ("Eukaryota" in anno):
            memcolor="#0000FF"; #blue
        elif ("Archaea" in anno):
            memcolor="#00FF00"; #Green
        elif ("Bacteria" in anno):
            memcolor="#FF0000"; #red
        else:  
            memcolor="#808080"; #grey

    while i < lengthSeq:
        j=i
        while j < lengthSeq and toposeq[j] == toposeq[i]:
            j += 1
        lengthSegment = j-i

        if toposeq[i] == "M":
            if seq_typeTM[i] == "M":
                bg=g_params['memcolor_out_to_in_MSA']
            elif seq_typeTM[i] == "W":
                bg = g_params['memcolor_in_to_out_MSA']
            #bg = "red"
        elif toposeq[i] == "i":
            #bg = "#F2EABD"; # light yellow
            bg = g_params['loopcolor_in_MSA']
        elif toposeq[i] == "o":
            #bg="#CCFFFF"; # faded blue
            bg = g_params['loopcolor_out_MSA']
        elif toposeq[i] == "S":
            #bg="#228B22"; # forestgreen for signal peptide
            bg = g_params['spcolor']

        else:
            if g_params['isColorWholeTMbox'] and lcmp.IsWithinTMRegion(i, posTM):
                bg = memcolor #"#FF0000"
            else:
                bg = "#FFFFFF"  #white
        box=[x,y,x+fontWidth*lengthSegment,y+fontHeight]
        draw.rectangle(box, fill=bg)
        x+=fontWidth*lengthSegment
        i=j

# draw text, foreground
    if isDrawText:
        x = x0
        y = y0
        #ss = string.ljust(anno[0:widthAnnotation], widthAnnotation, " ")
        ss = label[0:widthAnnotation].ljust(widthAnnotation, " ")
#    print "ss=%s, anno=%s" %(ss, anno)
        fg="#000000";# black
        draw.text((x,y), ss, font=fnt, fill=fg)
        x += (widthAnnotation*fontWidth)
        x += (annoSeqInterval*fontWidthTMbox)

        fg = "#000000" #black
# it is much faster to draw a block of text than drawing characters one by one
        if g_params['isShowTMIndex']: # show TM1 TM2 TM3 ... in the middle of the TMbox
            tmpli = [' ']*lengthSeq
            for kk in range(len(posTM)):
                (bb, ee) = posTM[kk]
                tmpstr = "TM%d"%(kk+1)
                mid = (bb+ee)/2
                bb1 = min(mid - len(tmpstr)/2, lengthSeq - len(tmpstr)-1)
                for jj in range(len(tmpstr)):
                    tmpli[bb1+jj] = tmpstr[jj]
            seq = "".join(tmpli)
        else:
            if aaseq != "":
                seq = aaseq
            else:
                seq = toposeq
            seq = seq.replace('-',' ')
        draw.text((x,y), seq, font=fnt, fill=fg)
#         i=0
#         while i < lengthSeq:
#             j=i
#             while j < lengthSeq and seq[j] == seq[i]:
#                 j+=1
#             lengthSegment=j-i
#             if seq[i] == "-":
#                 fg="#FFFFFF"
#             else:
#                 fg="#000000"
#             draw.text((x,y), seq[i:j], font=fnt, fill=fg)
#             x+=(fontWidth*lengthSegment)
#             i=j
#}}}
def CalculateImageParameter(fontWidth, fontHeight, lengthAlignment, numSeq, numSeprationLine, sectionSepSpace, specialProIdxDict, posTMList, TMnameList,  widthAdjustRatio):# {{{
    """
    Calculate image parameters for the PIL method
    """
    (fontWidthScaleBar, fontHeightScaleBar) = g_params['fntScaleBar'].getsize("a")
    fontHeightDGProfileLegend = g_params['fntDGprofileLegend'].getsize("a")[1]
    (fontWidthTMbox, fontHeightTMbox) = AutoSizeFontTMBox(fontWidth, fontHeight, numSeq, specialProIdxDict, posTMList, TMnameList)

    width = ((g_params['widthAnnotation'] + lengthAlignment) * (fontWidth) +
            g_params['annoSeqInterval']*fontWidthTMbox + int(widthAdjustRatio*g_params['marginX'] * 2+0.5))

    dgprofileRegionWidth = lengthAlignment * fontWidth
    dgprofileRegionHeight = max(30, int(width*0.15+0.5), int(round(numSeq * fontHeight * 0.2)), g_params['heightTMbox']*fontHeightTMbox*2)

    histoRegionWidth = lengthAlignment * fontWidth
    histoRegionHeight = max(50, int(round(lengthAlignment*fontHeight*widthAdjustRatio* 0.1)), int(round(numSeq*fontHeight* 0.1)))


    AutoSizeFontDGProfileLabel(dgprofileRegionHeight)

    width = ((g_params['widthAnnotation'] + lengthAlignment) * (fontWidth) +
            g_params['annoSeqInterval']*fontWidthTMbox + g_params['marginX'] * 2)

    height = (
            g_params['marginY']* int(myfunc.maprange((0.5,1), (2,6), sig2(numSeq, scale=5))) +
            #g_params['marginY']* 6+
            len(specialProIdxDict['reppro'])*int(g_params['heightTMbox']*fontHeightTMbox*1.5+0.5) +
            g_params['isDrawScaleBar']*int(fontHeightScaleBar*2.5+0.5)+
            g_params['isDrawMSA']*numSeq*fontHeight +
            g_params['isDrawSeprationLine'] * numSeprationLine * g_params['scaleSeprationLine']* fontHeight +
            (len(specialProIdxDict['pdb'])+len(specialProIdxDict['final']) >0)*(int(g_params['heightTMbox']*fontHeightTMbox+0.5)+sectionSepSpace*fontHeightScaleBar)+
            len(specialProIdxDict['pdb'])*int(g_params['heightTMbox']*fontHeightTMbox*1.5+0.5)+
            len(specialProIdxDict['final'])*int(g_params['heightTMbox']*fontHeightTMbox*1.5+0.5)+
            g_params['isDrawDGprofile'] *(int(dgprofileRegionHeight*1.1+0.5)+sectionSepSpace*fontHeightScaleBar)+
            g_params['isDrawDGprofile'] *((len(specialProIdxDict['reppro'])>1)*fontHeightDGProfileLegend*4)+
            g_params['isDrawPerMDistribution'] * (histoRegionHeight)
            )
    return (width, height, fontWidthTMbox, fontHeightTMbox, dgprofileRegionWidth, dgprofileRegionHeight, histoRegionWidth, histoRegionHeight)
# }}}
def DrawMSATopo_PIL(inFile, g_params):#{{{
    """Draw multiple alignment of topologies using the PIL library"""
    lcmp.SetMakeTMplotColor_g_params(g_params)
    isDrawSeqLable = True
    (idList, annotationList, topoSeqList) = myfunc.ReadFasta(inFile)

    H2W_ratio = g_params['H2W_ratio']
    widthAdjustRatio = 1.0

    topoSeqList = lcmp.RemoveUnnecessaryGap(topoSeqList)
    numSeq = len(idList)
    if numSeq < 1:
        print("No sequence in the file %s. Ignore." %(inFile), file=sys.stderr)
        return 1

    CalculateImageScale(numSeq)
    logger.debug("image_scale=%g"%(g_params['image_scale']))
    g_params['marginX'] = int(g_params['marginX']*g_params['image_scale']+0.5)
    g_params['marginY'] = int(g_params['marginY']*g_params['image_scale']+0.5)

    seqIDIndexDict = {}
    for i in range(numSeq):
        seqIDIndexDict[idList[i]] = i

    (fontWidthScaleBar, fontHeightScaleBar) = g_params['fntScaleBar'].getsize("a")
    sectionSepSpace = int(g_params['image_scale']*g_params['heightTMbox']/2+0.5)

    rootname = os.path.basename(os.path.splitext(inFile)[0])
    aaSeqDict = GetAASeqDict(inFile)

    if g_params['outpath'] == "":
        outpath = myfunc.my_dirname(inFile)
    else:
        outpath = g_params['outpath']

    str_krbias = ""
    if g_params['isDrawKRBias'] == True:
        str_krbias = ".krbias"
    outFile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, g_params['outFormat'])



    # ============================
    # get index of special proteins
    #   * representative protein
    #   * topology of a PDB structure
    #   * the final topology
    specialProIdxDict = GetSpecialProIndex(seqIDIndexDict)
    specialProIdxList = specialProIdxDict['reppro'] + specialProIdxDict['pdb'] + specialProIdxDict['final']
    if len(specialProIdxList) == 0:
        if g_params['makeCleanPlot']:
            # make all as final proteins for Sudha's case
            specialProIdxDict['final'] = list(range(len(idList)))
#             print("FATAL ERROR! specialPro is not set, TM box can not be plotted. "\
#                     "Please check your input file %s"%(inFile), file=sys.stderr)
#             return 1
        else:
            specialProIdxDict['reppro'] = [0]

    specialProIdxList = specialProIdxDict['reppro'] + specialProIdxDict['pdb'] + specialProIdxDict['final']
    posTMList = [myfunc.GetTMPosition(x) for x in topoSeqList]
    TMnameList = [] # note that TMname is also a list, TMnameList is a list of list
    foldTypeList = []
    for i in range(numSeq):
        if i in specialProIdxList:
            TMname = myfunc.GetTMnameFromAnnotation(annotationList[i]) # annotation from the file topomsa
            foldType = myfunc.GetFoldTypeFromAnnotation(annotationList[i])
            TMnameList.append(TMname)
            foldTypeList.append(foldType)
        else:
            TMnameList.append([])
            foldTypeList.append("")

    origPosTMList = [myfunc.GetTMPosition(x) for x in topoSeqList]
    VerifyTerminalStatus(topoSeqList, origPosTMList)

# posindexmap: map of the residue position to the original MSA
# e.g. pos[0] = 5 means the first residue is actually the 6th residue position
# in the original MSA
#   backup original aligned topoSeqList
    origTopoSeqList = []
    for seq in topoSeqList:
        origTopoSeqList.append(seq)


    posindexmap = {}
    if g_params['isShrink']:
        if g_params['method_shrink'] == 0:
            posindexmap = ShrinkGapInMSA_0(idList, topoSeqList, specialProIdxList=[])
        elif g_params['method_shrink'] == 1:
            posindexmap = ShrinkGapInMSA_exclude_TMregion(idList, topoSeqList)

    # get posTMList for the shink version of MSA
    posTMList = [myfunc.GetTMPosition(x) for x in topoSeqList]

#     for i in range(len(topoSeqList)):
#         print ("%10s: %s" %(idList[i], topoSeqList[i]))

    g_params['widthAnnotation'] = GetSizeAnnotationToDraw(annotationList)
    if g_params['makeCleanPlot']:
        g_params['widthAnnotation'] = 4

    widthAnnotation = g_params['widthAnnotation']
    tagList = []

    for seqAnno in annotationList:
        tagList.append(GetSeqTag(seqAnno))
    numSeprationLine = len(set(tagList))

    lengthAlignment = len(topoSeqList[0])

    fnt = ImageFont.truetype(g_params['font_dir'] + g_params['font'],
            int(g_params['image_scale']*g_params['font_size']))
    (fontWidth, fontHeight) = fnt.getsize("a")
    (width, height, fontWidthTMbox, fontHeightTMbox, dgprofileRegionWidth, dgprofileRegionHeight, histoRegionWidth, histoRegionHeight) = CalculateImageParameter(fontWidth, fontHeight, lengthAlignment, numSeq,  numSeprationLine, sectionSepSpace, specialProIdxDict, posTMList, TMnameList, 1.0)

    if (H2W_ratio != None and height/float(width)!=H2W_ratio):
        widthAdjustRatio = height/float(width)/H2W_ratio
        fontWidth = int(fontWidth * widthAdjustRatio + 0.5)
        g_params['marginY'] += int(widthAdjustRatio*10+0.5)
        (width, height, fontWidthTMbox, fontHeightTMbox, dgprofileRegionWidth, dgprofileRegionHeight, histoRegionWidth, histoRegionHeight) = CalculateImageParameter(fontWidth, fontHeight, lengthAlignment, numSeq,  numSeprationLine, sectionSepSpace,specialProIdxDict, posTMList, TMnameList,  widthAdjustRatio)

    isDrawText = g_params['isDrawText']
    font_size = g_params['font_size']
    if g_params['isAutoSize']:
        while height*width > g_params['MAXIMAGESIZE']:
            if font_size > 3:
                font_size -= 1
                fnt = ImageFont.truetype(g_params['font_dir'] +
                        g_params['font'], font_size)
                (fontWidth, fontHeight) = fnt.getsize("a")
            else:
                if fontWidth > 1:
                    fontWidth -= 1
                if fontHeight > 1:
                    fontHeight -= 1

            (width, height, fontWidthTMbox, fontHeightTMbox, dgprofileRegionWidth, dgprofileRegionHeight, histoRegionWidth, histoRegionHeight) = CalculateImageParameter(fontWidth, fontHeight, lengthAlignment, numSeq,  numSeprationLine, sectionSepSpace,specialProIdxDict, posTMList, TMnameList,  1.0)

            if font_size < 3:
                isDrawText = False
            if fontWidth < 2 and fontHeight < 2:
                break
        logger.debug("height (%d) *width (%d) = %d"%(height, width, height*width))

        if (H2W_ratio != None and height/float(width)!=H2W_ratio):
            widthAdjustRatio = height/float(width)/H2W_ratio
            fontWidth = int(fontWidth * widthAdjustRatio + 0.5)
            g_params['marginY'] += int(widthAdjustRatio*10+0.5)
            (width, height, fontWidthTMbox, fontHeightTMbox, dgprofileRegionWidth, dgprofileRegionHeight, histoRegionWidth, histoRegionHeight) = CalculateImageParameter(fontWidth, fontHeight, lengthAlignment, numSeq,  numSeprationLine, sectionSepSpace, specialProIdxDict, posTMList, TMnameList,  widthAdjustRatio)

        if height*width > g_params['MAXIMAGESIZE']:
            msg = "%s: (fontWidth, fontHeight) have been reduced to (%d, %d)"\
                  ", but the image size is still too big (%dM)"%(
                          inFile, fontWidth, fontHeight, height*width/1024/1024)
            print(msg)
        else:
            msg = "font is autosized to %d, isDrawText = %s, "\
                    "(fontWidth, fontHeight) = (%d, %d)"%(
                            font_size,isDrawText, fontWidth, fontHeight)
            print(msg.format(font_size, isDrawText, fontWidth, fontHeight))


    bg_color="#FFFFFF"; # white
    if g_params['mode'] == "RGB":
        newImage = Image.new(g_params['mode'], (width, height),bg_color)
    elif g_params['mode'] == "P":
        newImage = Image.new(g_params['mode'], (width, height),255)

    draw = ImageDraw.Draw(newImage); # setup to draw on the main image
    x = g_params['marginX']
    y = g_params['marginY']

# Draw TM helices of the consensus.(or the representative topologies)
# Do not draw consensus if there is only final topologies in the topology
# alignment
    if not len(specialProIdxDict['final']) == len(specialProIdxList):
        if len(specialProIdxDict['reppro']) > 0:
            idxRepProList = specialProIdxDict['reppro']
        else:
            idxRepProList = [0] #if reppro is not set, take the first one
        cnt = 0
        for idx in idxRepProList:
            (posTM_rep,typeTM_rep) = lcmp.GetTMType(topoSeqList[idx])
            if isDrawSeqLable:
                xt = g_params['marginX'] + fontWidth*g_params['widthAnnotation']*0
                if not g_params['makeCleanPlot']:
                    if len(idxRepProList) == 1:
                        label = "Initial Topology"
                    else:
                        label = "Initial Topology %d"%(cnt+1)
                    label = rootname
                else:
                    label = ""
                #ss = string.ljust(label[0:widthAnnotation], widthAnnotation, " ")
                ss = label[0:widthAnnotation].ljust(widthAnnotation, " ")
                fg="#000000";# black
                draw.text((xt,y), ss, font=g_params['fntTMbox_label'], fill=fg)

            DrawTMOfConsensus2(topoSeqList[idx], posTM_rep, typeTM_rep,
                    TMnameList[idx], foldTypeList[idx], (x,y),
                    fontWidth, fontHeight, draw,lengthAlignment)
            y += int(g_params['heightTMbox'] * fontHeightTMbox*1.5+0.5)
            cnt += 1
        #y += sectionSepSpace*fontHeightTMbox

# Draw a scale bar of the residue position
    (fontWidthScaleBar, fontHeightScaleBar) = g_params['fntScaleBar'].getsize("a")
    if g_params['isDrawMSA']:
        if g_params['isDrawScaleBar']:
            DrawScale(lengthAlignment, posindexmap, (x,y), font_size, fontWidth,
                    fontHeight, draw)
            y += int(fontHeightScaleBar*2.5+0.5)

        maxDistKR = g_params['maxDistKR'] 
        isDrawKRBias = g_params['isDrawKRBias']

        tagFormer = tagList[0]
        for i in range(numSeq):
            tagCurrent = tagList[i]
            seqID = idList[i]

            if i in specialProIdxList: #do not draw the topology of special proteins 
                continue

            if tagCurrent != tagFormer:
                tagFormer = tagCurrent
                if g_params['isDrawSeprationLine'] == True:
                    box = [x, y+1, width-marginX, y+fontHeight*scaleSeprationLine-1]
                    draw.rectangle(box, fill="grey",outline="white")
                    y += fontHeight

            anno = annotationList[i]
            tag = tagList[i]
            toposeq = topoSeqList[i]
            aaseq = ""
            if seqID in aaSeqDict:
                aaseq = aaSeqDict[seqID]
                #print aaseq
                #print origTopoSeqList[i]
                # origTopoSeqList is the original (non-shinked MSA)
                aaseq = MatchToAlignedSeq(aaseq, origTopoSeqList[i], seqID)
                if posindexmap  != {}:
                    tmpaaseq = ""
                    for pp in range(len(posindexmap)):
                        aa = aaseq[posindexmap[pp]]
                        if (isDrawKRBias and  aa in ["K", "R"] and
                                IsOutofMaxDistKR(posTMList[i], posindexmap[pp],
                                    maxDistKR)):
                            aa = " "
                            if g_params['isPrintDebugInfo']:
                                print(seqID, aa, posindexmap[pp], posTMList[i])
                        tmpaaseq += (aa)
                    aaseq = tmpaaseq
                if isDrawKRBias:
                    aaseq = HideNonKRResidue(aaseq)
#        print aaseq
            DrawTopology(anno, tag, toposeq, aaseq, (x,y), fnt, fontWidth,
                    fontHeight, isDrawText, draw)
#       if tagCurrent == "Consensus":
#           y+=fontHeight
            y += fontHeight
            if y >= height:
                print(("Error! position y(%d) exceeds height (%d)" %
                        (y, height)), file=sys.stderr)

# Draw special topologies
    idxPDBList = specialProIdxDict['pdb']
    idxFinalProList = specialProIdxDict['final']
    if len(idxPDBList)+len(idxFinalProList) > 0:
        y += int(g_params['heightTMbox']*fontHeightTMbox+0.5)
        logger.debug("idxPDBList=%s, idxFinalProList=%s", str(idxPDBList), str(idxFinalProList))
        for idx in idxPDBList+idxFinalProList:
            if idx != -1:
                if isDrawSeqLable:
                    xt = g_params['marginX'] + fontWidth*g_params['widthAnnotation']*0
                    if idx in idxPDBList:
                        label = idList[idx].lstrip('pdb_')
                    else:
                        if g_params['makeCleanPlot']:
                            label = alphabet[idx]
                        else:
                            if len(idxFinalProList) == 1:
                                label = "Final Topology"
                            else:
                                label = "Final Topology "+ idList[idx].lstrip("final")
                    #ss = string.ljust(label[0:widthAnnotation], widthAnnotation, " ")
                    ss = label[0:widthAnnotation].ljust(widthAnnotation, " ")
                    fg="#000000";# black
                    draw.text((xt,y), ss, font=g_params['fntTMbox_label'], fill=fg)

                (posTM,typeTM) = lcmp.GetTMType(topoSeqList[idx])
                # draw topology of the representative protein
                DrawTMOfConsensus2(topoSeqList[idx], posTM, typeTM,
                        TMnameList[idx], foldTypeList[idx], (x,y), fontWidth,
                        fontHeight, draw, lengthAlignment)

                y += int(g_params['heightTMbox']*fontHeightTMbox*1.5+0.5)

        y += sectionSepSpace*fontHeightScaleBar
# Draw DGprofile
    if g_params['isDrawDGprofile']:
        dgpList = []
        for idx in specialProIdxDict['reppro']:
            seqID = idList[idx]
            toposeq = topoSeqList[idx]
            lengthAlignment = len(toposeq)
            DGProfileFile = GetDGProfileFileName(inFile, seqID)
            isDrawSeqID = True
            dgprofileDict = None
            logger.debug("seqID=%s, DGProfileFile=%s"%(seqID, DGProfileFile))
            if os.path.exists(DGProfileFile):
                dgprofileDict = ReadInDGProfile(DGProfileFile)
            dgp = None
            if dgprofileDict:
                if seqID in dgprofileDict:
                    dgp = dgprofileDict[seqID]
                elif 'query' in dgprofileDict:
                    dgp = dgprofileDict['query']
            elif seqID in aaSeqDict: #if dg profile file is not provided, calculate it
                aaseq = aaSeqDict[seqID]
                dgp = RunDGScan(aaseq, seqID)
            if dgp:
                idxmap_aligne2seq = lcmp.GetAlign2SeqMap(origTopoSeqList[idx],
                        origTopoSeqList[idx].replace(GAP,""))  
                aligned_dgp = MatchAlignedDGP(dgp, idxmap_aligne2seq, posindexmap, toposeq)
                dgpList.append([seqID, aligned_dgp])
        if len(dgpList) > 0:
            dgList = []
            for ii in range(len(dgpList)):
                dgp = dgpList[ii][1]
                ldg = [x[1] for x in dgp]
                dgList += ldg
            minDG = min(-1.0, min(dgList))
            maxDG = max(2, max(dgList))
            line_mode = "line"
            x = (g_params['marginX'] + g_params['widthAnnotation'] * fontWidth +  g_params['annoSeqInterval'] *
                    fontWidthTMbox)
            spaceToLeftBorder = (g_params['annoSeqInterval'] * fontWidthTMbox +
                    g_params['widthAnnotation'] * fontWidth)
            DrawDGProfile(dgpList, lengthAlignment, maxDG, minDG,
                    (x,y), dgprofileRegionWidth,
                    dgprofileRegionHeight, spaceToLeftBorder, isDrawSeqID,
                    line_mode, draw)

            # Add ylabel deltaG (kcal/mol)
            fnt_label = g_params['fntDGprofileLable']
            (fw_l, fh_l) = fnt_label.getsize(ylabel_DGprofile)
            image2 = Image.new('RGBA', (int(fw_l*1.1+0.5)+5, int(fh_l*1.2+0.5)+5))
            draw2 = ImageDraw.Draw(image2)
            draw2.text((5, 5), text=ylabel_DGprofile, font=fnt_label, fill="black")
            image2 = image2.rotate(90, resample=1, expand=1)
            logger.debug ("widthAdjustRatio=%f"%widthAdjustRatio)
            px, py = int(g_params['marginX']*1.5), y + 0*int(fw_l/2)
            sx, sy = image2.size
            newImage.paste(image2, (px, py), mask=image2)
            del image2

            y += dgprofileRegionHeight
        y += sectionSepSpace*fontHeightScaleBar

# draw distribution of 'M' percentage
    if g_params['isDrawPerMDistribution']:
        (cnt_i, cnt_o, cnt_M, cnt_GAP, per_i, per_o, per_M, per_GAP) =\
                lcmp.GetTopoStateFraction( topoSeqList)
# histoList is a list of 2-tuples (width, height) where height is a value
# ranging to 1, and sum of width equals 1.
        histoList = []
        for i in range(lengthAlignment):
            histoList.append((1.0/lengthAlignment, per_M[i]))
        x = (g_params['marginX'] + g_params['widthAnnotation'] * fontWidth +  g_params['annoSeqInterval'] *
                fontWidthTMbox)
        color_fill = 'red'
        color_outline = 'black'
        spaceToLeftBorder = (g_params['annoSeqInterval'] * fontWidthTMbox + g_params['widthAnnotation']
                * fontWidth)
        DrawHistogram(histoList, (x,y), histoRegionWidth, histoRegionHeight,
                color_fill, color_outline, spaceToLeftBorder,
                draw)

    if not g_params['isQuiet']:
        print(("Topology MSA is drawn and output to \"%s\"" %(outFile)))

    del draw
    newImage.save(outFile)
    del newImage
#}}}
def DrawMSATopo_SVG(inFile, g_params):#{{{
    (idList, annotationList, topoSeqList) = myfunc.ReadFasta(inFile)
    topoSeqList = lcmp.RemoveUnnecessaryGap(topoSeqList)
    numSeq = len(idList)
    if numSeq < 1:
        print("No sequence in the file %s. Ignore." %(inFile), file=sys.stderr)
        return 1

    marginX = g_params['marginX']
    marginY = g_params['marginY']
    annoSeqInterval = g_params['annoSeqInterval']
    widthAnnotation = g_params['widthAnnotation']
    heightScaleBar = g_params['heightScaleBar']
    heightTMbox = g_params['heightTMbox']
    scaleSeprationLine = g_params['scaleSeprationLine']
    font_size_scalebar = g_params['font_size_scalebar']
    fntScaleBar = g_params['fntScaleBar']
    (fontWidthScaleBar, fontHeightScaleBar) = fntScaleBar.getsize("a")

    rootname = os.path.basename(os.path.splitext(inFile)[0])
    aaSeqDict = GetAASeqDict(inFile)

    #rootname=rootname.split('.')[0]
    if g_params['outpath'] == "":
        outpath = myfunc.my_dirname(inFile)
    else:
        outpath = g_params['outpath']

    str_krbias = ""
    if g_params['isDrawKRBias'] == True:
        str_krbias = ".krbias"
    svgfile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'svg')
    pdffile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'pdf')

# posindexmap: map of the residue position to the original MSA
# e.g. pos[0] = 5 means the first residue is actually the 6th residue position
# in the original MSA
#   backup original aligned topoSeqList
    alignedTopoSeqList = []
    posTMList = []
    for seq in topoSeqList:
        alignedTopoSeqList.append(seq)
        posTMList.append(myfunc.GetTMPosition(seq))

    posindexmap = {}
    if g_params['isShrink']:
        posindexmap = ShrinkGapInMSA_0(idList, topoSeqList)

    posTM = myfunc.GetTMPosition(topoSeqList[0])
    g_params['widthAnnotation'] = GetSizeAnnotationToDraw(annotationList)
    widthAnnotation = g_params['widthAnnotation']
    tagList = []
    for seqAnno in annotationList:
        tagList.append(GetSeqTag(seqAnno))
#     print tagList
    numSeprationLine = len(set(tagList))
    lengthAlignment = len(topoSeqList[0])

    svg_document = pysvg.structure.svg()
#     shape_builder = pysvg.builders.ShapeBuilder()
# 
#     svg_document.addElement(shape_builder.createRect(0, 0,
#                                                  "200px", "100px",
#                                                  strokewidth = 1,
#                                                  stroke = "black",
#                                                  fill = "rgb(255, 255, 0)"))
    myStyle = pysvg.builders.StyleBuilder()
    myStyle.setFontFamily(fontfamily="monospace")
#    myStyle.setFontSize('5pt')
    myStyle.setFillOpacity('0.5')
    myStyle.setFilling(fill="grey")
#    myStyle.setStroke("green")
#    myStyle.fill = "blue"
    
    
    newAnnoList = []
    for i in range(len(topoSeqList)):
        newAnnoList.append("%s %s"%(idList[i], tagList[i]))
    maxSizeAnno = max([len(x) for x in newAnnoList])

    print("maxSizeAnno=",maxSizeAnno)

    x0 = 10
    y0 = 10
    for i in range(len(topoSeqList)):
        y = y0 + (i)*20
        string = "%-*s %s"%(maxSizeAnno+5, newAnnoList[i], topoSeqList[i])
        #print len(string)
        print(string)
        t1 = pysvg.text.text(string, x=x0, y=y)
        t1.set_style(myStyle.getStyle())
        svg_document.addElement(t1)
    svg_document.save(svgfile)
    cmd = "svg2pdf %s %s"%(svgfile, pdffile)
    os.system(cmd)
    print("%s output"%(pdffile))

#}}}
def DrawMSATopo_MAT(inFile, g_params):#{{{
    (idList, annotationList, topoSeqList) = myfunc.ReadFasta(inFile)
    numSeq = len(idList)
    if numSeq < 1:
        print("No sequence in the file %s. Ignore." %(inFile), file=sys.stderr)
        return 1
    topoSeqList = lcmp.RemoveUnnecessaryGap(topoSeqList)

    lengthAlignmentOriginal = len(topoSeqList[0])

    maxDistKR = g_params['maxDistKR']
    marginX = g_params['marginX']
    marginY = g_params['marginY']
    annoSeqInterval = g_params['annoSeqInterval']
    widthAnnotation = g_params['widthAnnotation']
    heightScaleBar = g_params['heightScaleBar']
    heightTMbox = g_params['heightTMbox']
    scaleSeprationLine = g_params['scaleSeprationLine']
    font_size_scalebar = g_params['font_size_scalebar']
    fntScaleBar = g_params['fntScaleBar']
    (fontWidthScaleBar, fontHeightScaleBar) = fntScaleBar.getsize("a")



    pdfcrop_margin_left = g_params['pdfcrop_margin_left']
    pdfcrop_margin_top = g_params['pdfcrop_margin_top']
    pdfcrop_margin_right = g_params['pdfcrop_margin_right']
    pdfcrop_margin_bottom = g_params['pdfcrop_margin_bottom']

    rootname = os.path.basename(os.path.splitext(inFile)[0])
    aaSeqDict = GetAASeqDict(inFile)

    if g_params['outpath'] == "":
        outpath = myfunc.my_dirname(inFile)
    else:
        outpath = g_params['outpath']

    isDrawKRBias = g_params['isDrawKRBias']
    str_krbias = ""
    if g_params['isDrawKRBias'] == True:
        str_krbias = ".krbias"
    svgfile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'svg')
    pdffile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'pdf')
    txtfile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'txtplot')

# posindexmap: map of the residue position to the original MSA
# e.g. pos[0] = 5 means the first residue is actually the 6th residue position
# in the original MSA
#   backup original aligned topoSeqList
    alignedTopoSeqList = []
    posTMList = []
    for seq in topoSeqList:
        alignedTopoSeqList.append(seq)
        posTMList.append(myfunc.GetTMPosition(seq))

    posindexmap = {}
    method_shrink = g_params['method_shrink']
    if g_params['isDrawDGprofile']:
        method_shrink = 1

    if g_params['isShrink']:
        if g_params['method_shrink'] == 0:
            posindexmap = ShrinkGapInMSA_0(idList, topoSeqList)
        elif g_params['method_shrink'] == 1:
            posindexmap = ShrinkGapInMSA_exclude_TMregion(idList, topoSeqList)

    posTM = myfunc.GetTMPosition(topoSeqList[0])
    g_params['widthAnnotation'] = GetSizeAnnotationToDraw(annotationList)
    widthAnnotation = g_params['widthAnnotation']
    tagList = []
    for seqAnno in annotationList:
        tagList.append(GetSeqTag(seqAnno))
#     print tagList
    numSeprationLine = len(set(tagList))
    lengthAlignment = len(topoSeqList[0])

    newAnnoList = []
    for i in range(len(topoSeqList)):
        newAnnoList.append("%s %s"%(myfunc.GetFirstWord(annotationList[i]), tagList[i]))
    maxSizeAnno = max([len(x) for x in newAnnoList])

    print("maxSizeAnno=",maxSizeAnno)
    fonttype = 'monospace'

    numSeq = len(topoSeqList)

    aaSeqList = [] # list of amino acid sequences, aligned and shrinked if enabled
    final2seq_idxMapList = [] # the index map from the final (shrinked or not)
                              # sequence to the orignal unaligned sequence.
    for i in range(numSeq):
        seqID = idList[i]
        idxmap_aligne2seq = lcmp.GetAlign2SeqMap(alignedTopoSeqList[i],
                alignedTopoSeqList[i].replace(GAP,""))  
        idxmap = {}
        if seqID in aaSeqDict:
            aaseq = aaSeqDict[seqID]
            #print aaseq
            #print alignedTopoSeqList[i]
            # alignedTopoSeqList is the original (non-shinked MSA)
            aaseq = MatchToAlignedSeq(aaseq, alignedTopoSeqList[i], seqID)
            if posindexmap  != {}:
                tmpaaseq = ""
                for pp in range(len(posindexmap)):
                    aa = aaseq[posindexmap[pp]]
                    idxmap[pp] = idxmap_aligne2seq[posindexmap[pp]]
                    if (isDrawKRBias and  aa in ["K", "R"] and
                            IsOutofMaxDistKR(posTMList[i], posindexmap[pp],
                                maxDistKR)):
                        aa = " "
                        if g_params['isPrintDebugInfo']:
                            print(seqID, aa, posindexmap[pp], posTMList[i])
                    tmpaaseq += (aa)
                aaseq = tmpaaseq
            else:
                idxmap = idxmap_aligne2seq
            if isDrawKRBias:
                aaseq = HideNonKRResidue(aaseq)
            aaSeqList.append(aaseq)
        else:
            aaSeqList.append("")
        final2seq_idxMapList.append(idxmap)

# setting font properties
    #ffam = "monospace"
    #ffam = "Fixed"
    ffam = "Courier New"
    fontpath = "%s/%s"%(g_params['font_dir'], "Courier_New.ttf")
    fp = matplotlib.font_manager.FontProperties(
        fname=fontpath, family=ffam, style='normal', size=12,
        weight='normal', stretch='normal')

    ffam = "Arial"
    fontpath = "%s/%s"%(g_params['font_dir'], "Arial.ttf")
    fp_anno = matplotlib.font_manager.FontProperties(
        fname=fontpath, family=ffam, style='normal', size=14,
        weight='normal', stretch='normal')
    
# get the text width and height in pixels
    x=0
    y=0
    linespaceInPixel = 12
    ss = "M"*1
    pth = matplotlib.textpath.TextPath((x, y), ss, prop=fp)
    bb = pth.get_extents(transform=None)
    widthSingleCharInPixel =  float(bb.width)/len(ss)
    heightSingleCharInPixel = float(bb.height)
    print("charwidth, charheight", widthSingleCharInPixel, heightSingleCharInPixel)



    widthAnnoInPixel = 0
    for i in range(numSeq):
        ss = newAnnoList[i]
        pth = matplotlib.textpath.TextPath((x, y), ss, prop=fp_anno)
        bb = pth.get_extents(transform=None)
        wtd =  float(bb.width)/len(ss.replace(" ",""))*len(ss)
        if wtd > widthAnnoInPixel:
            widthAnnoInPixel = wtd
    print("widthAnnoInPixel=", widthAnnoInPixel)

# set aspect ratio
    if g_params['isDrawDGprofile']:
        heightRatios = [numSeq,10]
        gs = gridspec.GridSpec(2, 1, height_ratios=heightRatios) 
    else:
        heightRatios = [1]

    print("heightRatios=", heightRatios)

    sumTextHeightInPixel = (heightSingleCharInPixel + linespaceInPixel)*(numSeq+1)
    sumTextWidthInPixel = (widthSingleCharInPixel)*(lengthAlignment+1)
#    sumTextWidthAnnotationInPixel = widthSingleCharInPixel*(maxSizeAnno+5)
    sumTextWidthAnnotationInPixel = widthAnnoInPixel+20

    print("lengthAlignment=", lengthAlignment)
    print("sumTextWidthAnnotationInPixel=", sumTextWidthAnnotationInPixel)
    print("sumTextWidthInPixel=", sumTextWidthInPixel)
    print("sumTextHeightInPixel=", sumTextHeightInPixel)


    widthUnitFigureInPixel = 8*80
    heightUnitFigureInPixel = 6*80

    adjust_right = 0.99
    #adjust_left = float(sumTextWidthInPixel)/(lengthAlignment*widthSingleCharInPixel)
    adjust_left = (adjust_right) * (sumTextWidthAnnotationInPixel/(sumTextWidthAnnotationInPixel+lengthAlignment*widthSingleCharInPixel))
    adjust_top = max(1.0 - float(2)/numSeq, 0.7)
    adjust_bottom = min(float(2)/numSeq,0.3)
    print("adjust_left=",adjust_left)
    print("adjust_right=",adjust_right)
    print("adjust_top=",adjust_top)
    print("adjust_bottom=",adjust_bottom)

    subplot1_width_ratio = (adjust_right-adjust_left)
    subplot1_height_ratio = float(heightRatios[0])/sum(heightRatios)*(adjust_top-adjust_bottom)
#subplot1_width_ratio = 1.0/(1.0+0.2+0.2+adjust_left)
    widthUnitSubplot1InPixel = widthUnitFigureInPixel*subplot1_width_ratio
    heightUnitSubplot1InPixel = heightUnitFigureInPixel*subplot1_height_ratio

    widthscale = float(sumTextWidthInPixel)/widthUnitSubplot1InPixel+0.00
    heightscale = float(sumTextHeightInPixel)/heightUnitSubplot1InPixel+0.02
    print("sumTextWidthInPixel, sumTextHeightInPixel=", (sumTextWidthInPixel, sumTextHeightInPixel))
    print("scale=",(widthscale, heightscale))

    widthSubplot1InPixel = widthUnitSubplot1InPixel * widthscale
    heightSubplot1InPixel = heightUnitSubplot1InPixel * heightscale

    print("widthSubplot1InPixel=",widthSubplot1InPixel)
    print("heightSubplot1InPixel", heightSubplot1InPixel)

    widthSingleCharInAxes = float(widthSingleCharInPixel)/widthSubplot1InPixel
    heightSingleCharInAxes = float(heightSingleCharInPixel)/heightSubplot1InPixel
    widthAnnotationInAxes = float(sumTextWidthAnnotationInPixel)/widthSubplot1InPixel
    linespaceInAxes = float(linespaceInPixel)/heightSubplot1InPixel

    print("widthSingleCharInAxes, heightSingleCharInAxes=", (widthSingleCharInAxes, heightSingleCharInAxes))
    

# create figure object
    figsize = (8*widthscale, 6*heightscale) # fig size in inches (width,height)
    fig = plt.figure(figsize = figsize) # set the figsize
    fig.subplots_adjust(left=adjust_left, right=adjust_right, top=adjust_top, bottom=adjust_bottom)

    if g_params['isDrawDGprofile']:
        ax = fig.add_subplot(gs[0])
    else:
        ax = fig.add_subplot(111)

    #ax.axis('off') 

    inv = ax.transAxes.inverted()

    loc_xtics = []
    label_xtics = []
    if posindexmap != {}:
        for jj in range(0, lengthAlignment, 10):
            loc_xtics.append(jj)
            label_xtics.append(posindexmap[jj])

    ax.set_xlim(0, lengthAlignment)
    #ax.set_xlabel("Sequence position", fontsize=16)
    if posindexmap != {}:
        plt.xticks(np.array(loc_xtics), np.array(label_xtics))
    ax2 = ax.twiny()
    ax2.set_xlabel("Alignment position", fontsize=16)
    ax2.set_xlim(0,lengthAlignment)
    if posindexmap != {}:
        plt.xticks(np.array(loc_xtics), np.array(label_xtics))
    ax.set_ylim(numSeq,0)
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
# make a plot of sequence indexes
    l1 = []
    l2 = []
    for j in range(lengthAlignment):
        if posindexmap != {}:
            idxAlignedSeq = posindexmap[j]
        else:
            idxAlignedSeq = j
        l1.append(idxAlignedSeq)
        l2.append(0)
    plt.plot(l1,l2, ' ')

    x0 = 0
    y0 = 1.0 - linespaceInAxes - heightSingleCharInAxes

    x = x0
    y = y0
    yshift=0

    for i in range(len(topoSeqList)):
        y -= yshift
        anno = "%-*s"%(maxSizeAnno+5, newAnnoList[i])
        x = x0 - widthAnnotationInAxes
        plt.text(x, y, anno, fontproperties=fp_anno, transform=ax.transAxes)

        xshift=0
        x = x0
        if aaSeqList[i] != "":
            seq = aaSeqList[i]
        else:
            seq = topoSeqList[i]
        topo = topoSeqList[i]
        posTM = myfunc.GetTMPosition(topo)
        if len(posTM) == 0:
            txt = seq
            plt.text(x, y, txt, fontproperties=fp, transform=ax.transAxes)
            width = widthSingleCharInAxes * len(txt)
            height = heightSingleCharInAxes + linespaceInAxes
            edgecolor = 'none'
            y2 = y - linespaceInAxes/2.0
            rec = matplotlib.patches.Rectangle((x, y2), width, height,
                    transform=ax.transAxes)
            ax.add_patch(rec)
            xshift = width
            yshift = height
        else:
            li = []
            for (b, e) in posTM:
                li.append(b)
                li.append(e)
            for j in range(len(li)+1):
                if j == 0:
                    begin = 0
                else:
                    begin = li[j-1]

                if j != len(li):
                    end = li[j]
                else:
                    end = lengthAlignment

                txt = seq[begin:end]
                txt = txt.replace("-", " ")

                txt_topo = topo[begin:end].replace(GAP," ")

                x_char = x
                for char in txt:
                    plt.text(x_char, y, char, fontproperties=fp, transform=ax.transAxes)
                    x_char += widthSingleCharInAxes


                #t = matplotlib.textpath.TextPath((0,0), txt, prop=fp)
                #bb = t.get_extents(transform=inv)

                if txt_topo.find('M')!=-1:
                    edgecolor = 'black'
                    facecolor = 'red'
                    color = 'red'
                elif txt_topo.find('i') != -1:
                    color = 'lightyellow'
                elif txt_topo.find('o') != -1:
                    color = 'lightblue'
                else:
                    facecolor = 'none'
                    edgecolor = 'none'
                width = widthSingleCharInAxes * len(txt)
                #width = bb.width
                #print "w1, w2=", widthSingleCharInAxes*len(txt), bb.width
                height = heightSingleCharInAxes + linespaceInAxes
                y2 = y - linespaceInAxes/2.0
#                 rec = matplotlib.patches.Rectangle((x, y2), width, height,
#                         facecolor=facecolor, edgecolor=edgecolor, alpha=0.5,transform=ax.transAxes)
                rec = matplotlib.patches.Rectangle((x, y2), width, height,
                        color=color, transform=ax.transAxes)
                ax.add_patch(rec)
                xshift = width
                x += xshift
                yshift = height
                cnttextlength += len(txt)
        #print "%3d: %4d %4d %s"%(i, cnttextlength, len(seq), seq), posTM

    if g_params['isDrawDGprofile']:
        dgprofileDict = {} #{{{
        if os.path.exists(g_params['DGProfileFile']):
            dgprofileDict = ReadInDGProfile(g_params['DGProfileFile'])
        for i in range(numSeq):
            seqID = idList[i]
            toposeq = topoSeqList[i]
            lengthAlignment = len(toposeq)
            if (not seqID in dgprofileDict) and (seqID in aaSeqDict):
                aaseq = aaSeqDict[seqID]
                #print "aaseq=", aaseq
                tmpaaseqfile = tempfile.mktemp()
                tmpdgpfile = tempfile.mktemp()
                tmpfp = open(tmpaaseqfile, 'w')
                tmpfp.write(">%s\n"%seqID)
                tmpfp.write("%s\n"%aaseq)
                tmpfp.close()
                os.system("%s %s -lmin 21 -lmax 21 -o %s"%(g_params['dgscanProg'],
                    tmpaaseqfile, tmpdgpfile))
                tmpDGPDict = ReadInDGProfile(tmpdgpfile)
                os.system("rm -f %s %s" %(tmpaaseqfile, tmpdgpfile))
                for seqid in tmpDGPDict:
                    dgprofileDict[seqid] = tmpDGPDict[seqid]
                    #print dgprofileDict[seqid]
        #}}}
        ax = fig.add_subplot(gs[1])
        #ax.set_xlim(0, lengthAlignmentOriginal)
        ax.set_xlim(0, lengthAlignment)
        ax.set_xlabel("Alignment position", fontsize=16)
        ax.set_ylabel(r"${\Delta}G$", fontsize=16)
        if posindexmap != {}:
            plt.xticks(np.array(loc_xtics), np.array(label_xtics))
        for i in range(numSeq):
            seqid = idList[i]
            alignedTopoOriginal = alignedTopoSeqList[i]
            align2seqMap = lcmp.GetAlign2SeqMap(alignedTopoOriginal,
                    alignedTopoOriginal.replace(GAP,""))
            #print "dgprofileDict[%s]"%seqid, dgprofileDict[seqid]
            print("align2seqMap=", align2seqMap)
            try: 
                dgp = dgprofileDict[seqid]
                dt = {}
                for tup in dgp:
                    dt[tup[0]] = tup[1]

                x = []
                y = []
                for j in range(lengthAlignment):
                    if posindexmap != {}:
                        try:
                            idxAlignedSeq = posindexmap[j]
                        except KeyError:
                            #print "j=%d not in posindexmap"%j
                            pass
                    else:
                        idxAlignedSeq = j
                    try:
                        idxSeq = align2seqMap[idxAlignedSeq]
                    except KeyError:
#                         if g_params['isPrintDebugInfo']:
#                             print "idxAlignedSeq=%d not in align2seqMap"%idxAlignedSeq
                        pass

                    if idxSeq in dt:
                        #x.append(idxAlignedSeq)
                        x.append(j)
                        y.append(dt[idxSeq])
                    else:
#                         if g_params['isPrintDebugInfo']:
#                             print "idxSeq=%d not in dgp, idxAlignedSeq=%d"%(idxSeq, idxAlignedSeq)
                        pass
                if i < len(colorList_DG_profile):
                    color = colorList_DG_profile[i]
                else:
                    color = 'none'
# plot by line
                plt.plot(x,y, label=seqid, color=color, linewidth=2.0)
# plot with '+' symbol
                # plt.plot(x,y, '+', label=seqid, color=color)
                plt.hlines(0, 0, lengthAlignmentOriginal)
                plt.legend()
            except KeyError:
                print("no dgprofile for %s"%(seqid))
                pass

    plt.savefig(pdffile)
    print("%s output"%(pdffile))
    cmd = "pdfcrop --margins '%d %d %d %d' --clip %s"%(pdfcrop_margin_left,
            pdfcrop_margin_top, pdfcrop_margin_right, pdfcrop_margin_bottom,
            pdffile)
    os.system(cmd)


#   Write Txtformat alignment
    #print final2seq_idxMapList
    htmlfile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'html')
#     WriteTXTAlignment(idList, newAnnoList, topoSeqList, alignedTopoSeqList,
#             aaSeqList, final2seq_idxMapList, txtfile)
    if len(idList) == 2:
        #WriteHTMLAlignment2(idList, newAnnoList, topoSeqList,
        #        alignedTopoSeqList, aaSeqList, final2seq_idxMapList, htmlfile)
        tmpmapList = []
        for i in range(len(alignedTopoSeqList)):
            tmpmap = {}
            for j in range(len(alignedTopoSeqList[i])):
                tmpmap[j] = j
            tmpmapList.append(tmpmap)
        WriteHTMLAlignment2(idList, newAnnoList, alignedTopoSeqList,
                alignedTopoSeqList, aaSeqList, tmpmapList, htmlfile)
    elif len(idList) > 2:
        WriteHTMLAlignment3(idList, newAnnoList, topoSeqList,
                alignedTopoSeqList, aaSeqList, final2seq_idxMapList, htmlfile)


#}}}
def DrawMSATopo_MAT2(inFile, g_params):#{{{
    """
    draw topology similar as
    blue upperhyphen    outside
    red  _____          inside 
    grey box            TM helix (in->out)
    white box           TM helix (out->in)
    """
    logger = logging.getLogger(__name__)
    (idList, annotationList, topoSeqList) = myfunc.ReadFasta(inFile)
    numSeq = len(idList)
    if numSeq < 1:
        logger.debug("No sequence in the file %s. Ignore." %(inFile))
        return 1
    topoSeqList = lcmp.RemoveUnnecessaryGap(topoSeqList)

    lengthAlignmentOriginal = len(topoSeqList[0])

    maxDistKR = g_params['maxDistKR']
    marginX = g_params['marginX']
    marginY = g_params['marginY']
    annoSeqInterval = g_params['annoSeqInterval']
    widthAnnotation = g_params['widthAnnotation']
    heightScaleBar = g_params['heightScaleBar']
    heightTMbox = g_params['heightTMbox']
    scaleSeprationLine = g_params['scaleSeprationLine']
    font_size_scalebar = g_params['font_size_scalebar']
    fntScaleBar = g_params['fntScaleBar']
    (fontWidthScaleBar, fontHeightScaleBar) = fntScaleBar.getsize("a")

    pdfcrop_margin_left = g_params['pdfcrop_margin_left']
    pdfcrop_margin_top = g_params['pdfcrop_margin_top']
    pdfcrop_margin_right = g_params['pdfcrop_margin_right']
    pdfcrop_margin_bottom = g_params['pdfcrop_margin_bottom']

    rootname = os.path.basename(os.path.splitext(inFile)[0])
    aaSeqDict = GetAASeqDict(inFile)

    #rootname=rootname.split('.')[0]
    if g_params['outpath'] == "":
        outpath = myfunc.my_dirname(inFile)
    else:
        outpath = g_params['outpath']

    isDrawKRBias = g_params['isDrawKRBias']
    str_krbias = ""
    if g_params['isDrawKRBias'] == True:
        str_krbias = ".krbias"
    svgfile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'svg')
    pdffile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'pdf')
    txtfile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'txtplot')

# posindexmap: map of the residue position to the original MSA
# e.g. pos[0] = 5 means the first residue is actually the 6th residue position
# in the original MSA
#   backup original aligned topoSeqList
    alignedTopoSeqList = []
    posTMList = []
    for seq in topoSeqList:
        alignedTopoSeqList.append(seq)
        posTMList.append(myfunc.GetTMPosition(seq))

    posindexmap = {}
    method_shrink = g_params['method_shrink']
#     if g_params['isDrawDGprofile']:
#         method_shrink = 1

    aaSeqAlignList = [] # list of original aligned aa seq list
    for i in range(numSeq):
        seqid = idList[i]
        try:
            aaseq = aaSeqDict[seqid]
            aaseq = MatchToAlignedSeq(aaseq, alignedTopoSeqList[i], seqid)
        except KeyError:
            aaseq = " "*lengthAlignmentOriginal
        aaSeqAlignList.append(aaseq)


    if g_params['isShrink']:
        if g_params['method_shrink'] == 0:
            posindexmap = ShrinkGapInMSA_0(idList, topoSeqList)
        elif g_params['method_shrink'] == 1:
            posindexmap = ShrinkGapInMSA_exclude_TMregion(idList, topoSeqList)
        elif g_params['method_shrink'] == 2:
            (idxmap_align2shrink, idxmap_shrink2align) =\
                    ShrinkMSA_Method_2(topoSeqList, aaSeqAlignList, posTMList,
                            g_params['shrinkrate_TM'], g_params['max_hold_loop'],
                            g_params['isDrawKRBias'])
            posindexmap = idxmap_shrink2align

    posTM = myfunc.GetTMPosition(topoSeqList[0])
    g_params['widthAnnotation'] = GetSizeAnnotationToDraw(annotationList)
    widthAnnotation = g_params['widthAnnotation']
    tagList = []
    for seqAnno in annotationList:
        tagList.append(GetSeqTag(seqAnno))
#     print tagList
    numSeprationLine = len(set(tagList))
    lengthAlignment = len(topoSeqList[0])

# added 2013-12-04
    if g_params['shrinkrate'] != None:
        shrinkrate = g_params['shrinkrate'] # shrink the sequence proportionally
    else:
        shrinkrate = lengthAlignment/120.0
    lengthAlignmentShrinked = int(lengthAlignment/shrinkrate+0.5)
    logger.debug("lengthAlignment=%d"%lengthAlignment)
    logger.debug("shrinkrate=%s"%(str(shrinkrate)))
    logger.debug("lengthAlignmentShrinked=%d"%(lengthAlignmentShrinked))


    newAnnoList = []
    for i in range(len(topoSeqList)):
        try:
            anno = myfunc.GetFirstWord(annotationList[i])
        except:
            anno = ""
        anno = anno[:g_params['MAX_SIZE_ANNOTATION']]
        newAnnoList.append("%s %s"%(anno, tagList[i]))
    maxSizeAnno = max([len(x) for x in newAnnoList])

    logger.debug("maxSizeAnno=%d"%(maxSizeAnno))
    fonttype = 'monospace'

    numSeq = len(topoSeqList)

    aaSeqList = [] # list of amino acid sequences, aligned and shrinked if enabled
    final2seq_idxMapList = [] # the index map from the final (shrinked or not)
                              # sequence to the orignal unaligned sequence.
    krSeqList = []  # seqlist for positively charged residues, KR
    for i in range(numSeq):
        seqID = idList[i]
        idxmap_aligne2seq = lcmp.GetAlign2SeqMap(alignedTopoSeqList[i],
                alignedTopoSeqList[i].replace(GAP,""))  
        idxmap = {}
        if seqID in aaSeqDict:
            aaseq = aaSeqDict[seqID]
            #print aaseq
            #print alignedTopoSeqList[i]
            # alignedTopoSeqList is the original (non-shinked MSA)
            aaseq = MatchToAlignedSeq(aaseq, alignedTopoSeqList[i], seqID)
            tmpaaseq = ""
            tmpKRseq = ""
            for pp in range(lengthAlignment):
                if posindexmap != {}:
                    jseq = posindexmap[pp]
                else:
                    jseq = pp
                idxmap[pp] = idxmap_aligne2seq[jseq]
                aa = aaseq[jseq]
                tmpaaseq += (aa)
                if isDrawKRBias:
                    if (aa in ["K", "R"] and not
                            IsOutofMaxDistKR(posTMList[i], jseq,
                                maxDistKR)):
                        tmpKRseq += aa
                        if g_params['isPrintDebugInfo']:
                            print(seqID, aa, jseq, posTMList[i])
                    else:
                        tmpKRseq += " "
            aaSeqList.append(tmpaaseq)
            krSeqList.append(tmpKRseq)
        else:
            aaSeqList.append("")
            krSeqList.append("")
        final2seq_idxMapList.append(idxmap)

    # debug
    if g_params['isPrintDebugInfo'] and g_params['isDrawKRBias']:
        print("print shrinked toposeq and krseq")
        for k in range(numSeq):
            print("%s\t%s"%(idList[k], topoSeqList[k]))
            print("%s\t%s"%(idList[k], krSeqList[k]))
        sys.stdout.flush()
# setting font properties
    #ffam = "monospace"
    #ffam = "Fixed"
    ffam = "Courier New"
    fontpath = "%s/%s"%(g_params['font_dir'], "Courier_New.ttf")
    fontsize = 18
    fp = matplotlib.font_manager.FontProperties(
        fname=fontpath, family=ffam, style='normal', size=fontsize,
        weight='normal', stretch='normal')

    ffam = "Arial"
    fontpath = "%s/%s"%(g_params['font_dir'], "Arial.ttf")
    fp_anno = matplotlib.font_manager.FontProperties(
        fname=fontpath, family=ffam, style='normal', size=20,
        weight='normal', stretch='normal')
    
# get the text width and height in pixels
    x=0
    y=0
    linespaceInPixel = 36
    ss = "M"*1
    pth = matplotlib.textpath.TextPath((x, y), ss, prop=fp)
    bb = pth.get_extents(transform=None)
    widthSingleCharInPixel =  float(bb.width)/len(ss)
    heightSingleCharInPixel = float(bb.height)
    logger.debug("charwidth=%d, charheight=%d"%( widthSingleCharInPixel, heightSingleCharInPixel))



    widthAnnoInPixel = 0
    for i in range(numSeq):
        ss = newAnnoList[i]+ "M"*2 # leave 3 letters' space
        pth = matplotlib.textpath.TextPath((x, y), ss, prop=fp_anno)
        bb = pth.get_extents(transform=None)
        wtd =  float(bb.width)/len(ss.replace(" ",""))*len(ss)
        if wtd > widthAnnoInPixel:
            widthAnnoInPixel = wtd
    logger.debug( "widthAnnoInPixel=%d"%(widthAnnoInPixel))


    sumTextHeightInPixel = (heightSingleCharInPixel + linespaceInPixel)*(numSeq+1)
    sumTextWidthInPixel = (widthSingleCharInPixel)*(lengthAlignmentShrinked+1)
#    sumTextWidthAnnotationInPixel = widthSingleCharInPixel*(maxSizeAnno+5)
    sumTextWidthAnnotationInPixel = widthAnnoInPixel

    logger.debug("lengthAlignment=%d"% lengthAlignment)
    logger.debug("sumTextWidthAnnotationInPixel=%d"% sumTextWidthAnnotationInPixel)
    logger.debug("sumTextWidthInPixel=%d"% sumTextWidthInPixel)
    logger.debug("sumTextHeightInPixel=%d"% sumTextHeightInPixel)

# set aspect ratio
    if g_params['isDrawDGprofile']:
        heightRatios = [numSeq, 5]
        gs = gridspec.GridSpec(2, 1, height_ratios=heightRatios) 
    else:
        heightRatios = [1]

    logger.debug( "heightRatios=%s"%str( heightRatios))


    widthUnitFigureInPixel = 8*80
    heightUnitFigureInPixel = 6*80

    adjust_left = float(maxSizeAnno+5)/lengthAlignmentShrinked
    adjust_right = 0.99
    #adjust_left = (adjust_right) * (sumTextWidthAnnotationInPixel/(sumTextWidthAnnotationInPixel+lengthAlignment*widthSingleCharInPixel))
    adjust_top = max(1.0 - float(2)/numSeq, 0.7)
    adjust_bottom = min(float(2)/numSeq,0.3)
    logger.debug( "adjust_left=%d"%adjust_left)
    logger.debug( "adjust_right=%d"%adjust_right)
    logger.debug( "adjust_top=%d"%adjust_top)
    logger.debug( "adjust_bottom=%d"%adjust_bottom)

    subplot1_width_ratio = (adjust_right-adjust_left)
    subplot1_height_ratio = float(heightRatios[0])/sum(heightRatios)*(adjust_top-adjust_bottom)
#subplot1_width_ratio = 1.0/(1.0+0.2+0.2+adjust_left)
    widthUnitSubplot1InPixel = widthUnitFigureInPixel*subplot1_width_ratio
    heightUnitSubplot1InPixel = heightUnitFigureInPixel*subplot1_height_ratio

    widthscale = float(sumTextWidthInPixel)/widthUnitSubplot1InPixel+0.00
    heightscale = float(sumTextHeightInPixel)/heightUnitSubplot1InPixel+0.02
    logger.debug( "sumTextWidthInPixel=%d, sumTextHeightInPixel=%d"% (sumTextWidthInPixel, sumTextHeightInPixel))
    logger.debug( "widthscale=%s, heightscale=%s"%(str(widthscale), str(heightscale)))

    widthSubplot1InPixel = widthUnitSubplot1InPixel * widthscale
    heightSubplot1InPixel = heightUnitSubplot1InPixel * heightscale

    logger.debug( "widthSubplot1InPixel=%d"%widthSubplot1InPixel)
    logger.debug( "heightSubplot1InPixel=%d"% heightSubplot1InPixel)

    widthSingleCharInAxes = float(widthSingleCharInPixel)/widthSubplot1InPixel
    heightSingleCharInAxes = float(heightSingleCharInPixel)/heightSubplot1InPixel
    widthAnnotationInAxes = float(sumTextWidthAnnotationInPixel)/widthSubplot1InPixel
    linespaceInAxes = float(linespaceInPixel)/heightSubplot1InPixel

    logger.debug("widthSingleCharInAxes=%d, heightSingleCharInAxes=%d"%(
        widthSingleCharInAxes, heightSingleCharInAxes))

    fontsize_tics = 18
    fontsize_label = 24
    
#     fontsize_tics = 18/3
#     fontsize_label = 24/3

# create figure object
    figsize = (8*widthscale, 6*heightscale) # fig size in inches (width,height)
    fig = plt.figure(figsize = figsize) # set the figsize
    fig.subplots_adjust(left=adjust_left, right=adjust_right, top=adjust_top, bottom=adjust_bottom)

    plt.rc('legend',**{'fontsize':fontsize_label})

    if g_params['isDrawDGprofile']:
        ax = fig.add_subplot(gs[0])
    else:
        ax = fig.add_subplot(111)

    #ax.axis('off') 

    inv = ax.transAxes.inverted()

    loc_xtics = []
    label_xtics = []
    if posindexmap != {}:
        for jj in range(0, lengthAlignmentShrinked, 10):
            loc_xtics.append(jj)
            label_xtics.append(posindexmap[int(jj*shrinkrate)])
    else:
        for jj in range(0, lengthAlignmentShrinked, 10):
            loc_xtics.append(jj)
            label_xtics.append(int(jj*shrinkrate))

    ax.set_xlim(0, lengthAlignmentShrinked)
    ax.xaxis.set_visible(False)
    #ax.set_xlabel("Sequence position", fontsize=16)
    plt.xticks(np.array(loc_xtics), np.array(label_xtics))
    ax2 = ax.twiny()
    #ax2.set_xlabel("Alignment position", fontsize=16)
    ax2.set_xlim(0,lengthAlignmentShrinked)
    plt.xticks(np.array(loc_xtics), np.array(label_xtics), fontsize=fontsize_tics)
    plt.tick_params(labelsize=fontsize_tics, direction='out', pad=10)
    ax.set_ylim(numSeq,0)
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.spines['right'].set_visible(True)

# # make a plot of sequence indexes
#     l1 = []
#     l2 = []
#     for j in xrange(lengthAlignment):
#         if posindexmap != {}:
#             idxAlignedSeq = posindexmap[j]
#         else:
#             idxAlignedSeq = j
#         l1.append(idxAlignedSeq)
#         l2.append(0)
#     plt.plot(l1,l2, ' ')

    x0 = 0
    y0 = 1.0 - linespaceInAxes - heightSingleCharInAxes

    IOSTATE_LIST = ["i","o"]
    row_height = heightSingleCharInAxes + linespaceInAxes


# plot symbol annotation above the alignment
# blue upper hyphen Outside loop
# red  under hyphen Inside loop
# grey box TM helix (In -> out)
# white box TM helix (Out -> In)

    line1 = Line2D(list(range(1)), list(range(1)), color="red", marker='', markersize=5, markerfacecolor="red")
    line2 = Line2D(list(range(1)), list(range(1)), color="blue", marker='', markersize=5, markerfacecolor="blue")
    line3 = Line2D(list(range(1)), list(range(1)), color="white", marker='s', markersize=20, markerfacecolor="grey")
    line4 = Line2D(list(range(1)), list(range(1)), color="white", marker='s', markersize=20,markerfacecolor="white")
    legend = plt.legend((line1,line2,line3,line4),('Inside loop','Outside loop', 'TM helix (In -> Out)',
        'TM helix (Out -> In)'),numpoints=1, loc='upper center', bbox_to_anchor=(0.5,
            2.0), ncol=4, fancybox=False, shadow=False)
    legend.draw_frame(False)


    for i in range(len(topoSeqList)):
        # write sequence description
        anno = "%-*s"%(maxSizeAnno+5, newAnnoList[i])
        x = x0 - widthAnnotationInAxes
        y = y0 - i*row_height
        plt.text(x, y, anno, fontproperties=fp_anno, transform=ax.transAxes)

        if isDrawKRBias:
            seq = krSeqList[i]
        else:
            seq = ""
#         if aaSeqList[i] != "":
#             seq = aaSeqList[i]
#         else:
#             seq = topoSeqList[i]
        topo = topoSeqList[i]
        posTM = myfunc.GetTMPosition(topo)
        NtermState = lcmp.GetNtermState(topo)
        if len(posTM) == 0: # if non TM protein, just draw the sequence if specified
            x = x0 
            y = y0 - row_height*i
            txt = seq
            plt.text(x, y, txt, fontproperties=fp, transform=ax.transAxes)
        else: # draw TM regions and loops
            # terminal gaps are ignored
            li = []
            for (b, e) in posTM: # get list of segment borders
                li.append(b)
                li.append(e)
            for j in range(len(li)+1):
                # for TM helices, j-1 is even, j is odd
                if j == 0:
                    # ignore the N-terminal gaps
                    m = re.search('^-*', topo)
                    begin = len(m.group(0))
                else:
                    begin = li[j-1]

                if j != len(li):
                    end = li[j]
                else:
                    # ignore the C-terminal gaps
                    m = re.search('-*$', topo)
                    end = lengthAlignment - len(m.group(0))

                if isDrawKRBias: # draw positively charged K and R if enabled
                    y = y0 - row_height*i
                    for jpos in range(begin, end):
                        char = seq[jpos]
                        if char in ["K", "R"]:
                            x = x0 + jpos*widthSingleCharInAxes/shrinkrate
                            plt.text(x, y, char, fontproperties=fp,
                                    transform=ax.transAxes)

                txt_topo = topo[begin:end].replace(GAP," ")
                type_topo_stat = ""  #the state can be [IN, OUT, TM_IN_OUT, TM_OUT_IN]
                # setting properties for loops and TM regions
                if txt_topo.find('M')!=-1:
                    if lcmp.Get_IOState_upstream(topo, li[j-1]) == 'i':
                        type_topo_stat = "TM_IN_OUT"
                        edgecolor = 'black'
                        facecolor = 'grey'
                    else:
                        type_topo_stat = "TM_OUT_IN"
                        edgecolor = 'black'
                        facecolor = 'white'
                elif txt_topo.find('i') != -1: #inside
                    type_topo_stat = "IN"
                    color = 'red'
                elif txt_topo.find('o') != -1: #outside
                    type_topo_stat = "OUT"
                    color = 'blue'
                else:
                    facecolor = 'none'
                    edgecolor = 'none'
                width = widthSingleCharInAxes * (end-begin)/shrinkrate
                height = heightSingleCharInAxes + linespaceInAxes/2.0
                if type_topo_stat.find("TM") != -1: # draw TM regions
                    x = x0 + begin*widthSingleCharInAxes/shrinkrate
                    y = y0 - row_height*i - linespaceInAxes/4.0
                    rec = matplotlib.patches.Rectangle((x, y), width,
                            height, facecolor=facecolor,
                            edgecolor=edgecolor, transform=ax.transAxes)
                    ax.add_patch(rec)
                elif type_topo_stat in ["IN", "OUT"]: # draw loops
                    if type_topo_stat == "IN":
                        x1 = x0 + begin*widthSingleCharInAxes/shrinkrate
                        y1 = y0 - row_height*i - linespaceInAxes/4.0
                        x2 = x1 + width
                        y2 = y1
                    else: #OUT
                        x1 = x0 + begin*widthSingleCharInAxes/shrinkrate
                        y1 = y0 - row_height*i + heightSingleCharInAxes + linespaceInAxes/4.0
                        x2 = x1 + width
                        y2 = y1
                    ax.plot([x1, x2], [y1, y2], color=color, linestyle='-',
                            linewidth=2, transform=ax.transAxes)

    if g_params['isDrawDGprofile']:
        dgprofileDict = {} #{{{
        if os.path.exists(g_params['DGProfileFile']):
            dgprofileDict = ReadInDGProfile(g_params['DGProfileFile'])
        for i in range(numSeq):
            seqID = idList[i]
            toposeq = topoSeqList[i]
            lengthAlignment = len(toposeq)
            if (not seqID in dgprofileDict) and (seqID in aaSeqDict):
                aaseq = aaSeqDict[seqID]
                #print "aaseq=", aaseq
                tmpaaseqfile = tempfile.mktemp()
                tmpdgpfile = tempfile.mktemp()
                tmpfp = open(tmpaaseqfile, 'w')
                tmpfp.write(">%s\n"%seqID)
                tmpfp.write("%s\n"%aaseq)
                tmpfp.close()
                cmd =  [g_params['dgscanProg'], tmpaaseqfile,  "-lmin", "21", "-lmax",
                        "21", "-o", tmpdgpfile]
                cmdline = " ".join(cmd)
                print(cmdline)
                try:
                    subprocess.check_output(cmd)
                except subprocess.CalledProcessError as e:
                    print(e)
                tmpDGPDict = ReadInDGProfile(tmpdgpfile)
                os.remove(tmpaaseqfile)
                os.remove(tmpdgpfile)
                for seqid in tmpDGPDict:
                    dgprofileDict[seqid] = tmpDGPDict[seqid]
                if g_params['isPrintDebugInfo']:
                    print(dgprofileDict)
        #}}}
        ax = fig.add_subplot(gs[1])
        #ax.set_xlim(0, lengthAlignmentOriginal)
        ax.set_xlim(0, lengthAlignmentShrinked)
        ax.set_xlabel("Alignment position", fontsize=fontsize_label, labelpad=20)
        ax.set_ylabel(r"${\Delta}G$ (kJ/mol)", fontsize=fontsize_label)
        plt.xticks(np.array(loc_xtics), np.array(label_xtics), fontsize=fontsize_tics)
        plt.tick_params(labelsize=fontsize_tics, direction='out', pad=10)
        for i in range(numSeq):
            seqid = idList[i]
            alignedTopoOriginal = alignedTopoSeqList[i]
            align2seqMap = lcmp.GetAlign2SeqMap(alignedTopoOriginal,
                    alignedTopoOriginal.replace(GAP,""))
            #print "align2seqMap=", align2seqMap
            print("posindexmap=", posindexmap)
            #print "dgprofileDict[%s]"%seqid, dgprofileDict[seqid]
            try: 
                dgp = dgprofileDict[seqid]
                dt = {}
                for tup in dgp:
                    dt[tup[0]] = tup[1]

                x = []
                y = []
                for j in range(lengthAlignmentShrinked):
                    if posindexmap != {}:
                        try:
                            idxAlignedSeq = posindexmap[int(j*shrinkrate+0.5)]
                        except KeyError:
                            print("j=%d not in posindexmap"%j)
                            pass
                    else:
                        idxAlignedSeq = int(j*shrinkrate)
                    try:
                        idxSeq = align2seqMap[idxAlignedSeq]
                    except KeyError:
                        #print "idxAlignedSeq=%d not in align2seqMap"%idxAlignedSeq
                        pass

                    if idxSeq in dt:
                        x.append(j)
                        y.append(dt[idxSeq])
                    else:
#                         print "idxSeq=%d not in dgp, idxAlignedSeq=%d"%(idxSeq, idxAlignedSeq)
                        pass
                if i < len(colorList_DG_profile):
                    color = colorList_DG_profile[i]
                else:
                    color = 'none'
                if g_params['isPrintDebugInfo']:
                    print("DG-x:",x)
                    print("DG-y:",y)
# plot by line
                plt.plot(x,y, label=seqid, color=color)
# plot with '+' symbol
                # plt.plot(x,y, '+', label=seqid, color=color)
                plt.hlines(0, 0, lengthAlignmentOriginal)
                plt.legend()
            except KeyError:
                print("no dgprofile for %s"%(seqid))
                pass

    plt.savefig(pdffile)
    print("%s output"%(pdffile))
    cmd = "pdfcrop --margins '%d %d %d %d' --clip %s"%(pdfcrop_margin_left,
            pdfcrop_margin_top, pdfcrop_margin_right, pdfcrop_margin_bottom,
            pdffile)
    os.system(cmd)
    pdf_cropfile = os.path.splitext(pdffile)[0]+"-crop.pdf"
    pngfile = os.path.splitext(pdffile)[0] + "-crop.png"
    thumb_pngfile = os.path.splitext(pdffile)[0] + "-crop.thumb.png"
    cmd = ["convert", pdffile, pngfile] 
    try:
        subprocess.check_output(cmd)
    except subprocess.CalledProcessError as e:
        print(e)

    cmd = ["convert", "-thumbnail", "100", pngfile, thumb_pngfile] 
    try:
        subprocess.check_output(cmd)
    except subprocess.CalledProcessError as e:
        print(e)

#   Write Txtformat alignment
    #print final2seq_idxMapList
    htmlfile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'html')
#     WriteTXTAlignment(idList, newAnnoList, topoSeqList, alignedTopoSeqList,
#             aaSeqList, final2seq_idxMapList, txtfile)
    if len(idList) == 2:
        #WriteHTMLAlignment2(idList, newAnnoList, topoSeqList,
        #        alignedTopoSeqList, aaSeqList, final2seq_idxMapList, htmlfile)
        tmpmapList = []
        for i in range(len(alignedTopoSeqList)):
            tmpmap = {}
            for j in range(len(alignedTopoSeqList[i])):
                tmpmap[j] = j
            tmpmapList.append(tmpmap)
        WriteHTMLAlignment2(idList, newAnnoList, alignedTopoSeqList,
                alignedTopoSeqList, aaSeqAlignList, tmpmapList, htmlfile)
    elif len(idList) > 2:
        WriteHTMLAlignment3(idList, newAnnoList, topoSeqList,
                alignedTopoSeqList, aaSeqList, final2seq_idxMapList, htmlfile)


#}}}
def GetAlignedRegion(annotationList, topoSeqList):#{{{
    """Get aligned region from the description line

return a list of tuples [(beg, end)], starting from 0

In the annotation line, starting and both indeces are included
    """
    li = []
    for i in range (len(annotationList)):
        aligned_region = myfunc.GetAlignedRegionFromAnnotation(annotationList[i], method=0)
        if aligned_region[0] == -1:
            logger.debug("seqNo %d, aligned_region bad format, desp=\"\""%(i+1, annotationList[i]))
            aligned_region = (0, len(topoSeqList[i]))
        li.append(aligned_region)
    return li
#}}}


def DrawMSATopo_MAT_Core_unalign_rainbow(inFile, g_params):#{{{
    """Draw topology alignment core unaligned with rainbow color

    """
    logger = logging.getLogger(__name__)
    logger.debug("method=%s"%(g_params['method']))
    (idList, annotationList, orig_topoSeqList) = myfunc.ReadFasta(inFile)
    lst_aligned_region = GetAlignedRegion(annotationList, orig_topoSeqList)
    numSeq = len(idList)
    if numSeq < 1:
        logger.debug("No sequence in the file %s. Ignore." %(inFile))
        return 1


    topoSeqList = []
    lst_unalignedNterm = []
    lst_unalignedCterm = []
    for i in range(numSeq):
        (beg, end) = lst_aligned_region[i]
        topoSeqList.append(orig_topoSeqList[i][beg:end])
        lst_unalignedNterm.append(orig_topoSeqList[i][0:beg])
        lst_unalignedCterm.append(orig_topoSeqList[i][end:len(orig_topoSeqList[i])])



    maxlen_unaligned_nterm = max([len(x) for x in lst_unalignedNterm])
    maxlen_unaligned_cterm = max([len(x) for x in lst_unalignedCterm])

    orig_posTMList = [ myfunc.GetTMPosition(topo) for topo in orig_topoSeqList]
    numTMList=[ len (posTM) for posTM in orig_posTMList]

    blue = Color("blue")
    red = Color("red")
    lst_color = [list(blue.range_to(red,numTM)) for numTM in numTMList]

    lengthAlignmentOriginal = len(topoSeqList[0])

    maxDistKR = g_params['maxDistKR']
    marginX = g_params['marginX']
    marginY = g_params['marginY']
    annoSeqInterval = g_params['annoSeqInterval']
    widthAnnotation = g_params['widthAnnotation']
    heightScaleBar = g_params['heightScaleBar']
    heightTMbox = g_params['heightTMbox']
    scaleSeprationLine = g_params['scaleSeprationLine']
    font_size_scalebar = g_params['font_size_scalebar']
    fntScaleBar = g_params['fntScaleBar']
    (fontWidthScaleBar, fontHeightScaleBar) = fntScaleBar.getsize("a")

    pdfcrop_margin_left = g_params['pdfcrop_margin_left']
    pdfcrop_margin_top = g_params['pdfcrop_margin_top']
    pdfcrop_margin_right = g_params['pdfcrop_margin_right']
    pdfcrop_margin_bottom = g_params['pdfcrop_margin_bottom']

    rootname = os.path.basename(os.path.splitext(inFile)[0])
    aaSeqDict = GetAASeqDict(inFile)

    if g_params['outpath'] == "":
        outpath = myfunc.my_dirname(inFile)
    else:
        outpath = g_params['outpath']

    isDrawKRBias = g_params['isDrawKRBias']
    str_krbias = ""
    if g_params['isDrawKRBias'] == True:
        str_krbias = ".krbias"
    svgfile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'svg')
    pdffile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'pdf')
    txtfile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'txtplot')

# posindexmap: map of the residue position to the original MSA
# e.g. pos[0] = 5 means the first residue is actually the 6th residue position
# in the original MSA
#   backup original aligned topoSeqList

    logger.info("before adding Ms")
    logger.info("topoSeqList")
    for i in range(numSeq):
        logger.debug("%s - numTM = %d"%(idList[i],len(myfunc.GetTMPosition(topoSeqList[i]))))

    # add 'M's if the terminal at TM helix of the aligned region is chopped
    # to be too short, topoSeqList is the aligned region
    # 1. for the N-terminal at the aligned region
    MIN_LENGTH_OF_TM_TO_DRAWN = 6
    posTMList = [myfunc.GetTMPosition(x) for x in topoSeqList]
    for posTM in posTMList:
        lst_length_of_first_TM_withgaps = []
        if len(posTM) > 0:
            lenFirstTM = posTM[0][1]-posTM[0][0]
        else:
            lenFirstTM = 21
        lst_length_of_first_TM_withgaps.append(lenFirstTM)

    min_length_of_first_TM_withgaps = min(lst_length_of_first_TM_withgaps)
    if min_length_of_first_TM_withgaps < MIN_LENGTH_OF_TM_TO_DRAWN: # if the TM is chopped
        num_Ms_to_add = MIN_LENGTH_OF_TM_TO_DRAWN - min_length_of_first_TM_withgaps
        newTopoSeqList = []
        for top in topoSeqList:
            if lcmp.Get_nongap_downstream(top, 0) == "M": 
                newTopoSeqList.append("M"*num_Ms_to_add+top)
            else:
                newTopoSeqList.append("-"*num_Ms_to_add+top)
        topoSeqList = newTopoSeqList

    # 2. dealing with the last TM in the aligned region
    posTMList = [myfunc.GetTMPosition(x) for x in topoSeqList]
    for posTM in posTMList:
        lst_length_of_last_TM_withgaps = []
        if len(posTM) > 0:
            lenLastTM = posTM[-1][1]-posTM[-1][0]
        else:
            lenLastTM = 21
        lst_length_of_last_TM_withgaps.append(lenLastTM)

    min_length_of_last_TM_withgaps = min(lst_length_of_last_TM_withgaps)
    if min_length_of_last_TM_withgaps < MIN_LENGTH_OF_TM_TO_DRAWN: # if the TM is chopped
        num_Ms_to_add = MIN_LENGTH_OF_TM_TO_DRAWN - min_length_of_last_TM_withgaps
        newTopoSeqList = []
        for top in topoSeqList:
            if lcmp.Get_nongap_uppstream(top, len(top)-1) == "M": 
                newTopoSeqList.append(top+ "M"*num_Ms_to_add)
            else:
                newTopoSeqList.append(top+ "-"*num_Ms_to_add)
        topoSeqList = newTopoSeqList

    posTMList = [myfunc.GetTMPosition(x) for x in topoSeqList]
    alignedTopoSeqList = topoSeqList


    logger.info("after adding Ms")
    logger.info("topoSeqList")
    for i in range(numSeq):
        logger.debug("%s - numTM = %d"%(idList[i],len(myfunc.GetTMPosition(topoSeqList[i]))))
#===========================================
    posindexmap = {}
    method_shrink = g_params['method_shrink']

    aaSeqAlignList = [] # list of original aligned aa seq list
    for i in range(numSeq):
        seqid = idList[i]
        try:
            aaseq = aaSeqDict[seqid]
            aaseq = MatchToAlignedSeq(aaseq, alignedTopoSeqList[i], seqid)
        except KeyError:
            aaseq = " "*lengthAlignmentOriginal
        aaSeqAlignList.append(aaseq)


    if g_params['isShrink']:
        if g_params['method_shrink'] == 0:
            posindexmap = ShrinkGapInMSA_0(idList, topoSeqList)
        elif g_params['method_shrink'] == 1:
            posindexmap = ShrinkGapInMSA_exclude_TMregion(idList, topoSeqList)
        elif g_params['method_shrink'] == 2:
            (idxmap_align2shrink, idxmap_shrink2align) =\
                    ShrinkMSA_Method_2(topoSeqList, aaSeqAlignList, posTMList,
                            g_params['shrinkrate_TM'], g_params['max_hold_loop'],
                            g_params['isDrawKRBias'])
            posindexmap = idxmap_shrink2align

    logger.info("222222222222222222222")
    logger.info("topoSeqList")
    for i in range(numSeq):
        logger.debug("%s - numTM = %d"%(idList[i],len(myfunc.GetTMPosition(topoSeqList[i]))))

    posTM = myfunc.GetTMPosition(topoSeqList[0])
    g_params['widthAnnotation'] = GetSizeAnnotationToDraw(annotationList)
    widthAnnotation = g_params['widthAnnotation']
    tagList = []
    for seqAnno in annotationList:
        tagList.append(GetSeqTag(seqAnno))
#     print tagList
    numSeprationLine = len(set(tagList))
    lengthAlignment = len(topoSeqList[0])

# added 2013-12-04
    if g_params['shrinkrate'] != None:
        shrinkrate = g_params['shrinkrate'] # shrink the sequence proportionally
    else:
        shrinkrate = lengthAlignment/120.0
    lengthAlignmentShrinked = int(lengthAlignment/shrinkrate+0.5)
    logger.debug("lengthAlignment=%d"%lengthAlignment)
    logger.debug("shrinkrate=%s"%(str(shrinkrate)))
    logger.debug("lengthAlignmentShrinked=%d"%(lengthAlignmentShrinked))

    maxlen_unaligned_cterm_shrinked = int(maxlen_unaligned_cterm/shrinkrate+0.5)
    maxlen_unaligned_nterm_shrinked = int(maxlen_unaligned_nterm/shrinkrate+0.5)

    newAnnoList = []
    for i in range(len(topoSeqList)):
        try:
            anno = myfunc.GetFirstWord(annotationList[i]).split('|')[0]
        except:
            anno = ""
        anno = anno[:g_params['MAX_SIZE_ANNOTATION']]
        newAnnoList.append("%s %s"%(anno, tagList[i]))
    maxSizeAnno = max([len(x) for x in newAnnoList])

    logger.debug("maxSizeAnno=%d"%(maxSizeAnno))
    fonttype = 'monospace'

    numSeq = len(topoSeqList)

    logger.info("after shrinking")
    logger.info("topoSeqList")
    for i in range(numSeq):
        logger.debug("%s - numTM = %d"%(idList[i],len(myfunc.GetTMPosition(topoSeqList[i]))))

    aaSeqList = [] # list of amino acid sequences, aligned and shrinked if enabled
    final2seq_idxMapList = [] # the index map from the final (shrinked or not)
                              # sequence to the orignal unaligned sequence.
    krSeqList = []  # seqlist for positively charged residues, KR
    for i in range(numSeq):
        seqID = idList[i]
        idxmap_aligne2seq = lcmp.GetAlign2SeqMap(alignedTopoSeqList[i],
                alignedTopoSeqList[i].replace(GAP,""))  
        idxmap = {}
        if seqID in aaSeqDict:
            aaseq = aaSeqDict[seqID]
            #print aaseq
            #print alignedTopoSeqList[i]
            # alignedTopoSeqList is the original (non-shinked MSA)
            aaseq = MatchToAlignedSeq(aaseq, alignedTopoSeqList[i], seqID)
            tmpaaseq = ""
            tmpKRseq = ""
            for pp in range(lengthAlignment):
                if posindexmap != {}:
                    jseq = posindexmap[pp]
                else:
                    jseq = pp
                idxmap[pp] = idxmap_aligne2seq[jseq]
                aa = aaseq[jseq]
                tmpaaseq += (aa)
                if isDrawKRBias:
                    if (aa in ["K", "R"] and not
                            IsOutofMaxDistKR(posTMList[i], jseq,
                                maxDistKR)):
                        tmpKRseq += aa
                        if g_params['isPrintDebugInfo']:
                            print(seqID, aa, jseq, posTMList[i])
                    else:
                        tmpKRseq += " "
            aaSeqList.append(tmpaaseq)
            krSeqList.append(tmpKRseq)
        else:
            aaSeqList.append("")
            krSeqList.append("")
        final2seq_idxMapList.append(idxmap)

    # debug
    if g_params['isPrintDebugInfo']:
        print("print shrinked toposeq and krseq")
        for k in range(numSeq):
            print("%s\t%s"%(idList[k], topoSeqList[k]))
            print("%s\t%s"%(idList[k], krSeqList[k]))
        sys.stdout.flush()
# setting font properties
    #ffam = "monospace"
    #ffam = "Fixed"
    ffam = "Courier New"
    fontpath = "%s/%s"%(g_params['font_dir'], "Courier_New.ttf")
    fontsize = 18
    fp = matplotlib.font_manager.FontProperties(
        fname=fontpath, family=ffam, style='normal', size=fontsize,
        weight='normal', stretch='normal')

    ffam = "Arial"
    fontpath = "%s/%s"%(g_params['font_dir'], "Arial.ttf")
    fp_anno = matplotlib.font_manager.FontProperties(
        fname=fontpath, family=ffam, style='normal', size=20,
        weight='normal', stretch='normal')
    
# get the text width and height in pixels
    x=0
    y=0
    linespaceInPixel = 36
    ss = "M"*1
    pth = matplotlib.textpath.TextPath((x, y), ss, prop=fp)
    bb = pth.get_extents(transform=None)
    widthSingleCharInPixel =  float(bb.width)/len(ss)
    heightSingleCharInPixel = float(bb.height)
    logger.debug("charwidth=%d, charheight=%d"%( widthSingleCharInPixel, heightSingleCharInPixel))



    widthAnnoInPixel = 0
    for i in range(numSeq):
        ss = newAnnoList[i]+ "M"*2 # leave 3 letters' space
        pth = matplotlib.textpath.TextPath((x, y), ss, prop=fp_anno)
        bb = pth.get_extents(transform=None)
        wtd =  float(bb.width)/len(ss.replace(" ",""))*len(ss)
        if wtd > widthAnnoInPixel:
            widthAnnoInPixel = wtd
    logger.debug( "widthAnnoInPixel=%d"%(widthAnnoInPixel))


    sumTextHeightInPixel = (heightSingleCharInPixel + linespaceInPixel)*(numSeq+1)
    sumTextWidthInPixel = (widthSingleCharInPixel)*(lengthAlignmentShrinked+maxlen_unaligned_nterm_shrinked+maxlen_unaligned_cterm_shrinked+1)
    #sumTextWidthInPixel = (widthSingleCharInPixel)*(lengthAlignmentShrinked+maxlen_unaligned_cterm_shrinked+1)
    sumTextWidthAnnotationInPixel = widthAnnoInPixel

    logger.debug("lengthAlignment=%d"% lengthAlignment)
    logger.debug("sumTextWidthAnnotationInPixel=%d"% sumTextWidthAnnotationInPixel)
    logger.debug("sumTextWidthInPixel=%d"% sumTextWidthInPixel)
    logger.debug("sumTextHeightInPixel=%d"% sumTextHeightInPixel)

# set aspect ratio
    if g_params['isDrawDGprofile']:
        heightRatios = [numSeq, 5]
        gs = gridspec.GridSpec(2, 1, height_ratios=heightRatios) 
    else:
        heightRatios = [1]

    logger.debug( "heightRatios=%s"%str( heightRatios))


    widthUnitFigureInPixel = 8*80
    heightUnitFigureInPixel = 6*80

    adjust_left = float(maxSizeAnno+5)/lengthAlignmentShrinked
    adjust_right = 0.99
    adjust_top = max(1.0 - float(2)/numSeq, 0.65)
    adjust_bottom = min(float(2)/numSeq,0.3)
    logger.debug( "adjust_left=%d"%adjust_left)
    logger.debug( "adjust_right=%d"%adjust_right)
    logger.debug( "adjust_top=%d"%adjust_top)
    logger.debug( "adjust_bottom=%d"%adjust_bottom)

    subplot1_width_ratio = (adjust_right-adjust_left)
    subplot1_height_ratio = float(heightRatios[0])/sum(heightRatios)*(adjust_top-adjust_bottom)
#subplot1_width_ratio = 1.0/(1.0+0.2+0.2+adjust_left)
    widthUnitSubplot1InPixel = widthUnitFigureInPixel*subplot1_width_ratio
    heightUnitSubplot1InPixel = heightUnitFigureInPixel*subplot1_height_ratio

    widthscale = float(sumTextWidthInPixel)/widthUnitSubplot1InPixel+0.00
    heightscale = float(sumTextHeightInPixel)/heightUnitSubplot1InPixel+0.02
    logger.debug( "sumTextWidthInPixel=%d, sumTextHeightInPixel=%d"% (sumTextWidthInPixel, sumTextHeightInPixel))
    logger.debug( "widthscale=%s, heightscale=%s"%(str(widthscale), str(heightscale)))

    widthSubplot1InPixel = widthUnitSubplot1InPixel * widthscale
    heightSubplot1InPixel = heightUnitSubplot1InPixel * heightscale

    logger.debug( "widthSubplot1InPixel=%d"%widthSubplot1InPixel)
    logger.debug( "heightSubplot1InPixel=%d"% heightSubplot1InPixel)

    widthSingleCharInAxes = float(widthSingleCharInPixel)/widthSubplot1InPixel
    heightSingleCharInAxes = float(heightSingleCharInPixel)/heightSubplot1InPixel
    widthAnnotationInAxes = float(sumTextWidthAnnotationInPixel)/widthSubplot1InPixel
    linespaceInAxes = float(linespaceInPixel)/heightSubplot1InPixel

    logger.debug("widthSingleCharInAxes=%d, heightSingleCharInAxes=%d"%(
        widthSingleCharInAxes, heightSingleCharInAxes))

    fontsize_tics = 18
    fontsize_label = 24

#     fontsize_tics = 18/3
#     fontsize_label = 24/3

# create figure object
    figsize = (8*widthscale, 6*heightscale) # fig size in inches (width,height)
    fig = plt.figure(figsize = figsize) # set the figsize
    fig.subplots_adjust(left=adjust_left, right=adjust_right, top=adjust_top, bottom=adjust_bottom)

    plt.rc('legend',**{'fontsize':fontsize_label})

    if g_params['isDrawDGprofile']:
        ax = fig.add_subplot(gs[0])
    else:
        ax = fig.add_subplot(111)

    #ax.axis('off') 

    inv = ax.transAxes.inverted()

# label_xtics is just for the aligned region

    loc_xtics = []
    label_xtics = []
    if posindexmap != {}:
        for jj in range(0, lengthAlignmentShrinked, 10):
            loc_xtics.append(jj)
            label_xtics.append(posindexmap[int(jj*shrinkrate)])
    else:
        for jj in range(0, lengthAlignmentShrinked, 10):
            loc_xtics.append(jj)
            label_xtics.append(int(jj*shrinkrate))

    ax.set_xlim(0, lengthAlignmentShrinked)
    ax.xaxis.set_visible(False)
    #ax.set_xlabel("Sequence position", fontsize=16)
    #plt.xticks(np.array(loc_xtics), np.array(label_xtics))
    #ax2 = ax.twiny()
    #ax2.set_xlabel("Alignment position", fontsize=16)
    #ax2.set_xlim(0,lengthAlignmentShrinked)
    #plt.xticks(np.array(loc_xtics), np.array(label_xtics), fontsize=fontsize_tics)
    plt.tick_params(labelsize=fontsize_tics, direction='out', pad=10)
    ax.set_ylim(numSeq,0)
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)


    x0 = 0
    y0 = 1.0 - linespaceInAxes - heightSingleCharInAxes

    IOSTATE_LIST = ["i","o"]
    row_height = heightSingleCharInAxes + linespaceInAxes


# plot symbol annotation above the alignment
# blue upper hyphen Outside loop
# red  under hyphen Inside loop
# grey box TM helix (In -> out)
# white box TM helix (Out -> In)

    line1 = Line2D(list(range(1)), list(range(1)), color="red", marker='', markersize=5, markerfacecolor="red")
    line2 = Line2D(list(range(1)), list(range(1)), color="blue", marker='', markersize=5, markerfacecolor="blue")
    #line3 = Line2D(range(1), range(1), color="white", marker='s', markersize=20, markerfacecolor="grey")
    #line4 = Line2D(range(1), range(1), color="white", marker='s', markersize=20,markerfacecolor="white")
#     legend = plt.legend((line1,line2,line3,line4),('Inside loop','Outside loop', 'TM helix (In -> Out)',
#         'TM helix (Out -> In)'),numpoints=1, loc='upper center', bbox_to_anchor=(0.5,
#             2.0), ncol=4, fancybox=False, shadow=False)
    legend = plt.legend((line1,line2),('Inside loop','Outside loop'),
        numpoints=1, loc='upper center', bbox_to_anchor=(0.5,
            1.5), ncol=2, fancybox=False, shadow=False)
    legend.draw_frame(True)


#=================================================
#draw topology of each sequence
#=================================================
    # draw a transparent box for the aligned region
    shiftx = (widthSingleCharInAxes*maxlen_unaligned_nterm)/shrinkrate
    x = x0 - widthSingleCharInAxes/shrinkrate + shiftx
    y = y0 - row_height*(numSeq) + linespaceInAxes*0.625
    width = widthSingleCharInAxes * (lengthAlignment+2)/shrinkrate
    height = row_height*(numSeq+1) - linespaceInAxes
    facecolor = "#CCCCEF"
    edgecolor = 'black'
    rec = matplotlib.patches.Rectangle((x, y), width,
            height, alpha=0.6, linewidth=3, linestyle="solid", facecolor=facecolor,
            edgecolor=edgecolor, transform=ax.transAxes)
    ax.add_patch(rec)
    for i in range(numSeq):
        # write sequence description
        anno = "%-*s"%(maxSizeAnno+5, newAnnoList[i])
        x = x0 - widthAnnotationInAxes
        y = y0 - i*row_height
        plt.text(x, y, anno, fontproperties=fp_anno, transform=ax.transAxes)

        if isDrawKRBias:
            seq = krSeqList[i]
        else:
            seq = ""
        #------------------------------
        # Draw aligned region and unaligned C-terminal
        #------------------------------

        idxTM_drawn_set = set([]) # set of the TM index already been drawn,
                                  # avoid duplicated drawn when cross border
        for item in ["unaligned_N_term", "aligned_region", "unaligned_C_term"]:
        #for item in ["aligned_region"]:
        #for item in ["unaligned_C_term"]:
        #for item in ["unaligned_N_term"]:
            if item == "aligned_region":
                topo = topoSeqList[i]
                shiftx = (widthSingleCharInAxes*maxlen_unaligned_nterm)/shrinkrate
            elif item == "unaligned_C_term":
                topo = lst_unalignedCterm[i]
                shiftx = (widthSingleCharInAxes*(len(topoSeqList[0])+maxlen_unaligned_nterm))/shrinkrate
            elif item == "unaligned_N_term":
                topo = lst_unalignedNterm[i]
                shiftx = (widthSingleCharInAxes*(maxlen_unaligned_nterm-len(lst_unalignedNterm[i])))/shrinkrate


            posTM = myfunc.GetTMPosition(topo)
            logger.debug("%s - %s: %s"%(anno, item, topo))
            if len(posTM) == 0: # if non TM protein, just draw the sequence if specified
                x = x0 + shiftx
                y = y0 - row_height*i
                txt = seq
                plt.text(x, y, txt, fontproperties=fp, transform=ax.transAxes)
            else: # draw TM regions and loops
                # terminal gaps are ignored
                li = []
                for (b, e) in posTM: # get list of segment borders
                    li.append(b)
                    li.append(e)
                for j in range(len(li)+1):
                    # for TM helices, j-1 is even, j is odd
                    numTM_unalignedNterm = len(myfunc.GetTMPosition(lst_unalignedNterm[i]))
                    if item == "unaligned_N_term":
                        idxTM = j/2
                    elif item == "aligned_region":
                        idxTM = j/2 + numTM_unalignedNterm
                        if (len(lst_unalignedNterm[i]) > 0 
                                and (lst_unalignedNterm[i][-1] == "M" and topoSeqList[i][0]=="M")):
                            idxTM -= 1
                    elif item == "unaligned_C_term":
                        idxTM = j/2 + len(myfunc.GetTMPosition(lst_unalignedNterm[i]+topoSeqList[i]))
                        if topoSeqList[i][-1] == "M" and lst_unalignedCterm[i][0] == "M":
                            idxTM -= 1
                    try:
                        color = lst_color[i][idxTM].get_hex_l()
                    except:
                        color = "red"
#                     elif item in ["unaligned_C_term", "unaligned_N_term"]:
#                         color = "#EEEEEE"
                    if j == 0:
                        # ignore the N-terminal gaps
                        m = re.search('^-*', topo)
                        begin = len(m.group(0))
                    else:
                        begin = li[j-1]

                    if j != len(li):
                        end = li[j]
                    else:
                        # ignore the C-terminal gaps
                        m = re.search('-*$', topo)
                        end = len(topo) - len(m.group(0))

                    if isDrawKRBias: # draw positively charged K and R if enabled
                        y = y0 - row_height*i
                        for jpos in range(begin, end):
                            char = seq[jpos]
                            if char in ["K", "R"]:
                                x = x0 + jpos*widthSingleCharInAxes/shrinkrate
                                plt.text(x, y, char, fontproperties=fp,
                                        transform=ax.transAxes)

                    txt_topo = topo[begin:end].replace(GAP," ")
                    type_topo_stat = ""  #the state can be [IN, OUT, TM_IN_OUT, TM_OUT_IN]
                    # setting properties for loops and TM regions
                    if txt_topo.find('M')!=-1:
                        if lcmp.Get_IOState_upstream(topo, li[j-1]) == 'i':
                            type_topo_stat = "TM_IN_OUT"
                            edgecolor = 'black'
                            #facecolor = 'grey'
                            facecolor = color
                        else:
                            type_topo_stat = "TM_OUT_IN"
                            edgecolor = 'black'
                            #facecolor = 'white'
                            facecolor = color
                    elif txt_topo.find('i') != -1: #inside
                        type_topo_stat = "IN"
                        color = 'red'
                    elif txt_topo.find('o') != -1: #outside
                        type_topo_stat = "OUT"
                        color = 'blue'
                    else:
                        facecolor = 'none'
                        edgecolor = 'none'
                    width = widthSingleCharInAxes * (end-begin)/shrinkrate
                    height = heightSingleCharInAxes + linespaceInAxes/2.0
                    if type_topo_stat.find("TM") != -1: # draw TM regions
                        x = x0 + begin*widthSingleCharInAxes/shrinkrate + shiftx
                        y = y0 - row_height*i - linespaceInAxes/4.0
                        rec = matplotlib.patches.Rectangle((x, y), width,
                                height, facecolor=facecolor,
                                edgecolor=edgecolor, transform=ax.transAxes)
                        ax.add_patch(rec)

                        # draw TM indeces
                        if not idxTM in idxTM_drawn_set:
                            x = x + width/2.5
                            y = y + height/4.0
                            plt.text(x, y, "%d"%(idxTM+1), fontproperties=fp_anno, transform=ax.transAxes)
                            idxTM_drawn_set.add(idxTM)
                            logger.debug("idxTM=%d"%(idxTM))

                    # draw loops 
                    elif type_topo_stat in ["IN", "OUT"]: # draw loops
                        if type_topo_stat == "IN":
                            x1 = x0 + begin*widthSingleCharInAxes/shrinkrate + shiftx
                            y1 = y0 - row_height*i - linespaceInAxes/4.0
                            x2 = x1 + width
                            y2 = y1
                        else: #OUT
                            x1 = x0 + begin*widthSingleCharInAxes/shrinkrate + shiftx
                            y1 = y0 - row_height*i + heightSingleCharInAxes + linespaceInAxes/4.0
                            x2 = x1 + width
                            y2 = y1
                        ax.plot([x1, x2], [y1, y2], color=color, linestyle='-',
                                linewidth=2, transform=ax.transAxes)

    if g_params['isDrawDGprofile']:
        dgprofileDict = {} #{{{
        if os.path.exists(g_params['DGProfileFile']):
            dgprofileDict = ReadInDGProfile(g_params['DGProfileFile'])
        for i in range(numSeq):
            seqID = idList[i]
            toposeq = topoSeqList[i]
            lengthAlignment = len(toposeq)
            if (not seqID in dgprofileDict) and (seqID in aaSeqDict):
                aaseq = aaSeqDict[seqID]
                #print "aaseq=", aaseq
                tmpaaseqfile = tempfile.mktemp()
                tmpdgpfile = tempfile.mktemp()
                tmpfp = open(tmpaaseqfile, 'w')
                tmpfp.write(">%s\n"%seqID)
                tmpfp.write("%s\n"%aaseq)
                tmpfp.close()
                cmd =  [g_params['dgscanProg'], tmpaaseqfile,  "-lmin", "21", "-lmax",
                        "21", "-o", tmpdgpfile]
                cmdline = " ".join(cmd)
                print(cmdline)
                try:
                    subprocess.check_output(cmd)
                except subprocess.CalledProcessError as e:
                    print(e)
                tmpDGPDict = ReadInDGProfile(tmpdgpfile)
                os.remove(tmpaaseqfile)
                os.remove(tmpdgpfile)
                for seqid in tmpDGPDict:
                    dgprofileDict[seqid] = tmpDGPDict[seqid]
                if g_params['isPrintDebugInfo']:
                    print(dgprofileDict)
        #}}}
        ax = fig.add_subplot(gs[1])
        #ax.set_xlim(0, lengthAlignmentOriginal)
        ax.set_xlim(0, lengthAlignmentShrinked)
        ax.set_xlabel("Alignment position", fontsize=fontsize_label, labelpad=20)
        ax.set_ylabel(r"${\Delta}G$ (kJ/mol)", fontsize=fontsize_label)
        plt.xticks(np.array(loc_xtics), np.array(label_xtics), fontsize=fontsize_tics)
        plt.tick_params(labelsize=fontsize_tics, direction='out', pad=10)
        for i in range(numSeq):
            seqid = idList[i]
            alignedTopoOriginal = alignedTopoSeqList[i]
            align2seqMap = lcmp.GetAlign2SeqMap(alignedTopoOriginal,
                    alignedTopoOriginal.replace(GAP,""))
            #print "align2seqMap=", align2seqMap
            print("posindexmap=", posindexmap)
            #print "dgprofileDict[%s]"%seqid, dgprofileDict[seqid]
            try: 
                dgp = dgprofileDict[seqid]
                dt = {}
                for tup in dgp:
                    dt[tup[0]] = tup[1]

                x = []
                y = []
                for j in range(lengthAlignmentShrinked):
                    if posindexmap != {}:
                        try:
                            idxAlignedSeq = posindexmap[int(j*shrinkrate+0.5)]
                        except KeyError:
                            print("j=%d not in posindexmap"%j)
                            pass
                    else:
                        idxAlignedSeq = int(j*shrinkrate)
                    try:
                        idxSeq = align2seqMap[idxAlignedSeq]
                    except KeyError:
                        #print "idxAlignedSeq=%d not in align2seqMap"%idxAlignedSeq
                        pass

                    if idxSeq in dt:
                        x.append(j)
                        y.append(dt[idxSeq])
                    else:
#                         print "idxSeq=%d not in dgp, idxAlignedSeq=%d"%(idxSeq, idxAlignedSeq)
                        pass
                if i < len(colorList_DG_profile):
                    color = colorList_DG_profile[i]
                else:
                    color = 'none'
                if g_params['isPrintDebugInfo']:
                    print("DG-x:",x)
                    print("DG-y:",y)
# plot by line
                plt.plot(x,y, label=seqid, color=color)
# plot with '+' symbol
                # plt.plot(x,y, '+', label=seqid, color=color)
                plt.hlines(0, 0, lengthAlignmentOriginal)
                plt.legend()
            except KeyError:
                print("no dgprofile for %s"%(seqid))
                pass

    plt.savefig(pdffile)
    print("%s output"%(pdffile))
    cmd = "pdfcrop --margins '%d %d %d %d' --clip %s"%(pdfcrop_margin_left,
            pdfcrop_margin_top, pdfcrop_margin_right, pdfcrop_margin_bottom,
            pdffile)
    os.system(cmd)
    pdf_cropfile = os.path.splitext(pdffile)[0]+"-crop.pdf"
    pngfile = os.path.splitext(pdffile)[0] + "-crop.png"
    thumb_pngfile = os.path.splitext(pdffile)[0] + "-crop.thumb.png"
    cmd = ["convert", pdffile, pngfile] 
    try:
        subprocess.check_output(cmd)
    except subprocess.CalledProcessError as e:
        print(e)

    cmd = ["convert", "-thumbnail", "100", pngfile, thumb_pngfile] 
    try:
        subprocess.check_output(cmd)
    except subprocess.CalledProcessError as e:
        print(e)

#   Write Txtformat alignment
    #print final2seq_idxMapList
    if len(aaSeqList) == numSeq and len(aaSeqList[0]) == len(orig_topoSeqList[0]):
        htmlfile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'html')
        if len(idList) == 2:
            tmpmapList = []
            for i in range(len(alignedTopoSeqList)):
                tmpmap = {}
                for j in range(len(alignedTopoSeqList[i])):
                    tmpmap[j] = j
                tmpmapList.append(tmpmap)
            WriteHTMLAlignment2(idList, newAnnoList, alignedTopoSeqList,
                    alignedTopoSeqList, aaSeqAlignList, tmpmapList, htmlfile)
        elif len(idList) > 2:
            WriteHTMLAlignment3(idList, newAnnoList, topoSeqList,
                    alignedTopoSeqList, aaSeqList, final2seq_idxMapList, htmlfile)


#}}}
def DrawMSATopo_PYX(inFile, g_params):#{{{
    (idList, annotationList, topoSeqList) = myfunc.ReadFasta(inFile)
    topoSeqList = lcmp.RemoveUnnecessaryGap(topoSeqList)
    numSeq = len(idList)
    if numSeq < 1:
        print("No sequence in the file %s. Ignore." %(inFile), file=sys.stderr)
        return 1

    marginX = g_params['marginX']
    marginY = g_params['marginY']
    annoSeqInterval = g_params['annoSeqInterval']
    widthAnnotation = g_params['widthAnnotation']
    heightScaleBar = g_params['heightScaleBar']
    heightTMbox = g_params['heightTMbox']
    scaleSeprationLine = g_params['scaleSeprationLine']
    font_size_scalebar = g_params['font_size_scalebar']
    fntScaleBar = g_params['fntScaleBar']
    (fontWidthScaleBar, fontHeightScaleBar) = fntScaleBar.getsize("a")

    rootname = os.path.basename(os.path.splitext(inFile)[0])
    aaSeqDict = GetAASeqDict(inFile)

    #rootname=rootname.split('.')[0]
    if g_params['outpath'] == "":
        outpath = myfunc.my_dirname(inFile)
    else:
        outpath = g_params['outpath']

    str_krbias = ""
    if g_params['isDrawKRBias'] == True:
        str_krbias = ".krbias"

    pdffile = "%s%s%s%s.%s"%(outpath, os.sep, rootname, str_krbias, 'pdf')
    stemname = "%s%s%s%s"%(outpath, os.sep, rootname, str_krbias)

# posindexmap: map of the residue position to the original MSA
# e.g. pos[0] = 5 means the first residue is actually the 6th residue position
# in the original MSA
#   backup original aligned topoSeqList
    alignedTopoSeqList = []
    posTMList = []
    for seq in topoSeqList:
        alignedTopoSeqList.append(seq)
        posTMList.append(myfunc.GetTMPosition(seq))

    posindexmap = {}
    if g_params['isShrink']:
        posindexmap = ShrinkGapInMSA_0(idList, topoSeqList)

    posTM = myfunc.GetTMPosition(topoSeqList[0])
    g_params['widthAnnotation'] = GetSizeAnnotationToDraw(annotationList)
    widthAnnotation = g_params['widthAnnotation']
    tagList = []
    for seqAnno in annotationList:
        tagList.append(GetSeqTag(seqAnno))
#     print tagList
    numSeprationLine = len(set(tagList))
    lengthAlignment = len(topoSeqList[0])

    newAnnoList = []
    for i in range(len(topoSeqList)):
        newAnnoList.append("%s %s"%(idList[i], tagList[i]))
    maxSizeAnno = max([len(x) for x in newAnnoList])

    print("maxSizeAnno=",maxSizeAnno)
    fonttype = 'monospace'

    x0 = 10
    y0 = 10
    pyx.unit.set(xscale=1)
    for i in range(len(topoSeqList)):
        y = y0 + (i)*20
        ss = r"%-*s %s"%(maxSizeAnno+5, newAnnoList[i], topoSeqList[i])
        ss = ss.replace("_", "\\_")
        print(ss)
        #ss = r"Boxed text"
        #ss = r"%s %s"%("Boxed text new strings -- -- -----", " ".join(["one"]*20))
        #ss = r"Q81HG1 TM2SEQ                  --    -o ---oMMMMMMMMMMMMMMMMMMMMMiiiiiMMMMMMMMMMMMMMMMMMoMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiiii - -- --"
        #ss = r"Q63UX3 TM2GAPANDM2SEQ          -M MM-MMMMMMMiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMoMMMMMMM-MMMMMMMMMMMMMiMMMMMMMMMMMMMMMMMMMoo- - -"
        tbox = pyx.text.text(x0, y, ss)
        tpath = tbox.bbox().enlarged(3*pyx.unit.x_pt).path()
        c = pyx.canvas.canvas()
        c.draw(tpath, [pyx.deco.filled([pyx.color.cmyk.Yellow]), pyx.deco.stroked()])
        c.insert(tbox)

    c.writePDFfile(stemname)
    print("%s output"%(pdffile))

#}}}

def main(g_params):#{{{
    logger = logging.getLogger(__name__)
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    filelist = []
    filelistfile = ""
    aaSeqFile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:#{{{
        if isNonOptionArg == True:
            isNonOptionArg=False
            i = i + 1
        elif argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-i", "--infile"]:
                filelist.append(argv[i+1])
                i = i + 2
            elif argv[i] in ["-l", "--l", "-list", "--list"]:
                filelistfile, i = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-outpath",  "--outpath"]:
                (g_params['outpath'],i)  = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-method",  "--method"]:
                (g_params['method'],i)  = myfunc.my_getopt_str(argv, i)
            elif argv[i] in  ["-aapath", "--aapath"]:
                (g_params['aapath'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-of" , "--of", "--outformat"]:
                (g_params['outFormat'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-fontsize", "--fontsize"]:
                (g_params['font_size'], i) = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-mode", "--mode"]:
                (g_params['mode'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-log", "--log"]:
                (g_params['log_config_file'],i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-dgpfile", "--dgpfile"]:
                (g_params['DGProfileFile'],i) =  myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-text", "--text"]:
                (tmpstr, i) =  myfunc.my_getopt_str(argv, i)
                if (tmpstr.lower())[0] == "y": 
                    g_params['isDrawText'] = True
                else:
                    g_params['isDrawText'] = False
            elif argv[i] in ["-krbias", "--krbias"]:
                g_params['isDrawKRBias'] = True; i = i + 1
            elif argv[i] in ["-maxdistkr", "--maxdistkr"]:
                (g_params['maxDistKR'], i) = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-htmlheader", "--htmlheader"]:
                g_params['htmlheader'], i  = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-colorhtml", "--colorhtml"]:
                g_params['colorhtml'] = True; i += 1
            elif argv[i] in ["-win", "--win"]:
                g_params['window_size'], i = myfunc.my_getopt_int(argv, i)
            elif (argv[i] in ["-sep", "--sep"]):
                if (argv[i+1].lower())[0] == "y": 
                    g_params['isDrawSeprationLine'] = True
                else:
                    g_params['isDrawSeprationLine'] = False
                i = i + 2
            elif (argv[i] in ["-pfm", "--pfm"]):
                if (argv[i+1].lower())[0] == "y": 
                    g_params['isDrawPerMDistribution'] = True
                else:
                    g_params['isDrawPerMDistribution'] = False
                i = i + 2
            elif (argv[i] in ["-pdg", "--pdg"]):
                if (argv[i+1].lower())[0] == "y": 
                    g_params['isDrawDGprofile'] = True
                else:
                    g_params['isDrawDGprofile'] = False
                i = i + 2
            elif (argv[i] in ["-pmsa", "--pmsa"]):
                if (argv[i+1].lower())[0] == "y": 
                    g_params['isDrawMSA'] = True
                else:
                    g_params['isDrawMSA'] = False
                i = i + 2
            elif (argv[i] in ["-pscale", "--pscale"]):
                if (argv[i+1].lower())[0] == "y": 
                    g_params['isDrawScaleBar'] = True
                else:
                    g_params['isDrawScaleBar'] = False
                i = i + 2
            elif (argv[i] in ["-ptag", "--ptag"]):
                if (argv[i+1].lower())[0] == "y": 
                    g_params['isDrawTagColumn'] = True
                else:
                    g_params['isDrawTagColumn'] = False
                i = i + 2
            elif (argv[i] in ["-shrink", "--shrink"]):
                if (argv[i+1].lower())[0] == "y": 
                    g_params['isShrink'] = True
                else:
                    g_params['isShrink'] = False
                i = i + 2
            elif (argv[i] in ["-shrinkrate", "--shrinkrate"]):
                g_params['shrinkrate'], i = myfunc.my_getopt_float(argv, i)
            elif (argv[i] in ["-shrinkrateTM", "--shrinkrateTM"]):
                g_params['shrinkrate_TM'], i = myfunc.my_getopt_float(argv, i)
            elif (argv[i] in ["-imagescale", "--imagescale"]):
                g_params['image_scale'], i = myfunc.my_getopt_float(argv, i)
            elif (argv[i] in ["-h2wratio", "--h2wratio"]):
                g_params['H2W_ratio'], i = myfunc.my_getopt_float(argv, i)
            elif (argv[i] in ["-max-hold-loop", "--max-hold-loop"]):
                g_params['max_hold_loop'], i = myfunc.my_getopt_int(argv, i)
            elif (argv[i] in ["-m-shrink", "--m-shrink"]):
                g_params['method_shrink'], i = myfunc.my_getopt_int(argv, i)
            elif (argv[i] in ["-autosize", "--autosize"]):
                if (argv[i+1].lower())[0] == "y": 
                    g_params['isAutoSize'] = True
                else:
                    g_params['isAutoSize'] = False
                i = i + 2
            elif argv[i] in ["--aaseq", "-aaseq"]:
                aaSeqFile, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in["-colorTMbox", "--colorTMbox"]:
                g_params['isColorWholeTMbox'] = True; i += 1
            elif argv[i] in["-advtopo", "--advtopo"]:
                g_params['isAdvTopo'] = True; i += 1
            elif argv[i] in["-TMname", "--TMname"]:
                g_params['isTMname'] = True
                g_params['TMname']=argv[i+1].split(',')
                i += 2
            elif argv[i] in["-colorkingdom", "--colorkingdom"]:
                g_params['isColorByKingdom'] = True; i += 1
            elif argv[i] in["-showTMidx", "--showTMidx"]:
                g_params['isShowTMIndex'] = True; i += 1
            elif argv[i] in["-cleanplot", "--cleanplot"]:
                g_params['makeCleanPlot'] = True; i += 1
            elif argv[i] in["-showgap", "--showgap"]:
                g_params['isShowGap'] = True; i += 1
            elif argv[i] in["-debug", "--debug"]:
                g_params['isPrintDebugInfo'] = True; i += 1
            elif argv[i] == "-q":
                g_params['isQuiet'] = True; i +=  1
            else:
                print(("Error! Wrong argument:%s" % argv[i]), file=sys.stderr)
                return 1
        else:
            filelist.append(argv[i])
            i=i+1
#}}}

    if g_params['isDrawDGprofile']:
        g_params['dgscanProg'] = "%s/myscanDG.pl"%(rundir)

    myfunc.setup_logging(g_params['log_config_file'])
    logger = logging.getLogger(__name__)

    if g_params['outpath'] != "" and not os.path.exists(g_params['outpath']):
        os.system("mkdir -p %s" % g_params['outpath'])

    if g_params['method'].lower()[0:2] == "sv":
        g_params['method'] = 'svg'
        g_params['outFormat'] = 'svg'
    elif g_params['method'].lower()[0:2] == 'ma':
        g_params['method'] = 'mat'
        g_params['outFormat'] = 'pdf'
    elif g_params['method'].lower() == 'core-rainbow':
        g_params['method'] = 'core-rainbow'
        g_params['outFormat'] = 'pdf'
    elif g_params['method'].lower()[0:2] == 'py':
        g_params['method'] = 'pyx'
        g_params['outFormat'] = 'pdf'
    else:
        g_params['method'] = 'pil'

    logger.debug("font_dir = %s"%(g_params['font_dir']))

    if filelistfile != "":
        try:
            fp = open(filelistfile,"r")
            filelist += fp.read().split()
            fp.close()
        except IOError:
            print("file %s does not exist." %(filelistfile), file=sys.stderr)
    if len(filelist) < 1:
        print("Error! Input file not set.", file=sys.stderr)

    g_params['fntScaleBar'] = ImageFont.truetype(g_params['font_dir'] +
            g_params['font'], g_params['font_size_scalebar'])
    g_params['fntTMbox'] = ImageFont.truetype(g_params['font_dir'] +
        g_params['font'], g_params['font_size_TMbox'])
    g_params['fntDGprofileLegend'] = ImageFont.truetype(g_params['font_dir'] +
        g_params['font'], g_params['font_size_scalebar'])

    if aaSeqFile != "" and os.path.exists(aaSeqFile):
        (idList, aaSeqList) = myfunc.ReadFasta_without_annotation(aaSeqFile)
        dd = {}
        for i in range(len(idList)):
            dd[idList[i]] = aaSeqList[i].replace("-", "") # gapless aaseq
        g_params['aaSeqDict'] = dd
        if g_params['isDrawKRBias'] and len(aaSeqList) < 1:
            print("aaSeq must be set when krbias is enabled", file=sys.stderr)
            return 1


    logger.debug("method=%s"%(g_params['method']))
    for inFile in filelist:
        if g_params['method'] == 'pil':
            DrawMSATopo_PIL(inFile, g_params)
        elif g_params['method'] == 'svg':
            DrawMSATopo_SVG(inFile, g_params)
        elif g_params['method'] == 'mat':
            #DrawMSATopo_MAT(inFile, g_params)
            DrawMSATopo_MAT2(inFile, g_params)
        elif g_params['method'] == 'core-rainbow':
            DrawMSATopo_MAT_Core_unalign_rainbow(inFile, g_params)
        elif g_params['method'] == 'matfull':
            DrawMSATopo_MAT(inFile, g_params)
        elif g_params['method'] == 'pyx':
            DrawMSATopo_PYX(inFile, g_params)
    return 0
#}}}
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['outpath'] = ""
    g_params['isQuiet'] = False
    g_params['outFormat'] = "png"
    g_params['mode'] = "P"
    g_params['image_scale'] = None   #overal image scale, set it to a higher value to increase the DPI
    g_params['font_dir'] = "%s/../fonts/truetype/ttf-dejavu/"%(rundir)
    g_params['font_size'] = 16

    # font size for scale bar always set to 11 for easy reading
    g_params['font_size_scalebar'] = 11

    g_params['heightTMbox'] = 3; # number of lines for the TM helices box
    g_params['font_size_TMbox'] = 36; # size of the font for text written in TM
                                      # box
    g_params['font'] = "DejaVuSansMono.ttf"
    g_params['isAutoSize'] = True
    g_params['isDrawSeprationLine'] = True
    g_params['isDrawText'] = True
    g_params['isDrawKRBias'] = False

    g_params['DGProfileFile'] = ''

    #draw distribution of percentage of M
    g_params['isDrawPerMDistribution'] = True
    g_params['isDrawDGprofile'] = False
    g_params['isDrawTagColumn'] = False
    g_params['maxDistKR'] = 12
    g_params['isShrink'] = True
    g_params['aapath'] = ""
    g_params['MAXIMAGESIZE'] = 50*1024*1024; #50M in pixels
    g_params['marginX'] = 20; # marginX in pixels
    g_params['marginY'] = 50; # marginY in pixels
    g_params['H2W_ratio'] = None # set the height/width ratio for PIL image

    # number of columns for the annotation text
    g_params['widthAnnotation'] = 30
    # number of columns between the annotation and alignment
    g_params['annoSeqInterval'] = 4
    # number of lines for the scale bar
    g_params['heightScaleBar'] = 3
    g_params['scaleSeprationLine'] = 1
    g_params['GAP'] = "-"
    g_params['aaSeqDict'] = {}
    g_params['window_size'] = 70
    g_params['htmlheader'] = ""
    g_params['colorhtml'] = False
    g_params['method_shrink'] = 1
    g_params['isColorWholeTMbox'] = False
    g_params['isAdvTopo'] = False
    g_params['isShowGap'] = False
    g_params['isTMname'] = False
    g_params['TMname'] = []
    g_params['isColorByKingdom'] = False
    g_params['isShowTMIndex'] = False
    g_params['shrinkrate'] = None # shrink the alignment proportionally
    g_params['shrinkrate_TM'] = 2.0 # shrink the aligned TM region proportionally
    g_params['max_hold_loop'] = 12 # maximal positions to keep for the loop region when doing shrinking

    g_params['log_config_file'] = "%s/default_log.yml"%(rundir)
    g_params['logger'] = None
    g_params['dgscanProg'] = ""
    g_params['makeCleanPlot'] = False

    g_params['pdfcrop_margin_left']    = 20
    g_params['pdfcrop_margin_top']     = 5
    g_params['pdfcrop_margin_right']   = 10
    g_params['pdfcrop_margin_bottom']  = 5
    g_params['MAX_SIZE_ANNOTATION'] = 30 # maximum width of annotation


    g_params['memcolor_out_to_in'] = "#DCDCDC"  #very light grey, type = M
    g_params['memcolor_in_to_out'] = "#808080"  #grey           , type = W
    g_params['memcolor_out_to_in_MSA'] = "#FF6666"  # light red, type M
    g_params['memcolor_in_to_out_MSA'] = "#CC0000"  # dark red, type W
    #g_params['loopcolor_in'] = "#FFBFB3"        # light red
    g_params['loopcolor_in'] = "#FFFF00"        # yellow
    g_params['loopcolor_in_MSA'] = "#F2EABD"        # yellow
    #g_params['loopcolor_out'] = "#87CEFA"       # light sky blue
    g_params['loopcolor_out'] = "#3399FF"       # blue
    g_params['loopcolor_out_MSA'] = "#CCFFFF"       # faded blue
    g_params['spcolor'] = "#000000"       # signal peptide, black

    g_params['method'] = 'pil' #pyx
    g_params['isPrintDebugInfo'] = False
    g_params['isDrawMSA'] = True
    g_params['isDrawScaleBar'] = True
    g_params['isDrawDGProfileLegend'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    outpath_log = "logs"
    if not os.path.exists(outpath_log):
        os.makedirs(outpath_log)
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
    #cProfile.run("main()")
    # Check argv
