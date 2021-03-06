#!/usr/bin/env python
# --Python 2.7+--

import os
import sys
import argparse
import hashlib
import myfunc
import Bio.SubsMat.MatrixInfo
import libtopologycmp as lcmp
import logging
import logging.config
import yaml

blosum62 = Bio.SubsMat.MatrixInfo.blosum62

GAP = '-'
alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

logger = logging.getLogger(__name__)

def CheckFileExt(choices):# {{{
    class Act(argparse.Action):
        def __call__(self,parser,namespace,fname,option_string=None):
            ext = os.path.splitext(fname)[1][1:]
            if ext not in choices:
                option_string = '({})'.format(option_string) if option_string else ''
                parser.error("file doesn't end with one of {}{}".format(choices,option_string))
            else:
                setattr(namespace,self.dest,fname)

    return Act
# }}}
def WriteHTMLHeader(htmlheader, fpout):# {{{
    header = """
<!DOCTYPE html>
<html>
<body>
<style>
.TM_Nin {
 border-style:solid;
 border-color:#287EC7;
 background-color: #808080;
}
.TM_Nout {
 border-style:solid;
 border-color:#287EC7;
 background-color: #DCDCDC;
}
</style>
"""
    if not g_params['makeCleanPlot']:
        header += "\n<h2>%s</h2>"%(htmlheader)
    fpout.write("%s\n"%(header))
# }}}
def WriteHTMLTail(fpout):# {{{
    tail = """
</body>
</html>
"""
    fpout.write("%s\n"%(tail))
# }}}

def WriteHTMLAlignment2(aln_name, idList, annoList, alignedTopoSeqList,#{{{
        originalAlignedTopoSeqList, aaSeqList, final2seq_idxMapList,
        fpout):
    logger = logging.getLogger(__name__)
# for two sequence pairwise alignment
# assign also the identical and similarity by using BLOSUM62
    #WIDTH = 90 # but do not break in the middle of a helix, adjust 1
    #WIDTH = 60 # but do not break in the middle of a helix
    WIDTH = g_params['window_size']
    base_outline_width = 1
    base_text_color = 'black'
    base_outline_color = 'black'

    maxSizeAnno = max([len(x) for x in idList])
    if g_params['makeCleanPlot']:
        maxSizeAnno = 2
    lengthAlignment = len(alignedTopoSeqList[0])
    numSeq = len(idList)
    posTMList = []
    typeTMList = []
    TMnameList = []
    foldTypeList = []
    for i in range(numSeq):
        (posTM,typeTM) = lcmp.GetTMType(alignedTopoSeqList[i])
        TMname = myfunc.GetTMnameFromAnnotation(annoList[i])
        foldType = myfunc.GetFoldTypeFromAnnotation(annoList[i])
        posTMList.append(posTM)
        typeTMList.append(typeTM)
        TMnameList.append(TMname)
        foldTypeList.append(foldType)
    #print(annoList)

    if g_params['makeCleanPlot']:
        g_params['colorhtml'] = False

    if g_params['colorhtml']:
        color_TM = 'red'
        color_nonTM = 'grey'
    else:
        color_TM = 'black'
        color_nonTM = '#646464'


    fpout.write("<p>\n")
    if not g_params['makeCleanPlot']:
        fpout.write("<h4>Alignment for %s</h4>\n"%(aln_name))
    fpout.write("<pre>\n")
    j = 0 # iterator for the alignment position
    isStart = True
    cnt = 0
    jL = [0]*numSeq # iterator for each sequence
    while j < lengthAlignment:
        #print("j=%d, lengthAlignment=%d"%(j, lengthAlignment))
        if isStart: # at beginning of each line
        # writing annotation and start index
            strs = [""]*(numSeq+1)
            for i in range(numSeq):

                seqlabel = annoList[i]
                if g_params['makeCleanPlot']:
                    seqlabel = alphabet[i]
                try:
                    strs[i] += "%-*s %4d "%(maxSizeAnno, seqlabel,
                            final2seq_idxMapList[i][j]+1)
                except KeyError:
                    logger.debug("final2seq_idxMapList error  i=%d, j=%d"%(i,j))
                    pass
            strs[2] += "%-*s %4s "%(maxSizeAnno, "", "")
            isStart = False


        # Writing alignment
        isWithinTMregion = False
        for i in range(numSeq):
            posTM = posTMList[i]
            typeTM = typeTMList[i]
            TMname = TMnameList[i]

            idxTM = lcmp.GetTMIndex(j, posTMList[i])
            if idxTM >= 0:
                (b, e) = posTMList[i][idxTM]
                (text_TMlabel, text_color, outline_color, outline_width) = lcmp.SetMakeTMplotColor(
                        idxTM, TMnameList[i],  foldTypeList[i],
                        base_outline_width, base_text_color,
                        base_outline_color)
                aa_segment = aaSeqList[i][b:e].upper()
                if typeTM[idxTM] == 'M':
                    bgcolor = g_params['memcolor_out_to_in']
                elif typeTM[idxTM] == 'W':
                    bgcolor = g_params['memcolor_in_to_out']
                elif typeTM[idxTM] == 'R':
                    bgcolor = g_params['loopcolor_in']
                elif typeTM[idxTM] == 'r':
                    bgcolor = g_params['loopcolor_out']
                isWithinTMregion = True # if hit TM region of any sequence, set as TRUE

                if g_params['isBreakTM']:
                    jL[i] += 1
                    aa = aaSeqList[i][j].upper()
                    strs[i] += "<b><font style=\"background-color:%s; border-style:solid none solid none; border-color:%s; \" color=\"%s\">%s</font></b>"%(
                            bgcolor, outline_color, color_TM, aa)
                else:
                    if j >= jL[i]:
                        #print ("outline_color=", outline_color)
                        strs[i] += "<b><font style=\"background-color:%s; border-style:solid none solid none; border-color:%s; \" color=\"%s\">%s</font></b>"%(
                                bgcolor, outline_color, color_TM, aa_segment)
                        jL[i] = e
            else: #loop
                #posLoop = lcmp.GetLoopBeginEnd(j, posTMList[i], lengthAlignment)
                #(b, e) = posLoop
                #aa_segment = aaSeqList[i][b:e].lower()
                #jL[i] = e
                aa = aaSeqList[i][j].lower()
                top_aa = alignedTopoSeqList[i][j].lower()
                if top_aa == 'i':
                    loopcolor = g_params['loopcolor_in']
                elif top_aa == 'o':
                    loopcolor = g_params['loopcolor_out']
                else:
                    loopcolor = "none"
                #print(loopcolor)
                strs[i] += "<font style=\"background-color:%s\" color=\"%s\">%s</font>"%(loopcolor, color_nonTM, aa)
                jL[i] += 1

        #Writing alignment relation
        while j < min(jL) and j < lengthAlignment:
            aa1 = aaSeqList[0][j].upper()
            aa2 = aaSeqList[1][j].upper()
            strs[2] += myfunc.GetAlignmentRelationship(aa1, aa2, GAP, blosum62)
            j += 1
            cnt += 1
            #print("++ i=%d, j=%d"%(i,j))

        if ((cnt >= WIDTH and (g_params['isBreakTM'] or isWithinTMregion == False)) 
                or (j >= lengthAlignment)
                ):
            for i in range(numSeq):
                #print("i=%d, j=%d"%(i,j))
                strs[i] += " %4d"%(final2seq_idxMapList[i][j-1]+1)

            fpout.write("%s\n"%(strs[0]))
            if g_params['showRelationship']:
                fpout.write("%s\n"%(strs[2])) #relationship
            fpout.write("%s\n"%(strs[1]))
            fpout.write("\n\n")

            # init
            strs = [""]*(numSeq+1)
            isStart = True
            cnt = 0

    fpout.write("</pre>\n")
    fpout.write("</p>\n")

#}}}

def WriteSeqAlnHTML(seqAlnFileList, extTopoMSA, outfile):# {{{
    try:
        fpout = open(outfile,"w")
    except IOError:
        print("Failed to write to %s"%(outfile), file=sys.stderr)
        return 1
    WriteHTMLHeader('Alignment highlighted by <font color=%s>TM regions</font>'%('red'), fpout)
    print("Processed alignments:")
    for alnfile in seqAlnFileList:
        rootname_alnfile = os.path.basename(os.path.splitext(alnfile)[0])
        topomsafile = '.'.join([os.path.splitext(alnfile)[0], extTopoMSA])
        if not (os.path.exists(alnfile) and os.path.exists(topomsafile)):
            if not os.path.exists(alnfile):
                sys.stderr.write('alnfile %s does not exist\n'%(alnfile))
            if not os.path.exists(topomsafile):
                sys.stderr.write('topomsafile %s does not exist\n'%(topomsafile))
            continue
        (seqIDList, seqAnnoList, seqList) = myfunc.ReadFasta(alnfile)
        #print(seqIDList)
        (topoIDList, topoAnnoList, topoList) = myfunc.ReadFasta(topomsafile)
        #print(topoIDList)
        if g_params['removeUnnecessaryGap']:
            seqList = lcmp.RemoveUnnecessaryGap(seqList)
            topoList = lcmp.RemoveUnnecessaryGap(topoList)

        # since there is no shrinking, index map is always p->p
        final2seq_idxMapList = []
        for i in range(len(seqIDList)):
            seqlength = len(seqList[i])
            idxmap = {}
            for j in range(seqlength):
                idxmap[j] = j
            final2seq_idxMapList.append(idxmap)

        print(('\t'+rootname_alnfile))
        WriteHTMLAlignment2(rootname_alnfile, topoIDList, topoAnnoList, topoList,
                topoList, seqList, final2seq_idxMapList,
                fpout)


    WriteHTMLTail(fpout)

    fpout.close()
    return 0
# }}}

def main(g_params):#{{{
#metavar='' is the text shown after then option argument
    parser = argparse.ArgumentParser(
            description='Write sequence alignment highlighted by TM region',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='''\
Created 2019-12-30, updated 2019-12-30, Nanjiang Shu

Examples:

    %s fam1.aln fam2.aln

# do not show alignment relationship

    %s fam1.aln fam2.aln -norel

'''%(sys.argv[0], sys.argv[0]))
    parser.add_argument('seqAlnFileList', metavar='seqAlnFileList', type=str,
            nargs='+', help='Provide sequence alignment files in FASTA format')
    parser.add_argument('-ext-topomsa', dest='extTopoMSA',
            metavar='extTopoMSA', type=str, default='topomsa',
            help='Set the file extension for topology msa')
    parser.add_argument('-o', metavar='OUTFILE', dest='outfile', 
            action=CheckFileExt({'html'}), required=True,
            help='Output the result to THML file')
    parser.add_argument('-ws', dest='window_size',
            metavar='window_size', type=int, default=60,
            help='Set the file extension for topology msa')
    parser.add_argument('-norel', action='store_true', 
            help='Do not show alignment relationship')
    parser.add_argument('-rmgap', action='store_true', 
            help='Remove Unnecessary gap')
    parser.add_argument('-breakTM', action='store_true', 
            help='Break the TM helices to make equal length of each aligned line')
    parser.add_argument('-cleanplot', action='store_true', 
            help='Make clean plot, no header text, sequences labeled as A, B')

    args = parser.parse_args()

    seqAlnFileList = args.seqAlnFileList
    extTopoMSA = args.extTopoMSA
    g_params['window_size'] = args.window_size
    outfile = args.outfile
    if args.norel:
        g_params['showRelationship'] = False
    if args.rmgap:
        g_params['removeUnnecessaryGap'] = True
    if args.cleanplot:
        g_params['makeCleanPlot'] = True
    if args.breakTM:
        g_params['isBreakTM'] = True

    lcmp.SetMakeTMplotColor_g_params(g_params)
    WriteSeqAlnHTML(seqAlnFileList, extTopoMSA, outfile)


#}}}
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['window_size'] = 60
    g_params['colorhtml'] = True
    g_params['showRelationship'] = True
    g_params['removeUnnecessaryGap'] = False
    g_params['makeCleanPlot'] = False
    g_params['isBreakTM'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
