#!/usr/bin/env python
# --Python 2.7+--

import os
import sys
import argparse
import hashlib
import myfunc

def WriteHTMLHeader(htmlheader, fpout):# {{{
    header = """
<!DOCTYPE html>
<html>
<body>
<h3>%s</h3>
<pre>
"""%(htmlheader)
    fpout.write("%s\n"%(header))
# }}}
def WriteHTMLTail(htmltail, fpout):# {{{
    tail = """
</pre>
</body>
</html>
"""
    fpout.write("%s\n"%(tail))
# }}}

def WriteHTMLAlignment2(idList, annoList, alignedTopoSeqList,#{{{
        originalAlignedTopoSeqList, aaSeqList, final2seq_idxMapList,
        fpout):
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


    strs = [""]*numSeq
    j = 0 # iterator for the alignment position
    isStart = True
    cnt = 0
    while j < lengthAlignment:
        if isStart:
            strs = [""]*(numSeq+1)
            for i in xrange(numSeq):
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

        for i in xrange(numSeq):
            if IsWithinTMRegion(j, posTMList[i]):
                aa = aaSeqList[i][j].upper()
                isWithinTMregion = True # if hit TM region of any sequence, set as TRUE
                strs[i] += "<b><font color=\"%s\">%s</font></b>"%(color_TM, aa)
            else:
                aa = aaSeqList[i][j].lower()
                strs[i] += "<font color=\"%s\">%s</font>"%(color_nonTM, aa)
        if ((cnt >= WIDTH and isWithinTMregion == False) 
                or (j >= lengthAlignment-1)
                or j == 190):
            for i in xrange(numSeq):
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

#}}}

def WriteSeqAlnHTML(seqAlnFileList, extTopoMSA, outfile):
    try:
        fpout = open(outfile,"w")
    except IOError:
        print >> sys.stderr, "Failed to write to %s"%(outfile)
        return 1

    for alnfile in seqAlnFileList:


    fpout.close()
    return 0


def main(g_params):#{{{
#metavar='' is the text shown after then option argument
    parser = argparse.ArgumentParser(
            description='Write sequence alignment highlighted by TM region',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='''\
Created 2019-12-30, updated 2019-12-30, Nanjiang Shu

Examples:
%s seqs.aln
'''%(sys.argv[0]))
    parser.add_argument('seq-aln-files', metavar='seqAlnFileList', type=str, nargs='+',
                    help='Provide sequence alignment files')
    parser.add_argument('-ext-topomsa', metavar='extTopoMSA', type=str,
            default='topomsa', help='Set the file extension for topology msa')
    parser.add_argument('-o, --outfile', metavar='OUTFILE', dest='outfile',
            help='Output the result to file')

    args = parser.parse_args()

    seqAlnFileList = args.seqAlnFileList
    extTopoMSA = args.extTopoMSA
    outfile = args.outfile

    WriteSeqAlnHTML(seqAlnFileList, extTopoMSA, outfile)


#}}}
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
