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
import tempfile
import shutil
import platform
import distro
alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

if sys.platform == "linux" or sys.platform == "linux2":
    python_exec = "python"
elif sys.platform == "darwin":
    python_exec = "pythonw"

os_dist = distro.linux_distribution()[0]

rundir = os.path.dirname(os.path.realpath(__file__))

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
def CheckPrerequisite():# {{{
    """Check necessary installations
    """
    # TO BE IMPLEMENTED
    return True
# }}}
def MakeTMplot(seqAlnFile, topAlnFile, outpath, tmpdir):# {{{
    """Make topology plot for TM family.
    """
    rootname = os.path.basename(os.path.splitext(seqAlnFile)[0])
    basename_seqAlnFile = os.path.basename(seqAlnFile)
    basename_topAlnFile = os.path.basename(topAlnFile)
    ext_topAlnFile = os.path.splitext(topAlnFile)[1].lstrip('.')

    shutil.copy2(seqAlnFile, os.path.join(tmpdir, basename_seqAlnFile))
    shutil.copy2(topAlnFile, os.path.join(tmpdir, basename_topAlnFile))
    cwd = os.getcwd()

    os.chdir(tmpdir)
    # generate topology one line plot
    cmd = [python_exec, os.path.join(rundir, "drawMSATopo.py"), "-m-shrink",
        str(0), "-method", "pil",  "-pfm", "no", "-text", "n",  "-pdg", "n",
        "-pfm", "n",  "-pmsa", "y", "-ptag", "y", "-showTMidx", "-sep", "n",
        "--advtopo",   "-cleanplot", "-h2wratio", str(g_params["H2W_ratio"]),
        "-shrink", "no", "-showgap", basename_topAlnFile]

    if g_params['verbose']:
        print(("Generating toplogy alignment figure for %s"%(rootname)))

    (isCmdSuccess, t_runtime, t_msg) = myfunc.RunCmd(cmd)
    if not isCmdSuccess:
        print(t_msg)
        return 1
    topalnfigure = "%s.png"%(rootname)
    if not os.path.exists(topalnfigure):
        return 1
    # resize the figure file
    resized_topalnfigure = "%s.s%d.png"%(rootname, g_params['figure_resize'])
    shutil.copy2(topalnfigure, resized_topalnfigure)
    cmd = ["mogrify", "-resize", str(g_params['figure_resize']), resized_topalnfigure]
    if g_params['verbose']:
        print(("Resizing the topology alignment figure for %s"%(rootname)))
    (isCmdSuccess, t_runtime, t_msg) = myfunc.RunCmd(cmd)
    if not isCmdSuccess:
        print(t_msg)
        return 1

    # generate seqaln figure
    seqaln_htmlfigure = "%s.%s"%(rootname, "seqaln.html")
    cmd = [python_exec, os.path.join(rundir, "write_seqaln_colorTM.py"),
            basename_seqAlnFile, "-ext-topomsa", ext_topAlnFile, "-ws",
            str(g_params['window_size']), "-o",
            seqaln_htmlfigure, "-cleanplot", "-rmgap"]
    if g_params['isBreakTM']:
        cmd += ["-breakTM"]

    if g_params['verbose']:
        print(("Generating sequence alignment highlighted by TM regions for %s"%(rootname)))
    (isCmdSuccess, t_runtime, t_msg) = myfunc.RunCmd(cmd)
    if not isCmdSuccess:
        print(t_msg)
        return 1

    # convert html to pdf
    seqaln_pdffigure = "%s.%s"%(rootname, "seqaln.pdf")
    cmd = ["wkhtmltopdf",  seqaln_htmlfigure, seqaln_pdffigure]
    if os_dist.lower() in ["debian", "ubuntu"]:
        cmd = ["xvfb-run"] + cmd
    if g_params['verbose']:
        print("Convert the html figure to PDF for sequence alignment")
    (isCmdSuccess, t_runtime, t_msg) = myfunc.RunCmd(cmd)
    if not isCmdSuccess:
        print(t_msg)
        return 1
    # crop the PDF figure
    cmd = ["pdfcrop", seqaln_pdffigure]
    (isCmdSuccess, t_runtime, t_msg) = myfunc.RunCmd(cmd)
    if not isCmdSuccess:
        print(t_msg)
        return 1

    seqaln_pdffigure_crop =  "%s.%s"%(rootname, "seqaln-crop.pdf")

    # merge figures
    (seqIDList, seqAnnoList, seqList) = myfunc.ReadFasta(basename_topAlnFile)
    outfile = "%s.seqtopaln.pdf"%(rootname)
    cmd = ["bash", os.path.join(rundir, "merge_tmplot.sh"),
            resized_topalnfigure, seqaln_pdffigure_crop, "-cap", rootname, "-o", outfile]
    capList = []
    for i in range(len(seqIDList)):
        capList += ["-cap", "%s: %s"%(alphabet[i], seqIDList[i])]
    cmd += capList
    if g_params['verbose']:
        print(("Merging the topology alignment figure and sequence alignment figure for %s"%(rootname)))
    (isCmdSuccess, t_runtime, t_msg) = myfunc.RunCmd(cmd)
    if not isCmdSuccess:
        print(t_msg)
        return 1

    # copy the pdf figure generated by latex to a tmp file (a hack for the
    # PDFcrop
    tmpoutfile = "tt1.pdf"
    shutil.copy2(outfile, tmpoutfile)

    # crop the merged PDF figure
    cmd = ["pdfcrop", tmpoutfile]
    (isCmdSuccess, t_runtime, t_msg) = myfunc.RunCmd(cmd)
    if not isCmdSuccess:
        print(t_msg)
        return 1

    outfile_crop =  "tt1-crop.pdf"

    if os.path.exists(outfile_crop):
        final_targetfile =  os.path.join(outpath, "%s.seqtopaln.pdf"%(rootname))
        shutil.copy2(outfile_crop, final_targetfile)

    if g_params['verbose']:
        print(("Copy the result to final target %s"%(os.path.join(outpath, outfile))))


    os.chdir(cwd)

    return 0
# }}}


def main(g_params):#{{{
#metavar='' is the text shown after then option argument
    parser = argparse.ArgumentParser(
            description='Make plot for TM proteins',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='''\
Created 2020-06-26, updated 2020-06-26, Nanjiang Shu

Examples:

    %s -seqaln fam1.aln -topaln fam1.topomfa -outpath outdir

'''%(sys.argv[0]))
    parser.add_argument('-seqaln', dest='seqAlnFile',
            metavar='seqalnFile', type=str, required=True,
            help='Provide sequence alignment file')
    parser.add_argument('-topaln', dest='topAlnFile',
            metavar='topalnFile', type=str, required=True,
            help='Provide topology alignment file')
    parser.add_argument('-outpath', '--outpath', metavar='OUTPATH', dest='outpath',
            required=True,
            help='Output the result to outpath')
    parser.add_argument('-breakTM', action='store_true', 
            help='Break the TM helices to make equal length of each aligned line')
    parser.add_argument('-h2wratio', dest='H2W_ratio',
            metavar='H2W_ratio', type=float, default=0.08,
            help='Set the H2W ratio')

    args = parser.parse_args()

    seqAlnFile = args.seqAlnFile
    topAlnFile = args.topAlnFile
    g_params['H2W_ratio'] = args.H2W_ratio
    outpath = args.outpath
    outpath = os.path.abspath(outpath)
    if args.breakTM:
        g_params['isBreakTM'] = True

    if not CheckPrerequisite():
        return 1

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    tmpdir = tempfile.mkdtemp()
    if MakeTMplot(seqAlnFile, topAlnFile, outpath, tmpdir) == 0:
        shutil.rmtree(tmpdir)
    else:
        print(("makeplot failed for (%s, %s) "%(seqAlnFile, topAlnFile)))
        print(("temporary results can be found at %s"%(tmpdir)))

    return 0

#}}}
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['figure_resize'] = 5000
    g_params['window_size'] = 100
    g_params['isBreakTM'] = False
    g_params['verbose'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
