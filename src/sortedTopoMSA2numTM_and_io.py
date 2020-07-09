#!/usr/bin/env python

import sys
import os
import myfunc
import libtopologycmp as lcmp

infile = sys.argv[1]
red="#FF0000"
green="#00FF00"
blue="#0000FF"
GAP = '-'

# read in taxonomy def
if not os.path.exists(infile):
    print("Error! file infile (%s) does not exist." %infile, file=sys.stderr);
    sys.exit(1)


(idList, seqList) = myfunc.ReadFasta_without_annotation(infile)

# write out taxdef
fpout = sys.stdout
numSeq = len(idList)

# write settings
dataset_settings="""\
DATASET_MULTIBAR

SEPARATOR TAB


#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#label is used in the legend table (can be changed later)
DATASET_LABEL\tTree of TM evolution

#dataset color (can be changed later)
COLOR\t#ff0000

#define colors for each individual field column (use hexadecimal, RGB or RGBA notation; if using RGB/RGBA, COMMA cannot be used as SEPARATOR)
FIELD_COLORS\t#ff0000\t#0000ff

#field labels
FIELD_LABELS\tNterm-inside\tNterm-outside


#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#

#maximum width
WIDTH\t1000

#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.
#MARGIN\t0

#always show internal values; if set, values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
#SHOW_INTERNAL\t0

#bar height factor; Default bar height will be slightly less than the available space between leaves, but you can set a multiplication factor here to increase/decrease it (values to 1 will decrease it, values above 1 will increase it)
#HEIGHT_FACTOR\t1

#Bars are aligned to the node lines by default. Using BAR_SHIFT, you can move them all up/down by a fixed amount
#BAR_SHIFT\t0

#align individual fields; if set to 1, individual bar charts will not be stacked
#ALIGN_FIELDS\t0

#border width; if set above 0, a border of specified width (in pixels) will be drawn around the bars
#BORDER_WIDTH\t0

#border color; used when BORDER_WIDTH is above 0
#BORDER_COLOR\t#0000ff

#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA

"""

fpout.write(dataset_settings)


for i in range(numSeq):
    gid = idList[i]
    if gid != 'Consensus':
        n_i = 0
        n_o = 0
        NtermState = lcmp.GetNtermState(seqList[i])
        numTM = myfunc.CountTM(seqList[i])
        if NtermState == 'o':
            n_i = 0
            n_o = numTM
        else:
            n_i = numTM
            n_o = 0
        fpout.write("%s\t%d\t%d\n"%(gid, n_i,n_o ))
fpout.write("\n")

if fpout != sys.stdout:
    fpout.close()
