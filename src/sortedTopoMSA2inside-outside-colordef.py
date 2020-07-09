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
numSeq = len(idList)
for i in range(numSeq):
    gid = idList[i]
    if gid != 'Consensus':
        color = red
        NtermState = lcmp.GetNtermState(seqList[i])
        if NtermState == 'o':
            color = blue
        sys.stdout.write("%s,%s\n"%(gid, color ))
sys.stdout.write("\n")
