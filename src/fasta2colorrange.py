#!/usr/bin/env python

# Generate ITOL (3) color range definition file from FASTA annotation line,
# with taxonomy annotation
import os
import sys
import myfunc
import libtopologycmp as lcmp
from ete3 import Tree

fastafile = sys.argv[1]
treefile = sys.argv[2]
red="#ff0000"
green="#00ff00"
blue="#0000ff"
purple="#800080"
GAP = '-'

# read in taxonomy def
if not os.path.exists(fastafile):
    print >> sys.stderr, "Error! file fastafile (%s) does not exist." %fastafile;
    sys.exit(1)
if not os.path.exists(treefile):
    print >> sys.stderr, "Error! file treefile (%s) does not exist." %treefile;
    sys.exit(1)


t = Tree(treefile)
leaves = t.get_leaves()
leafNameList = [x.name for x in leaves]
leafNameSet = set(leafNameList)
(idList, annotationList, seqList) = myfunc.ReadFasta(fastafile)

# write out taxdef
fpout = sys.stdout
numSeq = len(idList)

# write settings
dataset_settings="""\
TREE_COLORS
#use this template to define branch colors and styles, colored ranges and label colors/font styles/backgrounds
#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file (except in the SEPARATOR line, which uses space).

SEPARATOR TAB

#First 3 fields define the node, type and color
#Possible types are:
#'range': defines a colored range (colored background for labels/clade)
#'clade': defines color/style for all branches in a clade
#'branch': defines color/style for a single branch
#'label': defines font color/style for the leaf label
#'label_background': defines the leaf label background color

#The following additional fields are required:
#for 'range', field 4 defines the colored range label (used in the legend)

#The following additional fields are optional:
#for 'label', field 4 defines the font style ('normal',''bold', 'italic' or 'bold-italic') and field 5 defines the numeric scale factor for the font size (eg. with value 2, font size for that label will be 2x the standard size)
#for 'clade' and 'branch', field 4 defines the branch style ('normal' or 'dashed') and field 5 defines the branch width scale factor (eg. with value 0.5, branch width for that clade will be 0.5 the standard width)

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA

"""

fpout.write(dataset_settings)
colorDict = {
        "Archaea": blue,
        "Bacteria": purple,
        "Eukaryota": green
        }
taxoList = ["Archaea", "Bacteria", "Eukaryota"]


for i in xrange(numSeq):
    seqID = idList[i]
    anno = annotationList[i]
    taxo = anno.split(',')[0].split('|')[-1]
    if seqID in leafNameSet and taxo in taxoList:
        color = colorDict[taxo]
        #fpout.write("%s\t%s\t%s\n"%(idList[i], "label_background", color ))
        fpout.write("%s\t%s\t%s\t%s\t%s\n"%(idList[i], "label", color , "bold", "2"))
fpout.write("\n")

if fpout != sys.stdout:
    fpout.close()
