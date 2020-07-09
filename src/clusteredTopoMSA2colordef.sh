#!/bin/bash

# Created 2012-03-16, updated 2012-03-16, Nanjiang Shu 
# create color strips dataset for ITOL phylogenetic tree drawing from 
# clustered.orig.topomsa.fa

infile=$1

if [ "$infile" == "" ]; then 
    echo "usage:   clusteredTopoMSA2colordef.sh infile"
    exit
fi

red=#FF0000
green=#00FF00
blue=#0000FF
yellow=#FFFF00
pink=#FF00FF
cyan=#00FFFF
black=#000000

# write settings
dataset_settings="\
DATASET_COLORSTRIP
#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#

SEPARATOR TAB

#label is used in the legend table (can be changed later)
DATASET_LABEL\tClusters by the number of TMs

#dataset color (can be changed later)
COLOR\t#ff0000

#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#


#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#

#width of the colored strip
#STRIP_WIDTH 25

#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.
#MARGIN 0

#border width; if set above 0, a border of specified width (in pixels) will be drawn around the color strip 
#BORDER_WIDTH 0

#border color; used when BORDER_WIDTH is above 0
#BORDER_COLOR #0000ff

#always show internal values; if set, values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
#SHOW_INTERNAL 0


#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages

#=================================================================#
#       Actual data follows after the DATA keyword              #
#=================================================================#
DATA
"

echo -e "$dataset_settings"

awk -v red=$red -v green=$green -v blue=$blue -v yellow=$yellow -v pink=$pink -v cyan=$cyan -v black=$black '
/^>/ {
sub("^>*", "");
seqid=$1
gsub(",","",seqid)
tag=$3;
if (tag=="ClusterNo=1"){
    print seqid "\t" "#0000FF"
}
else if (tag=="ClusterNo=2"){
    print seqid "\t" "#6363FF"
}
else if (tag=="ClusterNo=3"){
    print seqid "\t" "#8080FF"
}
else if (tag=="ClusterNo=4"){
    print seqid "\t" "#ABABFF"
}
else if (tag=="ClusterNo=5"){
    print seqid "\t" "#9B84FF"
}
else if (tag=="ClusterNo=6"){
    print seqid "\t" "#A0DFFF"
}
else if (tag=="ClusterNo=7"){
    print seqid "\t" "#DEECFF"
}
else{
    print seqid "\t" "#EDFFFF"
}
}
END{
print ""
}
' $infile
