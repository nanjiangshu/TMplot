#!/bin/bash

# Filename: merge_tmplot.sh
# Description: Merge TM plots
progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname [-l LISTFILE] [FILE [FILE...]]  [-o OUTFILE]
Options:
  -l       FILE    Set the fileListFile, one filename per line
  -cap     STR     Set captions for the figure, a new cap flag will be written from a new line.
  -o       FILE    Output the PDF file to FILE
  -q               Quiet mode
  -h, --help       Print this help message and exit

Created 2020-06-26, updated 2020-06-26, Nanjiang Shu 

Examples:
    $progname file1 file2 -cap caption
"
PrintHelp(){ #{{{
    echo "$usage"
}
#}}}
AddAbsolutePath(){ #$fileordir#{{{
    #convert a path to absolute path, 
    #return "" if the file or path does not exist
    if [ ! -e "$1" ]; then 
        return 1
    fi
    curr="$PWD"
    if [ -d "$1" ]; then
        cd "$1" && echo "$PWD"       
    else 
        file="$(basename "$1")"
        cd "$(dirname "$1")" && dir="$PWD"
        echo "$dir/$file"
    fi
    cd "$curr"
    return 0
}
#}}}
IsProgExist(){ #{{{
    # usage: IsProgExist prog
    # prog can be both with or without absolute path
    type -P $1 &>/dev/null \
        || { echo The program \'$1\' is required but not installed. \
        Aborting $0 >&2; exit 1; }
    return 0
}
#}}}
IsPathExist(){ #{{{
# supply the effective path of the program 
    if ! test -d "$1"; then
        echo Directory \'$1\' does not exist. Aborting $0 >&2
        exit 1
    fi
}
#}}}
WriteHead(){ #{{{
echo '\documentclass[12pt]{article}
\usepackage{booktabs}
\usepackage{amsmath}
\usepackage{hyperref}
\usepackage{color}
\usepackage{enumerate}
\pagestyle{plain}
\usepackage[small,bf]{caption}
\usepackage{subfig}
\usepackage{graphicx}
\usepackage{grffile}
\usepackage[margin=0.2in]{geometry}
\usepackage{float}
\captionsetup{justification=raggedright,singlelinecheck=false}
'
} 
#}}}
WriteTail(){ #{{{
echo '\end{document}'
} 
#}}}
WriteTexFile() { #{{{
    local texfile="$1"
    local file=
    local absfile=
    local modfilename=
    local NROW=
    local pd=
    local cnt=0
    local title=
    local caption=

    for ((i=0;i<numCaption;i++));do
        ss=${captionList[$i]}
        ss=$(echo $ss|sed 's/_/\\_/g')
        if ((i==0)); then
            caption=${ss}
        else
            caption="$caption\\\\${ss}"
        fi
    done

    (
    WriteHead
        echo "
\\begin{document}
\\begin{table} [H]
\\caption*{$caption}
"
    for ((i=0;i<numFile;i++));do
        file=${fileList[$i]}
        absfile=`AddAbsolutePath $file`
        modfilename=`echo $file | sed 's/_/\\\\_/g'`
        echo "\\subfloat{\\includegraphics[width=$width\\textwidth]{$absfile}}
\\
        "
    done 
    echo "\\end{table}"
    WriteTail
    ) > $texfile
}
#}}}
GetWidth() { #{{{
    n=$1
    case $n in 
        1) echo "0.95";;
        2) echo "0.49";;
        3) echo "0.33";;
        4) echo "0.24";;
        5) echo "0.19";;
        6) echo "0.12";;
    esac
}
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

isQuiet=0
outfile=
fileListFile=
fileList=()
titleList=()
captionList=()
numFigPerRow=1

isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        fileList+=("$1")
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -o|--o) outfile=$2;shift;;
            -cap*|--cap*) captionList+=("$2");shift;;
            -l|--l|-list|--list) fileListFile=$2;shift;;
            -q|-quiet|--quiet) isQuiet=1;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        fileList+=("$1")
    fi
    shift
done

if [ "$fileListFile" != ""  ]; then 
    if [ -s "$fileListFile" ]; then 
        while read line
        do
            fileList+=("$line")
        done < $fileListFile
    else
        echo listfile \'$fileListFile\' does not exist or empty. >&2
    fi
fi

numTitle=${#titleList[@]}
numCaption=${#captionList[@]}
numFile=${#fileList[@]}
if [ $numFile -eq 0  ]; then
    echo Input not set! Exit. >&2
    exit 1
fi

basedir=/tmp/viewpdfs_latex
if [ ! -d $basedir ]; then
    mkdir -p $basedir
fi
tmpdir=$(mktemp -d $basedir/tmpdir.merge_tmplot.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }  
texfile=$tmpdir/view.tex
width=`GetWidth $numFigPerRow`

WriteTexFile "$texfile"

osname=`uname -s`

currdir=$PWD
cd $tmpdir
pdflatex $texfile > $tmpdir/pdflatex.err 2>&1
cd $currdir

if [ "$outfile" != "" ];then
    /bin/cp -f $tmpdir/view.pdf $outfile
    /bin/cp -f $tmpdir/view.tex $outfile.tex
    echo "$outfile output"
    rm -rf $tmpdir
    pdffile=$outfile
else
    echo "pdffile can be found at $tmpdir/$pdffile"
    pdffile=$tmpdir/view.pdf
fi

case $osname in 
    Linux*) myopen $pdffile;;
    Darwin*) open $pdffile;;
esac
