#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import click
import libtopologycmp as lcmp
import myfunc

@click.command()
@click.argument('alnfile')
@click.option('--outfile', default="", help='write the output to file')
@click.option('--method', default=1, help='0 for slower method, 1 for faster method')

def action(method, alnfile, outfile):
    (seqidList, seqAnnoList, seqList) = myfunc.ReadFasta(alnfile)
    if (method == 0):
        newSeqList = lcmp.RemoveUnnecessaryGap_old(seqList)
    else:
        newSeqList = lcmp.RemoveUnnecessaryGap(seqList)
    try:
        if outfile == "":
            fpout = sys.stdout
        else:
            fpout = open(outfile,"w")
        for i in range(len(seqidList)):
            fpout.write(">%s\n"%(seqAnnoList[i]))
            fpout.write("%s\n"%(newSeqList[i]))
        if fpout and fpout != sys.stdout:
            fpout.close()
        return 0
    except IOError:
        click.echo("Failed to write to file %s"%(outfile))
        return 1


if __name__ == '__main__':
    action()

