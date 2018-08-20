#!/usr/bin/env python

import os
import sys
from ete2 import Tree
import math
import myfunc
from colour import Color
blue = Color("blue")
red = Color("red")

# itol_path = os.environ['HOME'] + os.sep + ".local/lib/python2.7/itol"
# #itol_path = "/data3/downloads/itol/"
# if not os.path.exists(itol_path):
#     itol_path = "/afs/pdc.kth.se/home/n/nanjiang" + os.sep +  "Public/lib/python2.7/itol" 
# sys.path.append(itol_path)

from itolapi import Itol
from itolapi import ItolExport

usage="""
USAGE:  itol_pfamtree.py [-datapath DIR] -l pfamidlist [ID [ID ...]]
    Visualize phylogenetic tree of Pfam family, highlighting important features
    of membrane proteins
OPTIONS:
  -m, -method STR Method for visualization, (default: 0)
                  0:
                  1:
                  sd1: phylogenetic tree showing the number of the TM helices and with orientation
                       represented by either red or blue, and also color subfamilies differently
                       in the branches.
                  sd2: phylogenetic tree with domain architectures
                  sd3: phylogenetic tree with the first level (kingdom) colored in branches.
                       and the next three levels are colored in outer circles,(three circles) 

  -datapath DIR   Set datapath, (default: ./)
  -outpath  DIR   Set outpath, (default: ./)
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-03-13, updated 2017-08-13, Nanjiang Shu 
"""

def PrintHelp():#{{{
    print usage
    #}}}
def GetFontSize(numLeave):#{{{
    fontsize = 500/math.sqrt(numLeave)
    fontsize = max(30, fontsize)
    fontsize = min(200, fontsize)

    return fontsize
#}}}
def WriteDomainColorDefFile(domain_colordef_file, domain_dict, domain_idlist, seqlen_dict, leaves_name_set):#{{{
    """Write domain color definition file for iTOL given domain file
    """
    default_color = "#008000"
    default_shape = "HH"
    lst_color = list(blue.range_to(red,len(domain_idlist)))
    color_dict = {}
    for i in xrange(len(domain_idlist)):
        domainid = domain_idlist[i]
        color = lst_color[i].get_hex_l()
        color_dict[domainid] = lst_color[i].get_hex_l()

    try:
        fpout = open(domain_colordef_file, "w")
        # write the domain colordef file
        for seqid in domain_dict:
            if not seqid in leaves_name_set:
                continue
            seqlen = seqlen_dict[seqid]
            fpout.write("%s,%d"%(seqid, seqlen))
            shape = default_shape
            dlist = domain_dict[seqid]
            for dm in dlist:
                domainid = dm[2]
                color = color_dict[domainid]
                #color = default_color
                fpout.write(",%s|%d|%d|%s|%s"%(shape, dm[0]+1, dm[1]+1, color, domainid))
            fpout.write("\n")
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%(domain_colordef_file)
        return 1

#}}}
def Itol_Tree_m0(pfamid, datapath, outpath):#{{{
#Create the Itol class
    itl = Itol.Itol()
#Set the tree file
    tree = datapath + os.sep + pfamid + '.kalignp.tree'
    (dataset1, dataset2, dataset3, dataset4) = ("", "", "", "")
    if not os.path.exists(tree):
        print >> sys.stderr, "tree file %s does not exist. Ignore" %(tree)
        return 1
    t = Tree(tree)
    leaves = t.get_leaves()
    numLeave = len(leaves)

    fontsize = GetFontSize(numLeave)

    dataset1 = datapath + os.sep + pfamid + '.numTM_and_io.txt'
    dataset2 = datapath + os.sep + pfamid + '.cmpclass.colordef.txt'
#    dataset3 = datapath + os.sep + pfamid + '.ntermstate.colordef.txt'
    dataset4 = datapath + os.sep + pfamid + '.cluster.colordef.txt'

    colordeffile = datapath + os.sep + pfamid + '.pfam.colordef.txt'
    branchlabelfile = datapath + os.sep + pfamid + '.branchlabel.txt'

#===================================
    itl.add_variable('treeFile',tree)
    itl.add_variable('treeName','PF00854')
    itl.add_variable('treeFormat','newick')
    if os.path.exists(colordeffile):
        itl.add_variable('colorDefinitionFile', colordeffile)
    if os.path.exists(branchlabelfile):
        itl.add_variable('branchLabelsFile', branchlabelfile)

#===================================
    if os.path.exists(dataset1):
        itl.add_variable('dataset1File',dataset1)
        itl.add_variable('dataset1Label','numTM_and_io')
        itl.add_variable('dataset1Separator','comma')
        itl.add_variable('dataset1Type','multibar')
        itl.add_variable('dataset1PreventOverlap','1')
        itl.add_variable('dataset1Color','#FF0000')

#===================================
    if os.path.exists(dataset2):
        itl.add_variable('dataset2File', dataset2)
        itl.add_variable('dataset2Label', 'cmpclass')
        itl.add_variable('dataset2Separator','comma')
        itl.add_variable('dataset2Type','colorstrip')
        itl.add_variable('dataset2StripWidth','200')
        itl.add_variable('dataset2PreventOverlap','1')
        itl.add_variable('dataset2ColoringType','both')

#===================================
    if os.path.exists(dataset3):
        itl.add_variable('dataset3File', dataset3)
        itl.add_variable('dataset3Label', 'NtermState')
        itl.add_variable('dataset3Separator','comma')
        itl.add_variable('dataset3Type','colorstrip')
        itl.add_variable('dataset3StripWidth','200')
        itl.add_variable('dataset3PreventOverlap','1')
        itl.add_variable('dataset3ColoringType','both')

#===================================
    if os.path.exists(dataset4):
        itl.add_variable('dataset4File', dataset4)
        itl.add_variable('dataset4Label', 'cluster')
        itl.add_variable('dataset4Separator','comma')
        itl.add_variable('dataset4Type','colorstrip')
        itl.add_variable('dataset4StripWidth','200')
        itl.add_variable('dataset4PreventOverlap','1')
        itl.add_variable('dataset4ColoringType','both')
        itl.add_variable('dataset4BranchColoringType','both')
        
#itl.add_variable('dataset1BarSizeMax','1')

#===================================
# Check parameters
# itl.print_variables()


#Submit the tree
    print ''
    print 'Uploading the tree.  This may take some time depending on how large the tree is and how much load there is on the itol server'
    good_upload = itl.upload()
    if good_upload == False:
        print 'There was an error:'+itl.comm.upload_output
        sys.exit(1)

#Read the tree ID
    print 'Tree ID: '+str(itl.comm.tree_id)

#Read the iTOL API return statement
    print 'iTOL output: '+str(itl.comm.upload_output)

#Website to be redirected to iTOL tree
    print 'Tree Web Page URL: '+itl.get_webpage()

# Warnings associated with the upload
    print 'Warnings: '+str(itl.comm.warnings)

#Export to pdf
    print 'Exporting to pdf'
    itol_exporter = itl.get_itol_export()
#itol_exporter = itolexport.ItolExport()
#itol_exporter.set_export_param_value('tree','18793532031912684633930')
    itol_exporter.set_export_param_value('format', 'eps')
    itol_exporter.set_export_param_value('displayMode',"circular")
    itol_exporter.set_export_param_value('showBS',"0")
    itol_exporter.set_export_param_value('fontSize',fontsize)
    itol_exporter.set_export_param_value('alignLabels',"1")
    itol_exporter.set_export_param_value('datasetList','dataset1,dataset2,dataset3,dataset4')
    epsfile = outpath + os.sep + pfamid + '-itol.eps'
    pdffile = outpath + os.sep + pfamid + '-itol.pdf'
    jpgfile = outpath + os.sep + pfamid + '-itol.jpg'
    thumbfile = outpath + os.sep + "thumb." + pfamid + '-itol.jpg'
    itol_exporter.export(epsfile)
    os.system("epstopdf %s" % epsfile )
    os.system("convert %s %s" % (epsfile, jpgfile) )
    os.system("convert -thumbnail 200 %s %s" % (jpgfile, thumbfile))
    print 'exported tree to ',pdffile
#}}}
def Itol_Tree_m1(pfamid, datapath, outpath):#{{{
# TM helices are treated as domains
#Create the Itol class
    itl = Itol.Itol()
#Set the tree file
    tree = datapath + os.sep + pfamid + '.kalignp.tree'
    (dataset1, dataset2, dataset3, dataset4) = ("", "", "", "")
    if not os.path.exists(tree):
        print >> sys.stderr, "tree file %s does not exist. Ignore" %(tree)
        return 1
    t = Tree(tree)
    leaves = t.get_leaves()
    numLeave = len(leaves)

    fontsize = GetFontSize(numLeave)

    dataset1 = datapath + os.sep + pfamid + '.numTM_and_io.txt'
    dataset2 = datapath + os.sep + pfamid + '.cmpclass.colordef.txt'
#    dataset3 = datapath + os.sep + pfamid + '.ntermstate.colordef.txt'
    dataset4 = datapath + os.sep + pfamid + '.cluster.colordef.txt'

    colordeffile = datapath + os.sep + pfamid + '.pfam.colordef.txt'
    branchlabelfile = datapath + os.sep + pfamid + '.branchlabel.txt'

#===================================
    itl.add_variable('treeFile',tree)
    itl.add_variable('treeName','PF00854')
    itl.add_variable('treeFormat','newick')
    if os.path.exists(colordeffile):
        itl.add_variable('colorDefinitionFile', colordeffile)
    if os.path.exists(branchlabelfile):
        itl.add_variable('branchLabelsFile', branchlabelfile)

#===================================
    if os.path.exists(dataset1):
        itl.add_variable('dataset1File',dataset1)
        itl.add_variable('dataset1Label','numTM_and_io')
        itl.add_variable('dataset1Separator','comma')
        itl.add_variable('dataset1Type','multibar')
        itl.add_variable('dataset1PreventOverlap','1')
        itl.add_variable('dataset1Color','#FF0000')

#===================================
    if os.path.exists(dataset2):
        itl.add_variable('dataset2File', dataset2)
        itl.add_variable('dataset2Label', 'cmpclass')
        itl.add_variable('dataset2Separator','comma')
        itl.add_variable('dataset2Type','colorstrip')
        itl.add_variable('dataset2StripWidth','200')
        itl.add_variable('dataset2PreventOverlap','1')
        itl.add_variable('dataset2ColoringType','both')

#===================================
    if os.path.exists(dataset3):
        itl.add_variable('dataset3File', dataset3)
        itl.add_variable('dataset3Label', 'NtermState')
        itl.add_variable('dataset3Separator','comma')
        itl.add_variable('dataset3Type','colorstrip')
        itl.add_variable('dataset3StripWidth','200')
        itl.add_variable('dataset3PreventOverlap','1')
        itl.add_variable('dataset3ColoringType','both')

#===================================
    if os.path.exists(dataset4):
        itl.add_variable('dataset4File', dataset4)
        itl.add_variable('dataset4Label', 'cluster')
        itl.add_variable('dataset4Separator','comma')
        itl.add_variable('dataset4Type','colorstrip')
        itl.add_variable('dataset4StripWidth','200')
        itl.add_variable('dataset4PreventOverlap','1')
        itl.add_variable('dataset4ColoringType','both')
        itl.add_variable('dataset4BranchColoringType','both')
        
#itl.add_variable('dataset1BarSizeMax','1')

#===================================
# Check parameters
# itl.print_variables()


#Submit the tree
    print ''
    print 'Uploading the tree.  This may take some time depending on how large the tree is and how much load there is on the itol server'
    good_upload = itl.upload()
    if good_upload == False:
        print 'There was an error:'+itl.comm.upload_output
        sys.exit(1)

#Read the tree ID
    print 'Tree ID: '+str(itl.comm.tree_id)

#Read the iTOL API return statement
    print 'iTOL output: '+str(itl.comm.upload_output)

#Website to be redirected to iTOL tree
    print 'Tree Web Page URL: '+itl.get_webpage()

# Warnings associated with the upload
    print 'Warnings: '+str(itl.comm.warnings)

#Export to pdf
    print 'Exporting to pdf'
    itol_exporter = itl.get_itol_export()
#itol_exporter = itolexport.ItolExport()
#itol_exporter.set_export_param_value('tree','18793532031912684633930')
    itol_exporter.set_export_param_value('format', 'eps')
    itol_exporter.set_export_param_value('displayMode',"circular")
    itol_exporter.set_export_param_value('showBS',"0")
    itol_exporter.set_export_param_value('fontSize',fontsize)
    itol_exporter.set_export_param_value('alignLabels',"1")
    itol_exporter.set_export_param_value('datasetList','dataset1,dataset2,dataset3,dataset4')
    epsfile = outpath + os.sep + pfamid + '-itol.eps'
    pdffile = outpath + os.sep + pfamid + '-itol.pdf'
    jpgfile = outpath + os.sep + pfamid + '-itol.jpg'
    thumbfile = outpath + os.sep + "thumb." + pfamid + '-itol.jpg'
    itol_exporter.export(epsfile)
    os.system("epstopdf %s" % epsfile )
    os.system("convert %s %s" % (epsfile, jpgfile) )
    os.system("convert -thumbnail 200 %s %s" % (jpgfile, thumbfile))
    print 'exported tree to ',pdffile
#}}}
def Itol_Tree_m_sd1(pfamid, datapath, outpath):#{{{
    """Phylogenetic tree with numTM_io and subfamilies branch coloring
    """
    tree = datapath + os.sep + pfamid + '.kalignp.tree'
    t = Tree(tree)
    leaves = t.get_leaves()
    lst_leaves_name = []
    for leaf in leaves:
        lst_leaves_name.append(leaf.name)
    numLeave = len(lst_leaves_name)
# read subfamily definition
    subfamfile = "%s/%s.subfamilies"%(datapath, pfamid)
    subfam_idlist = []
    subfamDict = myfunc.Read_subfamily(subfamfile, subfam_idlist)
    numSubFam = len(subfam_idlist)

# create subfamily branch color definition file
    subfam_colordef_file = "%s/%s.subfamilies.colordef.txt"%(outpath, pfamid)
    lst_color = list(blue.range_to(red,numSubFam))
    color_dict = {}
    for i in xrange(numSubFam):
        famid = subfam_idlist[i]
        color = lst_color[i].get_hex_l()
        color_dict[famid] = lst_color[i].get_hex_l()
    myfunc.WriteSubFamColorDef(subfam_colordef_file, subfamDict, lst_leaves_name, color_dict)

#Create the Itol class
    itl = Itol.Itol()
#Set the tree file
    (dataset1, dataset2, dataset3, dataset4) = ("", "", "", "")
    if not os.path.exists(tree):
        print >> sys.stderr, "tree file %s does not exist. Ignore" %(tree)
        return 1

    fontsize = GetFontSize(numLeave)

    dataset1 = datapath + os.sep + pfamid + '.numTM_and_io.txt'
    colordeffile = subfam_colordef_file

    if os.path.exists(colordeffile):
        itl.add_variable('colorDefinitionFile', colordeffile)
        itl.add_variable('colorDefinitionLabel', "Subfamilies")

#===================================
    itl.add_variable('treeFile',tree)
    itl.add_variable('treeName','SD1')
    itl.add_variable('treeFormat','newick')
#===================================
    if os.path.exists(dataset1):
        itl.add_variable('dataset1File',dataset1)
        itl.add_variable('dataset1Label','numTM_and_io')
        itl.add_variable('dataset1Separator','comma')
        itl.add_variable('dataset1Type','multibar')
        itl.add_variable('dataset1PreventOverlap','1')
        itl.add_variable('dataset1Color','#FF0000')
#===================================
# Check parameters
# itl.print_variables()

#Submit the tree
    print ''
    print 'Uploading the tree.  This may take some time depending on how large the tree is and how much load there is on the itol server'
    good_upload = itl.upload()
    if good_upload == False:
        print 'There was an error:'+itl.comm.upload_output
        sys.exit(1)

#Read the tree ID
    print 'Tree ID: '+str(itl.comm.tree_id)

#Read the iTOL API return statement
    print 'iTOL output: '+str(itl.comm.upload_output)

#Website to be redirected to iTOL tree
    print 'Tree Web Page URL: '+itl.get_webpage()

# Warnings associated with the upload
    print 'Warnings: '+str(itl.comm.warnings)

#Export to pdf
    itol_exporter = itl.get_itol_export()
#itol_exporter = itolexport.ItolExport()
#itol_exporter.set_export_param_value('tree','18793532031912684633930')
    itol_exporter.set_export_param_value('format', 'eps')
    itol_exporter.set_export_param_value('displayMode',"circular")
    itol_exporter.set_export_param_value('showBS',"0")
    itol_exporter.set_export_param_value('fontSize',fontsize)
    itol_exporter.set_export_param_value('alignLabels',"1")
    itol_exporter.set_export_param_value('datasetList','dataset1')
    extname = "-itol-sd1"
    epsfile = outpath + os.sep + pfamid + extname + '.eps'
    pdffile = outpath + os.sep + pfamid + extname + '.pdf'
    jpgfile = outpath + os.sep + pfamid + extname + '.jpg'
    pngfile = outpath + os.sep + pfamid + extname + '.png'
    thumbfile = outpath + os.sep + "thumb." + pfamid + extname + '.jpg'
    itol_exporter.export(epsfile)
    os.system("epstopdf %s" % epsfile )
    os.system("convert %s %s" % (epsfile, jpgfile) )
    os.system("convert -thumbnail 200 %s %s" % (jpgfile, thumbfile))
    print 'exported tree to ',pdffile
#}}}
def Itol_Tree_m_sd2(pfamid, datapath, outpath):#{{{
    tree = datapath + os.sep + pfamid + '.kalignp.tree'
    t = Tree(tree)
    leaves = t.get_leaves()
    lst_leaves_name = []
    for leaf in leaves:
        lst_leaves_name.append(leaf.name)
    numLeave = len(lst_leaves_name)
    leaves_name_set = set(lst_leaves_name)
# read seqlen file
    seqlenfile = "%s/%s.seqlen.txt"%(datapath, pfamid)
    seqlen_dict = myfunc.ReadSeqLengthDict(seqlenfile)

# read subfamily definition
    domain_idlist = []
    domainfile = "%s/%s.mdp"%(datapath, pfamid)
    domain_dict = myfunc.Read_domain_sd(domainfile, domain_idlist)
    domain_colordef_file = "%s/%s.mdp.colordef.txt"%(datapath, pfamid)
    WriteDomainColorDefFile(domain_colordef_file, domain_dict, domain_idlist, seqlen_dict, leaves_name_set)

#Create the Itol class
    itl = Itol.Itol()
#Set the tree file
    (dataset1, dataset2, dataset3, dataset4) = ("", "", "", "")
    if not os.path.exists(tree):
        print >> sys.stderr, "tree file %s does not exist. Ignore" %(tree)
        return 1

    fontsize = GetFontSize(numLeave)

    dataset1 = domain_colordef_file
#     colordeffile = subfam_colordef_file
#     if os.path.exists(colordeffile):
#         itl.add_variable('colorDefinitionFile', colordeffile)
#         itl.add_variable('colorDefinitionLabel', "Subfamilies")
#===================================
    itl.add_variable('treeFile',tree)
    itl.add_variable('treeName','SD2')
    itl.add_variable('treeFormat','newick')
#===================================
    if os.path.exists(dataset1):
        itl.add_variable('dataset1File',dataset1)
        itl.add_variable('dataset1Label','Domain architecture')
        itl.add_variable('dataset1Separator','comma')
        itl.add_variable('dataset1Type','domains')
        itl.add_variable('dataset1ProtSizeMax','1000')
        itl.add_variable('dataset1PreventOverlap','1')
        itl.add_variable('dataset1CirclesSpacing','100')
#===================================
# Check parameters
# itl.print_variables()

#Submit the tree
    print ''
    print 'Uploading the tree.  This may take some time depending on how large the tree is and how much load there is on the itol server'
    good_upload = itl.upload()
    if good_upload == False:
        print 'There was an error:'+itl.comm.upload_output
        sys.exit(1)

#Read the tree ID
    print 'Tree ID: '+str(itl.comm.tree_id)

#Read the iTOL API return statement
    print 'iTOL output: '+str(itl.comm.upload_output)

#Website to be redirected to iTOL tree
    print 'Tree Web Page URL: '+itl.get_webpage()

# Warnings associated with the upload
    print 'Warnings: '+str(itl.comm.warnings)

#Export to pdf
    itol_exporter = itl.get_itol_export()
#itol_exporter = itolexport.ItolExport()
#itol_exporter.set_export_param_value('tree','18793532031912684633930')
    itol_exporter.set_export_param_value('format', 'eps')
    itol_exporter.set_export_param_value('displayMode',"circular")
    itol_exporter.set_export_param_value('showBS',"0")
    itol_exporter.set_export_param_value('fontSize',fontsize)
    itol_exporter.set_export_param_value('alignLabels',"1")
    itol_exporter.set_export_param_value('datasetList','dataset1')
    extname = "-itol-sd2"
    epsfile = outpath + os.sep + pfamid + extname + '.eps'
    pdffile = outpath + os.sep + pfamid + extname + '.pdf'
    jpgfile = outpath + os.sep + pfamid + extname + '.jpg'
    pngfile = outpath + os.sep + pfamid + extname + '.png'
    thumbfile = outpath + os.sep + "thumb." + pfamid + extname + '.jpg'
    itol_exporter.export(epsfile)
    os.system("epstopdf %s" % epsfile )
    os.system("convert %s %s" % (epsfile, jpgfile) )
    os.system("convert -thumbnail 200 %s %s" % (jpgfile, thumbfile))
    print 'exported tree to ',pdffile
#}}}
def Itol_Tree_m_sd3(pfamid, datapath, outpath):#{{{
    """Phylogenetic tree with species definition
    the Kindom use branch colordefinition, and others using color strips
    """
    tree = datapath + os.sep + pfamid + '.kalignp.tree'
    t = Tree(tree)
    leaves = t.get_leaves()
    lst_leaves_name = []
    for leaf in leaves:
        lst_leaves_name.append(leaf.name)
    numLeave = len(lst_leaves_name)
    leaves_name_set = set(lst_leaves_name)

# read species definition
    speciesfile = "%s/%s.species"%(datapath, pfamid)
    speciesDict = myfunc.Read_species_sd(speciesfile)

# create branch color definition file for kingdom
    lst_kingdom = ["Archaea","Bacteria", "Eukaryota" ]
    lst_color_kingdom = ["#ff0000", "#0066ff","#cc6600"]
    species_colordef_file = "%s/%s.kingdom.colordef.txt"%(outpath, pfamid)
    color_dict_kingdom = {}
    this_speciesDict = {}
    for seqid in speciesDict:
        speciesname = speciesDict[seqid][0]
        this_speciesDict[seqid] = speciesname
    for i in xrange(len(lst_kingdom)):
        idd = lst_kingdom[i]
        color_dict_kingdom[idd] = lst_color_kingdom[i]
    myfunc.WriteKingdomColorDefFile(species_colordef_file, this_speciesDict, leaves_name_set, color_dict_kingdom)

# generate the next three levels of classification
    for level in [1,2,3]:
        outfile = "%s/%s.species.level_%d.txt"%(outpath, pfamid, level)
        this_speciesDict = {}
        speciesIDSet = set([])
        for seqid in speciesDict:
            try:
                speciesname = speciesDict[seqid][level]
                speciesIDSet.add(speciesname)
                this_speciesDict[seqid] = speciesname
            except IndexError, KeyError:
                pass
        color_dict = {}
        lst_color = list(blue.range_to(red,len(speciesIDSet)))
        lst_speciesID = list(speciesIDSet)
        for i in xrange(len(lst_speciesID)):
            idd = lst_speciesID[i]
            color_dict[idd] = lst_color[i].get_hex_l()
        myfunc.WriteSpeciesColorStripDefFile(outfile, this_speciesDict, leaves_name_set, color_dict)

#Create the Itol class
    itl = Itol.Itol()
#Set the tree file
    (dataset1, dataset2, dataset3, dataset4) = ("", "", "", "")
    if not os.path.exists(tree):
        print >> sys.stderr, "tree file %s does not exist. Ignore" %(tree)
        return 1

    fontsize = GetFontSize(numLeave)

    dataset1 = "%s/%s.species.level_%d.txt"%(outpath, pfamid, 1)
    dataset2 = "%s/%s.species.level_%d.txt"%(outpath, pfamid, 2)
    dataset3 = "%s/%s.species.level_%d.txt"%(outpath, pfamid, 3)
    colordeffile = species_colordef_file

    if os.path.exists(colordeffile):
        itl.add_variable('colorDefinitionFile', colordeffile)
        itl.add_variable('colorDefinitionLabel', "Kingdom")

#===================================
    itl.add_variable('treeFile',tree)
    itl.add_variable('treeName','SD3')
    itl.add_variable('treeFormat','newick')
#===================================
    if os.path.exists(dataset1):
        itl.add_variable('dataset1File',dataset1)
        itl.add_variable('dataset1Label','Phylum')
        itl.add_variable('dataset1Separator','comma')
        itl.add_variable('dataset1Type','colorstrip')
        itl.add_variable('dataset1StripWidth','100')
        itl.add_variable('dataset1ColoringType','both')
        itl.add_variable('dataset1PreventOverlap','1')
#===================================
    if os.path.exists(dataset2):
        itl.add_variable('dataset2File',dataset2)
        itl.add_variable('dataset2Label','Class')
        itl.add_variable('dataset2Separator','comma')
        itl.add_variable('dataset2Type','colorstrip')
        itl.add_variable('dataset2StripWidth','100')
        itl.add_variable('dataset2ColoringType','both')
        itl.add_variable('dataset2PreventOverlap','1')
#===================================
    if os.path.exists(dataset3):
        itl.add_variable('dataset3File',dataset3)
        itl.add_variable('dataset3Label','Order')
        itl.add_variable('dataset3Separator','comma')
        itl.add_variable('dataset3Type','colorstrip')
        itl.add_variable('dataset3StripWidth','100')
        itl.add_variable('dataset3ColoringType','both')
        itl.add_variable('dataset3PreventOverlap','1')
#===================================
# Check parameters
# itl.print_variables()

#Submit the tree
    print ''
    print 'Uploading the tree.  This may take some time depending on how large the tree is and how much load there is on the itol server'
    good_upload = itl.upload()
    if good_upload == False:
        print 'There was an error:'+itl.comm.upload_output
        sys.exit(1)

#Read the tree ID
    print 'Tree ID: '+str(itl.comm.tree_id)

#Read the iTOL API return statement
    print 'iTOL output: '+str(itl.comm.upload_output)

#Website to be redirected to iTOL tree
    print 'Tree Web Page URL: '+itl.get_webpage()

# Warnings associated with the upload
    print 'Warnings: '+str(itl.comm.warnings)

#Export to pdf
    itol_exporter = itl.get_itol_export()
#itol_exporter = itolexport.ItolExport()
#itol_exporter.set_export_param_value('tree','18793532031912684633930')
    itol_exporter.set_export_param_value('format', 'eps')
    itol_exporter.set_export_param_value('displayMode',"circular")
    itol_exporter.set_export_param_value('showBS',"0")
    itol_exporter.set_export_param_value('fontSize',fontsize)
    itol_exporter.set_export_param_value('alignLabels',"1")
    itol_exporter.set_export_param_value('datasetList','dataset1,dataset2,dataset3')
    extname = "-itol-sd3"
    epsfile = outpath + os.sep + pfamid + extname + '.eps'
    pdffile = outpath + os.sep + pfamid + extname + '.pdf'
    jpgfile = outpath + os.sep + pfamid + extname + '.jpg'
    pngfile = outpath + os.sep + pfamid + extname + '.png'
    thumbfile = outpath + os.sep + "thumb." + pfamid + extname + '.jpg'
    itol_exporter.export(epsfile)
    os.system("epstopdf %s" % epsfile )
    os.system("convert %s %s" % (epsfile, jpgfile) )
    os.system("convert -thumbnail 200 %s %s" % (jpgfile, thumbfile))
    print 'exported tree to ',pdffile
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    datapath = "."
    outpath = './'
    idList = []
    idListFile = ''

    i = 1;
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg=False;
            idList.append(sys.argv[i])
            i = i + 1;
        elif sys.argv[i] == "--":
            isNonOptionArg=True;
            i = i + 1;
        elif sys.argv[i][0] == "-":
            if sys.argv[i] in [ "-h", "--help"]:
                PrintHelp();
                return 1
            elif sys.argv[i] in [ "-datapath", "--datapath"]:
                datapath = sys.argv[i+1]
                i += 2;
            elif argv[i] in [ "-m", "--m", "-method", "--method"]:
                g_params['method'], i = myfunc.my_getopt_str(sys.argv, i)
            elif sys.argv[i] in [ "-l", "--l"]:
                idListFile = sys.argv[i+1]
                i = i + 2;
            elif sys.argv[i] in ["-outpath", "--outpath"]:
                outpath = sys.argv[i+1];
                i = i + 2;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                return 1
        else:
            idList.append(sys.argv[i]);
            i+=1;

    if idListFile != "":
        idList += myfunc.ReadIDList(idListFile)
    if len(idList) > 0:
        os.system("mkdir -p %s"%outpath)
        cnt = 0
        for pfamid in idList:
            print "================== ", cnt , pfamid, " ===================="
            if g_params['method'] == "0":
                Itol_Tree_m0(pfamid, datapath, outpath)
            elif g_params['method'] == "1":
                Itol_Tree_m1(pfamid, datapath, outpath)
            elif g_params['method'] == "sd1":
                Itol_Tree_m_sd1(pfamid, datapath, outpath)
            elif g_params['method'] == "sd2":
                Itol_Tree_m_sd2(pfamid, datapath, outpath)
            elif g_params['method'] == "sd3":
                Itol_Tree_m_sd3(pfamid, datapath, outpath)
            cnt += 1
#}}}
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['method'] = "0"
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))


