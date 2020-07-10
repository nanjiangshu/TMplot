import os
import sys
import myfunc
rundir = os.path.dirname(os.path.realpath(__file__))
# prints whether python is version 3 or not
python_version = sys.version_info.major
if python_version == 3:
    print("is python 3")
else:
    print("not python 3")

from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont


basedir = os.path.realpath("%s/../"%(rundir))

progname=os.path.basename(sys.argv[0])
general_usage = """ 
usage: %s TESTMODE options
"""%(sys.argv[0])

numArgv = len(sys.argv)
if numArgv <= 1:
    print(general_usage)
    sys.exit(1)
TESTMODE=sys.argv[1]

g_params = {}

if TESTMODE == "loadpil":
    g_params['font_dir'] = "%s/../fonts/truetype/ttf-dejavu/"%(rundir)
    g_params['font_size'] = 16
    fontpath = g_params['font_dir'] + "DejaVuSerif.ttf"
    print (fontpath)
    g_params['fntTMbox_label'] = ImageFont.truetype(fontpath, 10)

if TESTMODE == "getgapposition":
    topo = sys.argv[2]
    posGAP = myfunc.GetGapPosition(topo)
    print(posGAP)

if TESTMODE == "readfasta":
    seqfile = sys.argv[2]
    (idList, annoList, seqList) = myfunc.ReadFasta(seqfile)
    print(idList)
    print(seqList)


