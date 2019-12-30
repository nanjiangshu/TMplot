import os
import sys
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

g_params = {}
g_params['font_dir'] = "%s/fonts/truetype/ttf-dejavu/"%(rundir)
g_params['font_size'] = 16
fontpath = g_params['font_dir'] + "DejaVuSerif.ttf"
print (fontpath)
g_params['fntTMbox_label'] = ImageFont.truetype(fontpath, 10)
