## Visualize alignment of TM proteins with both TM alignment and sequence alignment

### Input file
* Pairwise alignment of membrane topology [glt_hct_full.topomfa](../test/glt_hct_full.topomfa) and corresponding sequences [glt_hct_full.aln](../test/glt_hct_full.aln)

### Requirements

In order to run the script, you will need to install the following tools apart
from the packages in the [requirements](../requirements.txt)

* `wkhtmltopdf` with patched QT
    for Debian, you will need to install `xvfb`
    for MacOS, it is fine to install wkhtmltopdf just by Homebrew

*  Imagemagick

   Note that on Ubuntu 18, you will need to change the default proxy, that is

    change
    ```
    <policy domain="resource" name="width" value="16KP"/>
    <policy domain="resource" name="height" value="16KP"/>
    <policy domain="resource" name="disk" value="1GiB"/>
    ```
    to

    ```
    <policy domain="resource" name="width" value="64KP"/>
    <policy domain="resource" name="height" value="64KP"/>
    <policy domain="resource" name="disk" value="4GiB”/>
    ```
    in the file `/etc/ImageMagick-6/policy.xml`

*  pdflatex

### Command

    python src/makeTMplot.py -seqaln test/glt_hct_full.aln -topaln test/glt_hct_full.topomfa -breakTM -outpath test/outdir


After successful run of the above command, the result will be output to `test/outdir/glt_hct_full.seqtopaln.pdf`.

<img src="../examples/example_images/glt_hct_full.seqtopaln.png">



