## Visualize multiple alignment of transmembrane protein family sulfate symporter transmembrane region

### Input file
    * Multiple alignment of membrane topology [antiport.fam_topomfaspall](../examples/antiport.fam_topomfaspall)
    * DeltaG values for the representative protein of the family [antiport_dg.txt](../examples/antiport_dg.txt)

### Command

    python drawMSATopo.py -m-shrink 0 -method pil -pfm no -text n  -pdg y -showTMidx  -ptag y -sep n  examples/antiport.fam_topomfaspall -h2wratio 0.67


The result will be output to `antiport.png`. The image might be too large to be
display and thus you may need to resize the image by using e.g. [mogrify](https://imagemagick.org/script/mogrify.php) 

<img src="../examples/example_images/antiport.s1600.jpg">



