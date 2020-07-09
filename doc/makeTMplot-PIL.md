## Visualize alignment of TM proteins with both TM alignment and sequence alignment

### Input file
    * Pairwise alignment of membrane topology [glt_hct_full.topomfa](../test/glt_hct_full.topomfa) and corresponding sequences [glt_hct_full.aln](../test/glt_hct_full.aln)

### Command

    python src/makeTMplot.py -seqaln test/glt_hct_full.aln -topaln test/glt_hct_full.topomfa -breakTM -outpath test/outdir


The result will be output to `test/outdir/glt_hct_full.seqtopaln.pdf`.

<img src="../examples/example_images/glt_hct_full.seqtopaln.png">



