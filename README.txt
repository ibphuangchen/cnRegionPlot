The R script plot the copy number alteration (log2 ratio) across genome regions.
You can specify the samples, chromosomes and genome regions that you want to display the copy number landscape.

You need to install the following packages to use it:
	ggplot2
	GenomicRanges
	argparse

First, try:
	 Rscript Rscript RegionPlot.R -h

Download the example segment files and chromosome length files and try the following commands.

Rscript RegionPlot.R -S 1,2 ~/mkWES_rm_HighThreshold19/segment-log2-ratio.txt tst.png

Rscript RegionPlot.R -S 1,2 ~/mkWES_rm_HighThreshold19/segment-log2-ratio.txt tst2.png

Rscript RegionPlot.R -S 1,2 -c 1,2,X -s 0,0,0, -e 1000000,2000000,10000000 ~/mkWES_rm_HighThreshold19/segment-log2-ratio.txt tst3.png

Rscript RegionPlot.R -S 1,2 -c 1,2,X -s 0,0,0, -e 1000000,2000000,10000000 ~/mkWES_rm_HighThreshold19/segment-log2-ratio.txt tst4.png

##############################################################################################################################

You are welcome to use this script or modify it, and I will appreciate if you can mention my name (Chen Huang) in the ackownledge part in you paper if you used it.

For help in using this script, contact ibphuangchen@gmail.com  

##############################################################################################################################
