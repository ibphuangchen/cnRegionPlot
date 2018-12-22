#!/usr/bin/env Rscript
#Author: Chen Huang ibphuangchen@gmail.com

##all the functions needed
seg2Ranges = function(totalChrLengthDf){
  #this is only for one samples
  rangeObj = GRanges(seqnames = totalChrLengthDf$chromosome, 
                     ranges = IRanges(start = 0, end = totalChrLengthDf$length),
                     strand = Rle(strand("*"), nrow(totalChrLengthDf)))
  rangeObj
}

fillwithNA = function(segDf, totalChrRangeObj){
  #this function fill the segments that are not shown in the segmentation result with NA
  #segDf is a df, should only come from one sample
  segDf$chromosome = as.character(segDf$chromosome)
  #harmonize the chromosome to be the style of "Y"
  if(grepl(pattern = "chr",segDf$chromosome[1]))
    segDf$chromosome = substr(segDf$chromosome, 4, nchar(segDf$chromosome))
  segDf$chromosome[segDf$chromosome == "23"] = "X"
  segDf$chromosome[segDf$chromosome == "24"] = "Y"
  
  segObj = GRanges(seqnames = segDf$chromosome, 
                   ranges = IRanges(start = segDf$start, end = segDf$end),
                   strand = Rle(strand("*"), nrow(segDf)),
                   log2 = segDf$log2)
  ##the segmenets not covered by the segObj
  naRanges = setdiff(totalChrRangeObj, segObj)
  naRanges$log2 = NA
  #combining naRanges to segObj
  segObj = c(segObj, naRanges)
  segObj = sort(segObj)
  segObj = as.data.frame(segObj)
  colnames(segObj)[1] = "chromosome"
  segObj
}

fillwithNA.2 = function(segDf, totalChrRangeObj){
  #this version works for multiple samples
  if(length(unique(segDf$sample)) == 1)
    fillwithNA(segDf, totalChrRangeObj)
  else{
    segList = split(segDf, f = segDf$sample)
    segList = lapply(segList, function(x)fillwithNA(segDf = x, totalChrRangeObj = totalChrRangeObj))
    newSegDf = do.call(rbind, segList)
    newSegDf$sample = rep(names(segList), sapply(segList, nrow))
    newSegDf
  }
}

selectRange = function(segDf, chr, start = NULL, end = NULL, chrLength){
  #chr, start, end are all vectors here
  #segDf should be from one sample
  require(GenomicRanges)
  chr = as.character(chr)
  if(is.null(start)){
    start = rep(0, length(chr))
  }
  if(is.null(end)){
    end = chrLength[chr,"length"]
  }
  if(any(start > chrLength[chr,"length"]) | any(start < 0) | any(end > chrLength[chr,"length"]))
    stop("Some displaying ranges exceed the chromosome length, please corret it.")
  
  #convert segDf to Grange Obj
  if(grepl(pattern = "chr",segDf$chromosome[1]))
    segDf$chromosome = substr(segDf$chromosome, 4, nchar(segDf$chromosome))
  segDf$chromosome[segDf$chromosome == "23"] = "X"
  segDf$chromosome[segDf$chromosome == "24"] = "Y"
  segObj = GRanges(seqnames = segDf$chromosome, 
                   ranges = IRanges(start = segDf$start, end = segDf$end),
                   strand = Rle(strand("*"), nrow(segDf)),
                   log2 = segDf$log2)
  
  findOverlapRegions = function(singleChr, start, end, queryRegion){
    refRangeObj = GRanges(seqnames = singleChr, 
                          ranges = IRanges(start = start, end = end),
                          strand = Rle(strand("*"), length(singleChr)))
    overlapIndex = findOverlaps(query = queryRegion, subject = refRangeObj, type = 'any')
    overlapIndex = as.data.frame(overlapIndex)
    overlapDf = as.data.frame(queryRegion[overlapIndex$queryHit])
    overlapDf[1,"start"] = start
    overlapDf[1, "width"] = overlapDf[1,"end"] -overlapDf[1,"start"]+1
    overlapDf[nrow(overlapDf), "end"] = end
    overlapDf[nrow(overlapDf), "width"] = overlapDf[nrow(overlapDf),"end"] -overlapDf[nrow(overlapDf),"start"]+1
    overlapDf
  }
  if(length(chr) == 1) subSegDf = findOverlapRegions(chr, start, end, segObj)
  else{
    subSegDf = lapply(1:length(chr), function(x) findOverlapRegions(chr[x], start[x], end[x], segObj))
    subSegDf = do.call(rbind, subSegDf)
  }
  colnames(subSegDf)[1] = "chromosome"
  subSegDf
}

selectRange.2 = function(segDf, samples = NULL, chr = NULL, start = NULL, end = NULL, chrLength){
  #this is the version that can do multiple samples
  if(is.null(chr))
    chr = chrLength[,"chromosome"]
  if(length(unique(segDf$sample)) == 1)
    selectRange(segDf, chr, start, end, chrLength)
  else{
    segList = split(segDf, f = segDf$sample)
    if(!is.null(samples))
      segList = segList[samples]
    segList = lapply(segList, function(x)selectRange(x, chr, start, end, chrLength))
    newSegDf = do.call(rbind, segList)
    newSegDf$sample = rep(names(segList), sapply(segList, nrow))
    newSegDf
  }
}

plotCNV = function(segDf, outputFile, fileFormat){
  require(ggplot2, quietly = TRUE)
  require(scales, quietly = TRUE)
  segDf$name = paste(segDf$sample, segDf$chromosome, segDf$start, sep = '-')
  segDf$width = segDf$end - segDf$start
  segDf$log2 = ifelse(segDf$log2 > 1, 1, segDf$log2)
  segDf$log2 = ifelse(segDf$log2 < -1, -1, segDf$log2)
  samples = unique(segDf$sample)
  plotCNV = ggplot(segDf, aes(1, width, fill=log2)) + 
    geom_bar(stat="identity")+
    scale_fill_gradient2(high = muted("red"), low = muted("blue"))+
    geom_hline(yintercept = 0, size = 1)+
    scale_y_continuous(expand = c(0,0))+
    facet_grid(sample ~ chromosome, space = "free",scales = "free") +
    theme(axis.title =element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.spacing = unit(0, "lines"),
          strip.background = element_rect(colour="black", fill="white", 
                                          size=5, linetype="blank"),
          strip.text.y = element_text(angle = 0, size = 0.1*length(samples)+5),
          strip.text.x = element_text(size = 0.1*length(samples)+5),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom",
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5, angle = 90, hjust = 0.1, vjust = 0.5)) +
    coord_flip()
  baseHigh = 1
  ggsave(plot = plotCNV, filename = outputFile, device = fileFormat, 
         dpi = 300, width = 20, units = "cm", height = baseHigh*(length(samples)+2))
}

# create parser object
suppressPackageStartupMessages(library("argparse"))
parser = ArgumentParser(description= "Author: Chen Huang ibphuangchen@gmail.com")
parser$add_argument("file", nargs=1, help="Segment file to be displayed")
parser$add_argument("output", nargs=1, help="File of the output plot")
parser$add_argument("-S", "--samples", help="Samples to be displayed [default ALL]",
                  metavar="sample1,sample2,...", default = NULL)
parser$add_argument("-c", "--chromosome", help="Chromosomes to be displayed [default ALL]",
                    metavar="chromosome1, chromsome2, ...", default = NULL)

parser$add_argument("-s", "--start", help="Start position of the chromsome [default NULL]",
                    metavar = "start1, start2, ...", default = NULL)

parser$add_argument("-e", "--end", help="End position of the chromsome [default NULL]",
                    metavar = "end1, end2, ...", default = NULL)

parser$add_argument("-g", "--genome", help="Genome version [default hg19]",
                    metavar = "hg19|hg38", default = "hg19")

args <- parser$parse_args()

#load libraries
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("ggplot2"))

#parse and examine all the parameters
checkFlag = FALSE
if(!is.null(args$start)){
  startVec = strsplit(args$start,",")[[1]]
  startVec = as.numeric(startVec)
  checkFlag = TRUE
}else
  startVec = NULL
if(!is.null(args$end)){
  endVec = strsplit(args$end,",")[[1]]
  endVec = as.numeric(endVec)
  checkFlag = TRUE
}else
  endVec = NULL
if (!is.null(args$samples)) {
  sampleVec = strsplit(args$samples,",")[[1]]
}else
  sampleVec = NULL
if(!is.null(args$chromosome)){
  chromVec = strsplit(args$chromosome,",")[[1]]
  checkFlag = TRUE
}else chromVec = NULL

ChrLengthDf = NULL
if(args$genome == "hg19"){
  ChrLengthDf = readRDS('./hg19ChrLength.RDS')
}else if(args$genome == "hg38"){
  ChrLengthDf = readRDS('./hg38ChrLength.RDS')
}else
  stop("Current only support genome version hg19 and hg38")
chrLengthObj = seg2Ranges(ChrLengthDf)

if(checkFlag == TRUE &(length(startVec)!=length(chromVec) | length(endVec)!=length(chromVec)))
  stop("Error: You need to specify equal number of starts and ends to the chromsomes")
output = args$output

if(grepl(output, pattern = "\\.")) {
  outputFormat = gsub(output, pattern = ".*\\.(.+)", replacement = "\\1")
  outputFormat = tolower(outputFormat)
  }else
    outputFormat = "png"

if(!outputFormat %in% c("pdf", "jpeg", "tiff", "png", "bmp", "svg"))
    stop("Unsupported output format. Currently support:\n
      \t\t\"pdf\", \"jpeg\", \"tiff\", \"png\", \"bmp\", \"svg\"")


cat("Reading in the segment file\n")
segFile = read.table(args$file, sep = "\t", header = TRUE, stringsAsFactors = F)

if(any(!c("sample","chromosome","start","end","log2") %in% colnames(segFile))){
  stop("Error: the segFiles should contain 5 columns:\n
       \t\"sample\",\"chromsome\",\"start\",\"end\",\"log2\"")
}
cat("Reading input file DONE\n")

cat("Preprocessing the segment file\n")
segDf = fillwithNA.2(segDf = segFile, totalChrRangeObj = chrLengthObj)
cat("Preprocessing the segment file DONE\n")

cat("Selecting genome range for plotting\n")
segDf = selectRange.2(segDf = segDf, samples = sampleVec, 
                      start = startVec, end = endVec, chr = chromVec, chrLength = ChrLengthDf)
cat("Selecting genome range for plotting DONE\n")
cat("Plotting figure\n")
ggplotObj = plotCNV(segDf, outputFile = output, fileFormat = outputFormat)
cat("Plotting figure DONE\n")

