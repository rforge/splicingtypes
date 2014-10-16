##################################                           
## SplicingTypesAnno
##  - global.r: global variables
##################################

library(WriteXLS)
library(gdata)
library(hwriter)
library(SortableHTMLTables)
library(Rsamtools)

## new env
splicingENV <- new.env()

## variables
assign("PROJECT_NAME", "Splicing Types", splicingENV)

assign("html.header", 
"<html xmlns:mml=\"http://www.w3.org/1998/Math/MathML\">
<head>
<title> Splicing Events </title>
<link rel=stylesheet href=\"img/report.css\" type=text/css>
</head>
<body ><br>
<div align=\"center\"><img src=\"img/splice.jpg\" WIDTH=\"25%\"></div>",
, splicingENV)

assign("outputRawHTML", "sample.html", splicingENV) 
assign("outputSummaryHTML", "summary.html", splicingENV)
assign("outputGeneHTML", "gene.html", splicingENV)
assign("outputReadHTML", "read.html", splicingENV)
assign("html.index.name", "index.html", splicingENV)

#assign("outputDir", "html", splicingENV)
assign("raw.dir", "raw", splicingENV)
assign("intron.dir", "intron", splicingENV)
assign("summary.dir", "summary", splicingENV)
assign("gene.dir", "gene", splicingENV)
assign("read.dir","read" , splicingENV)
assign("bed.dir", "bedfiles", splicingENV)
#assign("novelFind.dir", "novelFind", splicingENV)
#assign("novel.df", "NULL", splicingENV)

assign("RI1.dir", "intronRentionType1", splicingENV)
assign("RI2.dir", "intronRentionType2", splicingENV)
assign("ES1.dir","exonSkippingType1", splicingENV)
assign("ES2.dir","exonSkippingType2", splicingENV)
assign("ADLEFT1.dir", "alterDonorType1", splicingENV)
assign("ADLEFT2.dir", "alterDonorType2", splicingENV)
assign("ADRIGHT1.dir", "alterAcceptorType1", splicingENV)
assign("ADRIGHT2.dir", "alterAcceptorType2", splicingENV)
assign("ADBOTH1.dir", "alterBothSiteType1", splicingENV)
assign("ADBOTH2.dir", "alterBothSiteType2", splicingENV)

## important assignment
assign("table.list", list(), splicingENV)
assign("spData", list(), splicingENV)

assign("gene.GRange", NULL, splicingENV)
assign("exonGRange", NULL, splicingENV)
assign("intronGRange", NULL, splicingENV)
assign("splicingLink", NULL, splicingENV)

assign("dataSum", 
      list(ri.1=NULL, ri.2=NULL, es.1=NULL, es.2=NULL, 
      adleft.1=NULL, adleft.2=NULL, adright.1=NULL, adright.2=NULL, 
      adboth.1=NULL, adboth.2=NULL), 
      splicingENV)

#assign("geneSum", 
#      list(nonOverlap=NULL, overlap=NULL), 
#      splicingENV)

assign("readSum", 
      list(totalRead=NULL, 
           qualityJunRead=NULL,
           qualityNonJunRead=NULL,
           sampleLevel=NULL
           ), 
      splicingENV)

