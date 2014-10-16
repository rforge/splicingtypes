\name{translateGTF}
\Rdversion{0.99}

\alias{translateGTF}
\alias{translateGTF,character-method}

\title{Translating gtf/gff file into gene, exon and intron structure}

\description{
Translate gtf/gff file into gene, exon and intron structure as the \link[GenomicRanges]{GRanges}
objects.
}

\usage{
translateGTF(gtfFile,   
             gtfColnames=c("chr", "source", "feature", "start", "end",
                           "score", "strand", "frame", "geneName"),
             gtfChrFormat="chr", 
             gtfGeneLabel="gene_id",
             gtfGeneValue= "g", 
             gtfTranscriptLabel=NULL,                         
             exonString="exon",
             selectGenes=NULL,
             geneOverlap="both", 
             exonOverlap="yes",
             exonBoundary="no" 
             )
}

\arguments{
  \item{gtfFile}{A gtf or gff file.}
  \item{gtfColnames}{The column names for gtf or gff file. "chr", "start", "end", "strand", "geneName" must be provided and can not be changed.}  
  \item{gtfChrFormat}{The chromsome format in the gtf or gff file. It takes four different
  formats, including "chr" (default), "Chr", "chr0", "Chr0", "1" (number represents chromosome). }

  \item{gtfGeneLabel}{The gene label in the gtf or gff file. Generally, for gff3 format (UCSC),
  it is "gene_id" (default); for gtf format, it is "Parent".}
  \item{gtfGeneValue}{The gene value format, which matches the gene label in the gtf or gff file. Generally, for gff3 format (UCSC),
  it is "g" (default), such as "Pnpla3" for mouse; for gtf format, it is "g.1", such as 
  "AT1G01010.1" for arabidopsis. }  

  \item{gtfTranscriptLabel}{The transcript label for identifying isoforms. 
   Generally, for gff3 format (UCSC),
  it is "transcript_id" (default); for gtf format, it is ".1" after gene value, such as 
  "AT1G01010.1" for arabidopsis. ".1" is for different version of a gene product based on TAIR Nomenclature Guidelines.}  

  \item{exonString}{The string used in gtf or gff file to respresent exon. The default value is "exon".}  
  \item{selectGenes}{The selected gene symbols. It takes two values: NULL or a character vector. If NULL is selected,
  all genes will be processed. The select gene symbols must match \code{gtfGeneValue} in the gtf or gff file.}  

  \item{geneOverlap}{The status of gene overlapping. It takes three values:
    \code{no}: Only genes that do not overlap with other genes are processed; 
    \code{yes}: only genes that overlap with other genes are procssed;
    \code{both}: All genes are processed. }  
  \item{exonOverlap}{The status of exon overlapping. Some exons for the same genes are overlapped. It takes two values: TRUE or FALSE.
    \code{no}: The overlapped exons for the same gene are reduced to one single exon. This reduced exon covers the related overlapping exons;
    \code{yes}: The overlapped exons are not reduced.}  
  \item{exonBoundary}{The string of chrosome number, start or end, and strand. Optional parameter for \code{splicingReport}.
    \code{no}: default parameter. The exon boundaries are not calcuated;
    \code{yes}: The exon boundaires are calculated.}   
}

\details{
  gtf/gff format: With the same format, there are still some subtle annotation difference among varous species. 
  \code{gtfChrFormat}, \code{gtfGeneLabel}, \code{gtfGeneValue},
  \code{gtfTranscriptLabel}, \code{exonString} are provided to address these changes.
  The default values are based on the UCSC format.
  
  Intron: Many transcripts for the same gene have the same or overlapping exons.
  To reduce ambiguity, introns are those regions that do not overlap with any known exons.
  
  Output: the output includes gene, exon and intron for non-overlapping genes; 
  and gene2, exon2 and intron2 for overlapping genes if \code{geneOverlap} is set as "both".
   Splicing connections are
  also provided as a character vector when \code{gtfTranscriptLabel} is specified.
}

\value{
  A single character for gtf or gff file name.
  
  A character vector respresenting colname names of gtf or gff file.
  
  A single character representing chromosome format in gtf or gff file.
  
  A single character representing gene lable in gtf or gff file.
  
  A single character representing gene value, which matches gene label in gtf or gff file.
  
  A single character representing transcript label in gtf or gff file.
  
  A single character representing exon string in gtf or gff file.
  
  A character vector respresenting selected gene symbols in gtf or gff file.

  A single character representing gene overlapping status.
  
  A single character representing exon overlapping status.

}

\author{X. Sun}

\seealso{

  \code{\link{splicingGene}}
}

\examples{
library(SplicingTypesAnno)

## Simple use:
mm9.pnpla7.gtfFile <- system.file("extdata", "mm9_Pnpla7.gtf", package="SplicingTypesAnno")
results <- translateGTF(mm9.pnpla7.gtfFile)

## generate splicing links
results <- translateGTF(mm9.pnpla7.gtfFile, gtfTranscriptLabel="transcript_id")
gene <- results$gene
exon <- results$exon
intron <- results$intron
gene

# NOT RUN
## arabidopsis, TAIR 10
# url <- "ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff"
# arab.gff3 <- read.delim(url, header=TRUE)
# write.table(arab.gff3, file="arab.gff3", row.names=FALSE, quote=FALSE, sep="\t")
# results <- translateGTF("arab.gff3", gtfChrFormat="Chr", gtfGeneLabel="Parent",
#             gtfGeneValue= "g.1")
# results$gene
# results$exon
# results$intron 

}

\keyword{methods}
\keyword{classes}
