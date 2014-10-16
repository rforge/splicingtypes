\name{splicingCount}
\Rdversion{0.99}

\alias{splicingCount}
\alias{splicingCount,character-method}

\title{Gene-level count quantification}

\description{
Calculate read counts for the gene, exon and intron.
}

\usage{
splicingCount(selectGene,
                      bamFile = "accepted_hits.sort.bam",
                      select.list, 
                      sampleName="sample",
                      sampleID=1, type="any"
                      )
}

\arguments{
  \item{selectGene}{Gene name or gene symbol. This should match 
    \emph{gtfGeneLabel} of \code{\link{translateGTF}} }

  \item{bamFile}{A bam file for analysis, which is the alignment file from alignment software, such as bowtie or tophat.}

  \item{select.list}{A list from \code{\link{translateGTF}}. This list is 
    a collection of GRanges object, including 
    gene, exon, intron details. }
    
  \item{sampleName}{sample name for the bam file.}
  \item{sampleID}{sample ID for the bam file.}
  \item{type}{overlapping types. See \link[IRanges]{countOverlaps} for details.}

  
}

\details{
  The output includes three types of counts: \code{totalCount}: the total
  counts for the selected gene;  
  \code{exonCount}: the counts for all exons;
  \code{intronCount}: the counts for all introns. 
}

\value{
  A single character respresenting a gene symbol in gtf or gff file.

  A bam file from alignment.
  
  A list containing three \link[GenomicRanges]{GRanges} objects generated from \code{translateGTF}.
  
  A single character representing sample name for the bam file.
  
  A single character representing sample ID for the bam file.

  A single character representing overlap strategy between reads and gene/exon/intron. 
}


\author{X. Sun and F. Zuo}

\seealso{
  \code{\link{translateGTF}},
  \code{\link{combineCount}}
}

\examples{

library(SplicingTypesAnno)
## setup
mm9.pnpla7.gtfFile <- system.file("extdata", "mm9_Pnpla7.gtf", package="SplicingTypesAnno")

## bam files are only for gene: Pnpla7. 
## The raw data are from GSM653002 and GSM653003 
bam.1 <- system.file("extdata", "liver_ctr.sort.bam", package="SplicingTypesAnno")
bam.2 <- system.file("extdata", "liver_ko.sort.bam", package="SplicingTypesAnno")

## translate gtf file

result.list <- translateGTF(mm9.pnpla7.gtfFile)

selectGene <- "Pnpla7"
sample1.c <- splicingCount(selectGene, bam.1, result.list,
                      sampleName="liver_ctr",
                      sampleID=1
                      )

sample1.c$totalCount
sample1.c$exonCount
sample1.c$intronCount

}

\keyword{methods}
\keyword{classes}
