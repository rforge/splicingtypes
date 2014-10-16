\name{combineCount}
\Rdversion{0.99}

\alias{combineCount}
\alias{combineCount,character-method}

\title{sample-level count calculation}

\description{
Combine read counts (gene, intron and exon) for all samples.
}

\usage{
combineCount(sample.list)
}

\arguments{
  \item{sample.list}{a list cotaining results from \code{\link{splicingCount}}}.
  
}


\value{
  A list.

}

\author{X. Sun and F. Zuo}

\examples{

library(SplicingTypesAnno)
## setup
mm9.pnpla7.gtfFile <- system.file("extdata", "mm9_Pnpla7.gtf", package="SplicingTypesAnno")
mm9.pnpla7.gtfFile

## bam files are only for gene: Pnpla7. 
## The raw data are from GSM653002 and GSM653003 
bam.1 <- file.path(system.file("extdata", "liver_ctr.sort.bam", package="SplicingTypesAnno"))
bam.2 <- file.path(system.file("extdata", "liver_ko.sort.bam", package="SplicingTypesAnno"))

## translate gtf file
result.list <- translateGTF(gtfFile=mm9.pnpla7.gtfFile)

selectGene <- "Pnpla7"
sample1.c <- splicingCount(selectGene, bam.1, result.list,
                      sampleName="liver_ctr",
                      sampleID=1
                      )

sample2.c <- splicingCount(selectGene, bam.2, result.list,
                      sampleName="liver_ko",
                      sampleID=2
                      )
sList.c <- list(sample1.c, sample2.c)
ccount <- combineCount(sList.c )

}

\keyword{methods}
\keyword{classes}
