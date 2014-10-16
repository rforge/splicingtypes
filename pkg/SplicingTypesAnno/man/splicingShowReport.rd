\name{splicingShowReport}
\Rdversion{0.99}

\alias{splicingShowReport}
\alias{splicingShowReport,character-method}

\title{Browser the web report for alternative splicing types}

\description{
Browser the web report for alternative splicing types.
}

\usage{
splicingShowReport(outputDir="html")
}

\arguments{
  \item{outputDir}{the web-report folder name.}
}


\value{
  A single character for the folder name.
}


\author{X. Sun}

\seealso{
  \code{\link{splicingReport}},
  \code{\link{splicingCleanReport}}
}

\examples{
library(SplicingTypesAnno)

## ---------------------------------------------------------------------
## A. Setup path
## ---------------------------------------------------------------------
mm9.pnpla7.gtfFile <- system.file("extdata", "mm9_Pnpla7.gtf", package="SplicingTypesAnno")

## bam files from GSM653002 and GSM653003 
bam.1 <- system.file("extdata", "liver_ctr.sort.bam", package="SplicingTypesAnno")
bam.2 <- system.file("extdata", "liver_ko.sort.bam", package="SplicingTypesAnno")

## ---------------------------------------------------------------------
## B. Generate report
## ---------------------------------------------------------------------

#sampleList <- list(SampleName=c("liver_ctr", "liver_ko"),
#                    BamFiles=c(bam.1, bam.2),
#                    SampleID=c(1,2))
#splicingReport(inputData=sampleList, 
#               gtfFile=mm9.pnpla7.gtfFile, 
#               selectGenes=c("Pnpla7"))                      
## NOT RUN
## Show report
# splicingShowReport()
}

\keyword{methods}
\keyword{classes}
