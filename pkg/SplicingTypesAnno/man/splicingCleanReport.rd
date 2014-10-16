\name{splicingCleanReport}
\Rdversion{0.99}

\alias{splicingCleanReport}
\alias{splicingCleanReport,character-method}

\title{Clean the web report}

\description{
Delete the folder for web report.
}

\usage{
splicingCleanReport(outputDir="html")
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
  \code{\link{splicingShowReport}}
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

sampleList <- list(SampleName=c("liver_ctr", "liver_ko"),
                    BamFiles=c(bam.1, bam.2),
                    SampleID=c(1,2)
                     )
## NOT RUN
# splicingReport(inputData=sampleList, 
#               gtfFile=mm9.pnpla7.gtfFile, 
#               selectGenes=c("Pnpla7"))

# splicingReport(inputData=sampleList, 
#               gtfFile=mm9.pnpla7.gtfFile, 
#               selectGenes=c("Pnpla7"),
#               parallel=TRUE, cpus=2)

## Show report
# splicingShowReport()

## Delete report
# splicingCleanReport()

}

\keyword{methods}
\keyword{classes}
