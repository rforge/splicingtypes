\name{combineGene}
\Rdversion{0.99}

\alias{combineGene}
\alias{combineGene,character-method}

\title{Gene-level annotation of major alternative splicing types for all samples}

\description{
Combine all splicing annotations from \code{\link{splicingGene}}.
}

\usage{
combineGene(sampleList, splicingEvents=c("ri.1"), selectGene, all=TRUE)
}

\arguments{
  \item{sampleList}{a list of all samples annotation results from \code{\link{splicingGene}}. }
 
  \item{splicingEvents}{Splicing types for analysis. It includes: \emph{"ri.1", "ri.2", "es.1", "es.2", 
       "adleft.1", "adleft.2", "adright.1", "adright.2","adboth.1", "adboth.2"}. 
       ri.1: type I of intron retention; ri.2: type II of intron retention;
       es.1: type I of exon skipping; es.2: type II of exon skipping;
       adleft.1: type I of alternative 5' site; adleft.2: type II of alternative 5' site;
       adright.1: type I of alternative 3' site; adright.2: type II of alternative 3' site;
       adboth.1: type I of alternative both sites; adboth.2: type II of alternative both sites. 
       See vignette for details.
  }  
  \item{selectGene}{Gene name or gene symbol. This should match 
    \code{gtfGeneLabel} of \code{\link{translateGTF}} }

  \item{all}{logical value. See \code{all} of \code{merge} for details.}
}

\value{
  A list for storing results of all samples. 

  A character vector for splicing type.
  
  A character vector for gene symbol.
  
  A logical value for the approaches to merge data.

}


\author{X. Sun and F. Zuo}

\seealso{
  \code{\link{splicingGene}}
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
## B. Translate gtf/gff file to gene/exon/intron structure
## ---------------------------------------------------------------------
result.GRange <- translateGTF(mm9.pnpla7.gtfFile, gtfTranscriptLabel="transcript_id")

## ---------------------------------------------------------------------
## C. Annotate and quantify the major splicing types
## ---------------------------------------------------------------------

#selectGene <- c("Pnpla7")
#sample1 <- splicingGene(selectGene, bam.1, result.GRange,
#                        sampleName="liver_ctr", sampleID=1
#                        )

#sample2 <- splicingGene(selectGene, bam.2, result.GRange,
#                        sampleName="liver_KO", sampleID=2
#                        )

## ---------------------------------------------------------------------
## D. Combine results from all samples
## ---------------------------------------------------------------------
#sList <- list(sample1$type, sample2$type)
#ri.1 <- combineGene(sList, splicingEvents="ri.1", selectGene)
#ri.2 <- combineGene(sList, splicingEvents="ri.2", selectGene)
#es.1 <- combineGene(sList, splicingEvents="es.1", selectGene)
#es.2 <- combineGene(sList, splicingEvents="es.2", selectGene)
#adleft.1 <- combineGene(sList, splicingEvents="adleft.1", selectGene)
#adleft.2 <- combineGene(sList, splicingEvents="adleft.2", selectGene)
#adright.1 <- combineGene(sList, splicingEvents="adright.1", selectGene)
#adright.2 <- combineGene(sList, splicingEvents="adright.2", selectGene)
#adboth.1 <- combineGene(sList, splicingEvents="adboth.1", selectGene)
#adboth.2 <- combineGene(sList, splicingEvents="adboth.2", selectGene)

}

\keyword{methods}
\keyword{classes}
