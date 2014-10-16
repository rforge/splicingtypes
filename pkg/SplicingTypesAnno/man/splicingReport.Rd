\name{splicingReport}
\Rdversion{0.99}

\alias{splicingReport}
\alias{splicingReport,character-method}

\title{Gene-level and Sample-level annotation report of major alternative splicing types}

\description{
Generate web report for major alternative splicing types.
}

\usage{
splicingReport(inputData, 
                        gtfFile, 
                        gtfChrFormat="chr", 
                        gtfGeneLabel="gene_id", 
                        gtfGeneValue="g",
                        gtfTranscriptLabel=NULL, 
                        outputDir="html", 
                        eventType=c("ri.1", "ri.2", "es.1", "es.2", 
                                  "adleft.1", "adleft.2", 
                                  "adright.1", "adright.2",
                                  "adboth.1", "adboth.2"),
                        inputFile,
                        minReadCounts=10, 
                        selectGenes=NULL, 
                        novel="both", 
                        parallel=FALSE, 
                        cpus=2
                        )
}

\arguments{
  \item{inputData}{a list containing two required vectors: 
    \code{SampleNames} and \code{Bamfiles}. The third vector: SampleID is optional.

   }

  \item{gtfFile}{gtf or gff files.}

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

  \item{outputDir}{The folder name for web report. }
      
  \item{eventType}{Splicing types for analysis. It includes: \emph{"ri.1", "ri.2", "es.1", "es.2", 
       "adleft.1", "adleft.2", "adright.1", "adright.2","adboth.1", "adboth.2"}. 
       ri.1: type I of intron retention; ri.2: type II of intron retention;
       es.1: type I of exon skipping; es.2: type II of exon skipping;
       adleft.1: type I of alternative 5' site; adleft.2: type II of alternative 5' site;
       adright.1: type I of alternative 3' site; adright.2: type II of alternative 3' site;
       adboth.1: type I of alternative both sites; adboth.2: type II of alternative both sites. 
       See vignette for details.
  }

  \item{inputFile}{A csv file including two columns: 
      \emph{SampleName} and \emph{BamFiles}. The third column: SampleID is optional.}
  
  \item{minReadCounts}{The cutoff value for identifying the reads as splicing types.
    The default value is 10 reads.}
  \item{selectGenes}{The selected gene symbols. It takes two values: NULL or a character vector. If NULL is selected,
  all genes will be processed. The select gene symbols must match \code{gtfGeneValue} in the gtf or gff file.}  

  \item{novel}{A single character. This option only works for alternative 5' site, 
    alternative 3' site and alternative both sites. 
    It has three values: \emph{yes}, \emph{no}, \emph{both}.
      \emph{yes}: only annote the novel splicing types; 
      \emph{no}: only annote the known splicing types;
      \emph{both}: annote both novel and known splicing types.}
      
  \item{parallel}{Logical value for parallel computing. see \code{snowfall} for details.}
  \item{cpus}{The number of cpu for parallel computing. see \code{snowfall} for details.}  
  
}

\details{
  Input: gtf/gff file and bam file(s) for single or mutiple samples.
  
  Output: a web report.
}

\value{
  A list for the bam files.
  
  A character vector for gtf/gff file name.

  A character vector for chromosome format.
  
  A character vector for gene lables in the gtf/gff file.
  
  Null or a character vector for transcript ID.
  
  A character vector for output folder name.
  
  A chracter vector for major splicing types.
  
  A chracter vector for the input file, which should be csv format containing the information of all bam files.
  
  A single interger for minimum read number.
  
  A character vector for the selected genes.
  
  A character vector for the status of normalized ratio.
  
  A character vector for status of novel alternative sites.
  
  A character vector for status of parallel computing.
  
  A single interger for CPU numbers. 


}


\author{X. Sun and F. Zuo}

\seealso{
  \code{\link{splicingGene}},
  \code{\link{splicingShowReport}},
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

sampleList <- list(SampleName=c("liver_ctr", "liver_ko"),
                    BamFiles=c(bam.1, bam.2),
                    SampleID=c(1,2)
                     )
# NOT RUN
# splicingReport(inputData=sampleList, 
#               gtfFile=mm9.pnpla7.gtfFile, 
#               selectGenes=c("Pnpla7"))

# splicingReport(inputData=sampleList, 
#               gtfFile=mm9.pnpla7.gtfFile, 
#               selectGenes=c("Pnpla7"),
#               parallel=TRUE, cpus=2)

}

\keyword{methods}
\keyword{classes}
