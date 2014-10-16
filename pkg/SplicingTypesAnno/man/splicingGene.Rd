\name{splicingGene}
\Rdversion{0.99}

\alias{splicingGene}
\alias{splicingGene,character-method}

\title{Gene-level annotation of major alternative splicing types for single sample}

\description{
Extract reads for the selected gene(s) from bam file, and annotate major alternative splicing types.
}

\usage{
splicingGene(selectGene,
                      bamFile,
                      select.GRange,
                      eventType=c("ri.1", "ri.2", "es.1", "es.2", 
                                  "adleft.1", "adleft.2", 
                                  "adright.1", "adright.2",
                                  "adboth.1", "adboth.2"),
                      sampleName="sample",
                      sampleID=1, 
                      reportSummary=TRUE,
                      minReadCounts = 10, 
                      novel="both"
                      )
}

\arguments{
  \item{selectGene}{Gene name or gene symbol. This should match 
    \emph{gtfGeneLabel} of \code{\link{translateGTF}} }

  \item{bamFile}{Bam file for analysis, which is the alignment file from bowtie or tophat.}

  \item{select.GRange}{A list from \code{\link{translateGTF}}. This list includes 
    gene, exon, intron details }
  
  \item{eventType}{Splicing types for analysis. It includes: \emph{"ri.1", "ri.2", "es.1", "es.2", 
       "adleft.1", "adleft.2", "adright.1", "adright.2","adboth.1", "adboth.2"}. 
       ri.1: type I of intron retention; ri.2: type II of intron retention;
       es.1: type I of exon skipping; es.2: type II of exon skipping;
       adleft.1: type I of alternative 5' site; adleft.2: type II of alternative 5' site;
       adright.1: type I of alternative 3' site; adright.2: type II of alternative 3' site;
       adboth.1: type I of alternative both sites; adboth.2: type II of alternative both sites. 
       See vignette for details.
  }
  \item{sampleName}{Sample name for the bam file.}
  \item{sampleID}{Sample ID for the bam file.}
  \item{reportSummary}{todo}
  \item{minReadCounts}{The cutoff value for identifying the reads, either junction reads or nonjunection reads,
   as splicing types.
    The default value is 10 reads.The reads which do not have
    enough copies will be dicarded, and they are not used for 
    analysis and annotation.}

  \item{novel}{A single character. This option only works for alternative 5' site, 
    alternative 3' site and alternative both sites. 
    It has three values: \emph{yes}, \emph{no}, \emph{both}.
      \emph{yes}: only annote the novel splicing types; 
      \emph{no}: only annote the known splicing types;
      \emph{both}: annote both novel and known splicing types.}
  
}

\details{
  See vignette for the definition of splicing types.
  
  Output includes:
  \emph{mergeID}: the ID for splicing types.
  \emph{locus}: the exact coordinate for this splicing type. For intron retention,
  it should be the start and the end of the target intron; For exon skipping, it should be
  the start and the end of the skipped exon; For alternative sites, it should be the
  end of the left junction read and the start of the right junection read.
}

\value{
  A character vector representing the selected gene.

  A character vector representing the bam file name.
  
  A GRanges object including gene, exon, and intron structure. It is generated from
  \code{\link{translateGTF}}.
  
  A character vector for major splicing types.
  
  A character vector for sample name.
  
  A single integer for sample ID.
  
  A logical value for report status.
  
  A single interger for minimum read numbers.
  
  A logical value for normalized ratio.
  
  A character vector for status of novel alternative sites.  

}

\author{Xiaoyong Sun}

\seealso{
  \code{\link{translateGTF}},
  \code{\link{splicingReport}}
}

\examples{
## ---------------------------------------------------------------------
## A. splicingGene()
## ---------------------------------------------------------------------
## data
# mouse
# http://www.ncbi.nlm.nih.gov/sra/SRX000350?report=full
# http://www.ncbi.nlm.nih.gov/sra/SRX000351
# Mapping and quantifying mammalian transcriptomes by RNA-Seq
## data 2
## http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26561
# Transcriptome-wide study of alternative splicing in nonsense-mediated mRNA decay-deficient mouse tissues
# http://www.ncbi.nlm.nih.gov/sra?term=SRP005384

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

# type <- sample2$type
## intron retention
#type$ri.1
#type$ri.2

## exon skipping
#type$es.1

## alternative donor site
#type$adleft.1


}

\keyword{methods}
\keyword{classes}
