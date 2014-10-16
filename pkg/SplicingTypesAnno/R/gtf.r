##################################
##  SplicingTypesAnno
##  - gtf3.r: translate gtf/gff to gene annotation
##################################

## note:
##  - gtf: http://mblab.wustl.edu/GTF22.html
##  - gff3: http://www.sequenceontology.org/gff3.shtml
##  - gff2: http://www.sanger.ac.uk/resources/software/gff/spec.html

# convert gtf df to GRange
gtf2GRange <- function(gtf.dataframe, geneName)
{
    if (!all(c("chr", "strand", "start", "end") %in% colnames(gtf.dataframe)))
         stop("You need chr, strand, start, end as column names of data.frame")
    gtfGRange <- GRanges(seqnames=Rle(gtf.dataframe$chr),
                          strand=gtf.dataframe$strand,
                          ranges=IRanges(start=gtf.dataframe$start, end=gtf.dataframe$end),
                          geneName=geneName
                          )
    return(gtfGRange)

}

## Main goal: 
## function 1: generate gene.GRange; 2: add exon number to each exon
##           3: delete exon from different transcripts with same gene name
##            4: delete exon+genes or some geomic component with only ONE exon
## output: gene.GRange, exon.GRange
## Note: control mutiple exon with same gene name
geneGRangeFromGTF <- function(exon.gtf, sameExonFromIsoform.collapse=TRUE, 
                              keepOneExonGene=FALSE, 
                              gtfGeneLabel="gene_id",
                              gtfGeneValue="g" # g.1
                              )
{
    ## ID
    exonName <- exon.gtf$geneName
    exonName.str <- unlist(strsplit(exonName, ";"))
         
    ID.ind <- grep(gtfGeneLabel, exonName.str)
    if (length(ID.ind)==0) stop("gtf file does NOT have gtfGeneLabel as gene identifiers!")
  
    gtfGeneLabel.str <- exonName.str[ID.ind]
    if (gtfGeneValue != "g")
    {
        if (gtfGeneValue == "g.1") ## for arabidposis version
        {
            dot.ind <- grep("\\.", gtfGeneLabel.str)

            tmp.str <- gtfGeneLabel.str[dot.ind]
            tmp.split <- unlist(strsplit(tmp.str, "\\."))
            tmp.done <- tmp.split[grep(gtfGeneLabel, tmp.split)]
            gtfGeneLabel.str[dot.ind] <- tmp.done

        }
    }

    exon.ID <- paste(gtfGeneLabel.str, exon.gtf$chr, exon.gtf$strand, sep="_")
    exon.gtf$ID <- exon.ID

    ## annotated multiple ID with same name
    tmp.df <- exon.gtf[, c("chr", "strand", "ID")]
    tmp.ind <- which(!duplicated(tmp.df)) 
    nodup.df <- exon.gtf[tmp.ind,]
    
    target.ind <- which(duplicated(nodup.df$ID))
    multiID <- nodup.df$ID[target.ind]
    
    exon.1 <- exon.gtf[!(exon.gtf$ID %in% multiID),]
    exon.2 <- exon.gtf[(exon.gtf$ID %in% multiID),]
    exon.2$geneName <- paste(exon.2$geneName, exon.2$chr, exon.2$strand, sep="_")
    exon.2$ID <- paste(exon.2$ID, exon.2$chr, exon.2$strand, sep="_")
    
    exon.gtf <- rbind(exon.1, exon.2)

    ## TODO: collapse same exon from different isoform?
    ## collapse those isoform with exactly same start, end, ID(str, strand, geneID)
    if (sameExonFromIsoform.collapse)
    {
      exon.dup.str <- paste(exon.gtf$ID, exon.gtf$start, exon.gtf$end, sep=".")
      exon.dup.ind <- which(duplicated(exon.dup.str))
      if (length(exon.dup.ind) > 0) exon.gtf <- exon.gtf[-exon.dup.ind,]
    }

    if (!keepOneExonGene)
    {
       table.id <- table(exon.gtf$ID)
       oneGene.ind <- which(table.id==1)
       oneGene.name <- names(table.id)[oneGene.ind]
       exon.gtf <- exon.gtf[ !(exon.gtf$ID %in% oneGene.name),]
    }
    
    firstPick <- function(tmp.vec)
    {
        unique(tmp.vec)[1]
    }
    
    gene.start <- tapply(exon.gtf$start, exon.gtf$ID, min)
    gene.end <- tapply(exon.gtf$end, exon.gtf$ID, max)
    gene.strand <- tapply(exon.gtf$strand, exon.gtf$ID, firstPick) 
    gene.seqname <- tapply(exon.gtf$chr, exon.gtf$ID, firstPick)
    gene.geneName <- tapply(exon.gtf$ID, exon.gtf$ID, firstPick)

    gene.df <- data.frame(chr=gene.seqname, strand=gene.strand, start=gene.start,
                          end=gene.end)
                          
                          
    gene.GRange <- gtf2GRange(gene.df, geneName=gene.geneName)

    return(list(gene.GRange=gene.GRange, exon.gtf=exon.gtf))
    
}

## Main goal: generate exon & intron number for each gene
GRangeNumberAnno <- function(intron.GRange, gene.GRange, string)
{                     
      tmp <- as.matrix(findOverlaps(intron.GRange,gene.GRange))
      tmp.nodu <- tmp[!duplicated(tmp[,1]), ,drop=FALSE]
                              
      geneName <- unlist( values(gene.GRange) )
      tmpName <- as.character(rep("No", length(intron.GRange)))
      
      #check <- try(as.character(geneName[tmp.nodu[,2]]))
      #if (inherits(check, "try-error")) browser()
      
      tmpName[tmp.nodu[,1]] <- as.character(geneName[tmp.nodu[,2]])
    
      ## calculate intron number
      name.rle <- rle(tmpName)
      
      name.len <- name.rle$lengths
      intron.no <- NULL
          
      sapply(1:length(name.len), function(i)
            {
                intron.no <<- c(intron.no, seq(1, name.len[i], 1))
                return(invisible(NULL))
            })
    
      elementMetadata(intron.GRange)[, "geneName"] <-  paste(tmpName, string, intron.no, sep="_")      
      elementMetadata(intron.GRange)[, "ID"] <-  tmpName    
        
      return(intron.GRange)
}


## Main goal: select exon with selected geneID
getGRangeWithName <- function(exon.GRange, select.GRange, selectID)
{
    if (!missing(select.GRange)) selectID <- elementMetadata(select.GRange)[, "ID"]
    exonID <- elementMetadata(exon.GRange)[, "ID"]
    exon.ind <- which(exonID %in% selectID)
    
    if (length(exon.ind)!=0) return(exon.GRange[exon.ind])
    else return(NULL)
    
}

## Main goal: 
produceExonIntronNo <- function(gene.GRange, exon.GRange, ExInonReduce=TRUE)
{

        if (ExInonReduce) exon.GRange <- unique(exon.GRange)
    
        exon.GRange <- sort(exon.GRange)      
      
          ## exon number ID
        exon.ID <- elementMetadata(exon.GRange)[, "ID"]
        exon.Name <- elementMetadata(exon.GRange)[, "geneName"]

          ## calculate exon number
          order.ind <- order(exon.ID)
          exon.ID <- exon.ID[order.ind]
          exon.Name <- exon.Name[order.ind]
          exon.GRange <- exon.GRange[order.ind]
          
          name.rle <- rle(exon.ID)
          name.len <- name.rle$lengths
          exon.no <- NULL
          sapply(1:length(name.len), function(i)
                {
                    exon.no <<- c(exon.no, seq(1, name.len[i], 1))
                    return(invisible(NULL))
                })      

          elementMetadata(exon.GRange)[, "geneName"] <- paste(exon.Name, "Exon_", exon.no, sep="")
         
          ## intron
          if (ExInonReduce) 
          {
              intron.GRange <- setdiff(gene.GRange, exon.GRange)
              ## annotate intron no
              intron.GRange <- GRangeNumberAnno(intron.GRange, gene.GRange, "intron")
          }
          else 
          {
              ## mutually exclusive group
              gene.reduce <- reduce(gene.GRange)
              gene.m <- as.matrix(findOverlaps(gene.reduce, gene.GRange))
              gene.ind <- which(gene.m[,1] - gene.m[,2] != 0)
              
              target.m <- data.frame(gene.m[gene.ind,])     
              
              target.rle <- rle(target.m[,1])
              target.len <- target.rle$lengths
              type.no <- NULL
              sapply(1:length(target.len), function(i)
                    {
                        type.no <<- c(type.no, seq(1, target.len[i], 1))
                        return(invisible(NULL))
                    })                
              
              target.m$type <- type.no
              type.unique <- sort(unique(type.no))

              tmp.target <- target.m[target.m$type==type.unique[1],]
              select.gene <- gene.GRange[tmp.target$subject]
              select.exon <- getGRangeWithName(exon.GRange, select.gene)                        
              select.intron <- setdiff(select.gene, select.exon)
              select.intron <- sort(select.intron)
              select.intron <- GRangeNumberAnno(select.intron, select.gene, "intron")
              intron.GRange <- select.intron
                                        
              sapply(2:length(type.unique), function(i)
                    {
                        tmp.target <- target.m[target.m$type==type.unique[i],]
                        select.gene <- gene.GRange[tmp.target$subject]
                        select.exon <- getGRangeWithName(exon.GRange, select.gene)                        
                        select.intron <- setdiff(select.gene, select.exon)
                        select.intron <- sort(select.intron)
                        if (length(select.intron)>0)
                        {   
                          select.intron <- GRangeNumberAnno(select.intron, select.gene, "intron")                  
                          intron.GRange <<- c(intron.GRange, select.intron)
                        }
                                               
                        invisible()
                    })
              
          }
                    
          return(list(exon.GRange=exon.GRange, intron.GRange=intron.GRange))
}

## generate transcript binding string for all known splicingLink
transcriptStr <- function(exon.gtf, gtfTranscriptLabel)
{          
    tmp.str <- unlist(strsplit(exon.gtf$geneName, ";"))
    transID.ind <- grep(gtfTranscriptLabel, tmp.str)
    transID <- tmp.str[transID.ind]
    if (length(transID) != nrow(exon.gtf) ) stop("Some exon does NOT have transcript ID!")
    
    exon.gtf$transID <- transID

    
    bindTrans <- function(start, end, chr, strand)
    {
        if (length(start)==1) return(NULL)
        
        str.start <- end[-length(end)]
        str.end <- start[-1]
        str.chr <- chr[-1]
        str.strand <- strand[-1]
        final.str <- paste(str.chr, ":", str.start, "-", str.end, "_", str.strand, sep="")
        
        return(final.str)
        
    }
    
    t.str <- tapply(exon.gtf$start, as.factor(exon.gtf$transID), bindTrans, exon.gtf$end, exon.gtf$chr, exon.gtf$strand)
    done.str <- unlist(t.str)
    return(done.str)
}


## Main function
## Main goal: generate gene/exon/intrin GRange object
translateGTF <- function(
                              gtfFile,  
                              gtfColnames=c("chr", "source", "feature", "start", "end","score", "strand", "frame", "geneName"),                               
                              gtfChrFormat="chr", # "Chr", "Chr0", "chr0", "1" 
                              gtfGeneLabel="gene_id",
                              gtfGeneValue= "g", # "g.1"
                              gtfTranscriptLabel=NULL, # "Parent="  "transcript_id"                          
                              exonString="exon",
                              selectGenes=NULL,
                              geneOverlap="both", # "yes", "no", "both"
                              exonOverlap="yes", # "yes", "no"
                              exonBoundary="no"
                              )
{    

  if (any(!(c("chr", "start", "end", "strand", "geneName") %in% gtfColnames))) stop("gtfColnames must have chr, strand, start, end and geneName column!")
  if (!(gtfChrFormat %in% c("chr", "Chr", "Chr0", "chr0", "1"))) stop("gtfChrFormat must be chr, Chr, chr0, Chr0, or 1!")
  if (!(gtfGeneValue %in% c("g", "g.1"))) stop("gtfGeneValue must be g or g.1!")
  if ((!is.null(selectGenes))  && (geneOverlap != "both")) 
  {
      print("SelecteGenes are used. geneOverlap has been set as no mode!")
      geneOverlap <- "no"
  }

  if (!(geneOverlap %in% c("both", "yes", "no"))) stop("geneOverlap can only be both, yes or no!")
  if (!(exonOverlap %in% c("yes", "no"))) stop("exonOverlap can only be yes or no!")

  r1.gtf <- read.delim(gtfFile, header=FALSE, skip=0, stringsAsFactors =FALSE)
  colnames(r1.gtf) <- gtfColnames
  
  ## transfer to 0-based for output bed format
  #r1.gtf$start <- r1.gtf$start - 1
  #r1.gtf$end <- r1.gtf$end - 1

  # check chr format
  if (gtfChrFormat!="chr")
  {
      if (gtfChrFormat=="Chr") r1.gtf$chr <- gsub("Chr", "chr", as.character(r1.gtf$chr))
      else if (gtfChrFormat=="Chr0") r1.gtf$chr <- gsub("Chr0", "chr", as.character(r1.gtf$chr))
      else if (gtfChrFormat=="chr0") r1.gtf$chr <- gsub("chr0", "chr", as.character(r1.gtf$chr))
      else if (gtfChrFormat=="1") r1.gtf$chr <- paste("chr", as.character(r1.gtf$chr), sep="")
      else stop("gtfChrFormat is NOT supported!")
  }
  
  exon.gtf <- r1.gtf[r1.gtf$feature==exonString,]
  
  ## convert selectGenes format to match gtf
  
  if (!is.null(selectGenes))
  {
      select.ind <- NULL
      sapply(1:length(selectGenes), function(i)
            {
                tmp.ind <- grep(selectGenes[i], exon.gtf$geneName) # check bug if
                if (length(tmp.ind)>0) select.ind <<- c(select.ind, tmp.ind)
                invisible(NULL)
            })
      if (length(select.ind)==0) stop("There are no genes in the selected gene list!")

      exon.gtf <- exon.gtf[select.ind,]
  }

  ## store exon start/end for further usage
  if (exonBoundary == "yes")
  {
      exon.start <- paste(exon.gtf$chr, ":", exon.gtf$start, "_", exon.gtf$strand, sep="")
      exon.end <- paste(exon.gtf$chr, ":", exon.gtf$end, "_", exon.gtf$strand, sep="")
  }
  else
  {
      exon.start <- NULL
      exon.end <- NULL
  }

  ## generate unique transcript link for detect novel splicing
  ## only work for selected genes to limit memory usage
  if (!is.null(gtfTranscriptLabel) && !is.null(selectGenes))
  {
      splicingLink <- NULL
      sapply(1:length(selectGenes), function(i)
            {
               tmp.ind <- grep(selectGenes[i], exon.gtf$geneName) # check bug if  
               tmp.str <- transcriptStr(exon.gtf[tmp.ind,], gtfTranscriptLabel)   
               splicingLink <<- c(splicingLink, tmp.str)
            })
  }
  else splicingLink <- NULL

  ## remove duplicate exon with same strand, same geneName (transcript same...)
  dup.ind <- which(duplicated(paste(exon.gtf$chr, exon.gtf$strand, 
                                    exon.gtf$start, exon.gtf$end, 
                                    exon.gtf$geneName)))
  if (length(dup.ind) > 0) exon.gtf <- exon.gtf[-dup.ind,]


  total.GRange <- geneGRangeFromGTF(exon.gtf, gtfGeneLabel=gtfGeneLabel, gtfGeneValue=gtfGeneValue)
  gene.GRange <- total.GRange$gene.GRange
  elementMetadata(gene.GRange)[, "ID"] <- elementMetadata(gene.GRange)[, "geneName"]

  exon.gtf <- total.GRange$exon.gtf
  exon.GRange <- gtf2GRange(exon.gtf, geneName=exon.gtf$geneName)
                                                
  elementMetadata(exon.GRange)[, "ID"] <- exon.gtf$ID      

        ## genes include: non-overlap/unique + overlapping 
        gene.c <- countOverlaps(gene.GRange, gene.GRange)

        if (geneOverlap=="no" || geneOverlap=="both")
        {
            g1.GRange <- gene.GRange[gene.c <= 1]
            if (length(g1.GRange)==0) 
            {
                print("No genes are unique in gtf!")
                gene.final.1 <- NULL
                exon.final.1 <- NULL
                intron.final.1 <- NULL
            }
            else
            {
              e1.GRange <- getGRangeWithName(exon.GRange, g1.GRange)
    
              tmp <- produceExonIntronNo(g1.GRange, e1.GRange, TRUE)
              exon.final.1 <- sort(tmp$exon.GRange)
              intron.final.1 <- sort(tmp$intron.GRange)
              gene.final.1 <- sort(g1.GRange)
              if (exonOverlap!="yes")
              {
                exon.reduce <- GRangeNumberAnno(reduce(exon.GRange), gene.GRange, "exon")
              }
            }
        }
        
        if (geneOverlap=="yes"   || geneOverlap=="both")
        {
            g1.GRange <- gene.GRange[gene.c > 1]
            if (length(g1.GRange)==0) 
            {
                #print("No genes are overlaped in gtf!")
                gene.final.2 <- NULL
                exon.final.2 <- NULL
                intron.final.2 <- NULL                
            }
            else
            {                           
              e1.GRange <- getGRangeWithName(exon.GRange, g1.GRange)
    
              tmp <- produceExonIntronNo(g1.GRange, e1.GRange, FALSE)
              exon.final.2 <- sort(tmp$exon.GRange )
              intron.final.2 <- sort(tmp$intron.GRange)
              gene.final.2 <- sort(g1.GRange)
            }
        }
        
        print("Done...")
        if (geneOverlap=="no") 
          return(list(gene=gene.final.1, intron=intron.final.1, exon=exon.final.1, exonStart=exon.start, exonEnd=exon.end, splicingLink=splicingLink))
        else if (geneOverlap=="yes") 
          return(list(gene=gene.final.2, intron=intron.final.2, exon=exon.final.2, exonStart=exon.start, exonEnd=exon.end, splicingLink=splicingLink))
        else
          return(list(gene=gene.final.1, intron=intron.final.1, exon=exon.final.1, exonStart=exon.start, exonEnd=exon.end, splicingLink=splicingLink,
                      gene2=gene.final.2, intron2=intron.final.2, exon2=exon.final.2
          ))
      
}


