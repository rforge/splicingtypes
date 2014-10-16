##################################
##  SplicingTypesAnno
##  - algorithm.r: heuristic for detecting splicing types 
##################################

## TODO
## 0. summary statistics: all junction read, distribution...
# 1. all junction read NOT match annotation!!
# 2. all junction read has inside exon
## 3. all junction read has ..


## final output format for RI
riMessage <- function(intron.final, Type, sampleName, sampleID)
{               
            if (!is.null(intron.final))
            {
              print(paste("Intron Retention...", Type, "...FOUND...", sep=""))
              
              if (Type=="Type I")
                colnames(intron.final) <- c("mergeID", "geneName", "chr", "locus",
                                        paste(sampleName, sampleID, colnames(intron.final)[-c(1:4)], sep=".")
                                        #, paste(sampleName, sampleID, "score", sep=".")                             
                                        )
              else
                colnames(intron.final) <- c("mergeID", "geneName", "chr", "locus",
                                        paste(sampleName, sampleID, colnames(intron.final)[-c(1:4)], sep=".")             
                                        #, paste(sampleName, sampleID, "score", sep=".")                             
                                        )              
              return(intron.final)      
            }
            else
            {
                print(paste("Intron Retention...", Type, "...none...", sep=""))
                return(invisible(NULL))
            }
}

## final output format for ES
esMessage <- function(es.final, Type, sampleName, sampleID)
{
          if (!is.null(es.final))
          {
              print(paste("Exon skipping...", Type, "...FOUND...", sep=""))
              colnames(es.final) <- c("mergeID", "geneName", "chr", "locus",
                                  paste(sampleName, sampleID, colnames(es.final)[-c(1:4)], sep=".")
                                  #, paste(sampleName, sampleID, "score", sep=".")      
                                  )
              return(es.final)
          }
          else
          {
              print(paste("Exon skipping...", Type, "...none...", sep=""))
              return(invisible(NULL))
          }
}

## annotate novel boundary for ri and es
## fstart, fend are exon boundary
annotateNovel <- function(final, etype, novel, exonStart=NULL, exonEnd=NULL, 
                          splicingLink=NULL)
{

  if (!is.null(final) && nrow(final) > 0)
  {
        mergeID <- final$mergeID
        str.1 <- unlist(strsplit(mergeID, ":"))
        chr <- str.1[seq(1, length(str.1)-1,2)]
        str.1.left <- str.1[seq(2, length(str.1),2)]
        
        str.2 <- unlist(strsplit(str.1.left , "_"))
        strand <- str.2[seq(2, length(str.2), 2)]
        str.2.left <- str.2[seq(1, length(str.2)-1, 2)]
        
        str.3 <- unlist(strsplit(str.2.left, "-"))
        start <- str.3[seq(1, length(str.3)-1,2)]
        end <- str.3[seq(2, length(str.3),2)]
        
       if (etype=="ri.1" || etype=="ri.2")
       {
          start <- as.numeric(start)-1
          end <- as.numeric(end)+1        
       }
       
       start.str <- paste(chr, ":", start, "_", strand, sep="")
       end.str <- paste(chr, ":", end, "_", strand, sep="")     
       link.str <- paste(chr, ":", start, "-", end, "_", strand, sep="")
       
       ## annotate
       novel.ind <- NULL
       if (!is.null(exonStart))
       {                     
           novel.start <- rep(0, nrow(final))
           
           ## switch start match exon end because all events use it except Exon Skipping       
           start.ind <- which(!(start.str %in% exonEnd))
           novel.ind <- c(novel.ind, start.ind)
           
           if (length(start.ind)>0)
           {
              novel.start[start.ind] <- start[start.ind]
           }
           final$novelStart <- novel.start
       }
 
       if (!is.null(exonEnd))
       { 
           novel.end <- rep(0, nrow(final))   
           end.ind <- which(!(end.str %in% exonStart)) 
           novel.ind <- c(novel.ind, end.ind)
                  
           if (length(end.ind)>0)
           {
               novel.end[end.ind] <- end[end.ind]
           }
           final$novelEnd <- novel.end
       }
       
       if (!is.null(splicingLink))
       {    
          novel.link <- rep("0", nrow(final))
          link.ind <- which(!(link.str %in% splicingLink))
          novel.ind <- c(novel.ind, link.ind)
          if(length(link.ind)>0)
          {
              novel.link[link.ind] <- link.str[link.ind] 
          }
          final$novelSplicingLink <- novel.link
       }
      
        novel.ind <- unique(novel.ind)      
       if (novel=="yes")  
       {
          if (length(novel.ind) > 0) final <- final[novel.ind,]
          else return(NULL)
       }
       else if (novel=="no")
       {
           
           if (length(novel.ind) > 0) final <- final[-novel.ind,]
           else return(final)
       }
           
   }  
   else return(NULL)
       
   if (nrow(final)>0) return(final)
   else return(NULL)
}

## Swap donor/acceptor based on strand orientation
swapStrand <- function(adleftRaw.1, adrightRaw.1)
{
      if (is.null(adleftRaw.1) && is.null(adrightRaw.1)) return(list(adleft=NULL, adright=NULL))
       adleft.1 <- NULL
       adright.1 <- NULL
       adleft.a <- NULL 
       adleft.b <- NULL
       adright.a <- NULL 
       adright.b <- NULL
      
       if (!is.null(adleftRaw.1))
       { 
         istrandLeft.str <- unlist(strsplit(adleftRaw.1$mergeID, "_"))
         istrandLeft <- istrandLeft.str[seq(2, length(istrandLeft.str), 2)]
         if (length(istrandLeft) != nrow(adleftRaw.1)) stop("something wrong with strand orientation!")
         adleft.a <- adleftRaw.1[istrandLeft == "+",]
         adright.a <- adleftRaw.1[istrandLeft == "-",]
       }
       
       if (!is.null(adrightRaw.1))
       {
         istrandRight.str <- unlist(strsplit(adrightRaw.1$mergeID, "_"))
         istrandRight <- istrandRight.str[seq(2, length(istrandRight.str), 2)]
         if (length(istrandRight) != nrow(adrightRaw.1)) stop("something wrong with strand orientation!")
         adleft.b <- adrightRaw.1[istrandRight=="-",]
         adright.b <- adrightRaw.1[istrandRight=="+",]
       }
     
       adleft <- rbind(adleft.a, adleft.b)
       adright <- rbind(adright.a, adright.b)
       return(list(adleft=adleft, adright=adright))
  
}
  
## Main function for sample-level annotation
updateTRUEIntron <- function(
                      bamFile = "accepted_hits.sort.bam",
                      gtf.GRange,
                      sampleName,
                      sampleID,
                      eventType=c("ri.1", "ri.2", "es.1", "es.2", 
                                  "adleft.1", "adleft.2", 
                                  "adright.1", "adright.2",
                                  "adboth.1", "adboth.2"),
                      minReadCounts = 10,
                      #junctionRead.cutoff=2,
                      #nonJunRead.cutoff=2,
                      selectGenes=NULL, 
                      novel="yes"
                      )
{       
  result.list <- list()
  readSum <- list(totalRead=0,
           qualityJunRead=0,
           qualityNonJunRead=0,
           sampleLevel=NULL)
             
  if (!is.null(selectGenes))
  {
      select.GRange <- gtf.GRange
                            
      for (i in 1:length(selectGenes))
      {

          type.list <- splicingGene(selectGenes[i], bamFile, select.GRange, 
                       eventType=eventType,
                       sampleName=sampleName,
                      sampleID=sampleID,
                      reportSummary=TRUE,
                      #junctionRead.cutoff=junctionRead.cutoff,
                      #nonJunRead.cutoff=nonJunRead.cutoff,
                      minReadCounts = minReadCounts ,
                      novel=novel 
                       )
          
          tmp.list <- type.list$type
          ## TODO: select gene list?
          tmp.readSum <- type.list$summary
          readSum$totalRead <- readSum$totalRead + tmp.readSum$totalRead
          readSum$qualityJunRead<- readSum$qualityJunRead + tmp.readSum$qualityJunRead
          readSum$qualityNonJunRead<- readSum$qualityNonJunRead + tmp.readSum$qualityNonJunRead

          ## combine all
          result.list[["ri.1"]] <- rbind(result.list[["ri.1"]], tmp.list[["ri.1"]])
          result.list[["ri.2"]] <- rbind(result.list[["ri.2"]], tmp.list[["ri.2"]])                      
          result.list[["es.1"]] <- rbind(result.list[["es.1"]], tmp.list[["es.1"]])
          result.list[["es.2"]] <- rbind(result.list[["es.2"]], tmp.list[["es.2"]])                      
          result.list[["adleft.1"]] <- rbind(result.list[["adleft.1"]], tmp.list[["adleft.1"]])
          result.list[["adleft.2"]] <- rbind(result.list[["adleft.2"]], tmp.list[["adleft.2"]])                      
          result.list[["adright.1"]] <- rbind(result.list[["adright.1"]], tmp.list[["adright.1"]])
          result.list[["adright.2"]] <- rbind(result.list[["adright.2"]], tmp.list[["adright.2"]])                      
          result.list[["adboth.1"]] <- rbind(result.list[["adboth.1"]], tmp.list[["adboth.1"]])
          result.list[["adboth.2"]] <- rbind(result.list[["adboth.2"]], tmp.list[["adboth.2"]])                      
          
      }
      return(list(type=result.list, summary=readSum))
  } 

  ## check one gene range
  mydata <- readGAlignments(bamFile)
  readSum$sampleLevel <- TRUE


  if (any(!(seqlevels(mydata) %in% seqlevels(gtf.GRange$gene.GRange)) ))
  {
      print("some chrNo in bam file are NOT in gtf/gff3 file! The unmatched data will be overlooked...")
      #seqlevels(mydata, force=TRUE) <- seqlevels(splicingENV$gene.GRange)
  }
  if (any(!(seqlevels(gtf.GRange$gene.GRange) %in% seqlevels(mydata)) ))
  {
      print("some chrNo in gtf/gff3 file are NOT in bam file! The unmatched annotation will be overlooked...")
      seqlevels(gtf.GRange$gene.GRange, force=TRUE) <- seqlevels(mydata)
      seqlevels(gtf.GRange$exonGRange, force=TRUE) <- seqlevels(mydata)
      seqlevels(gtf.GRange$intronGRange, force=TRUE) <- seqlevels(mydata)
      #assign("seqlevels(gene.GRange, force=TRUE)", seqlevels(mydata), splicingENV)
  }

  ## check intron percentage of 5%
  #count.c <- countOverlaps(splicingENV$intronGRange, mydata)
  #count.sort <- sort(count.c)  
  #splicingENV$intronCutoffPerm <- c(splicingENV$intronCutoffPerm, count.sort[round((1-0.05)*length(count.sort))])

  ## gene summary
   readSum$totalRead <- length(mydata)

  ## PROCESS START 
  checkedData <- dataCheck(mydata, gtf.GRange$gene.GRange)
  if (is.null(checkedData)) 
  {
     readSum$qualityJunRead <- 0
     readSum$qualityNonJunRead <- 0
      
      return(list(type=NULL, summary=readSum))
  }
  left.jun <- checkedData$left.jun
  right.jun <- checkedData$right.jun
  middle.reduce <- checkedData$middle.reduce
  middle.jun.0 <- checkedData$middle.read
  nonjun.read <- checkedData$nonjun.read
  cis.read <- checkedData$cis.read
  
  ## read summary
  readSum$qualityJunRead <- length(cis.read)  
  readSum$qualityNonJunRead <- length(nonjun.read)
    
  ##########################################################################
  ## 1. intron retention
  ##########################################################################
  ## 1. check middle.jun has nonjun.read inside
  if ("ri.1" %in% eventType || "ri.2" %in% eventType)
  {
        
        if ("ri.1" %in% eventType) 
        {           
            intron.list <- riEvent(middle.reduce, middle.jun.0, nonjun.read, 
                                minReadCounts, 
                                eventType="ri.1",
                                sampleName, sampleID, total.read=readSum$totalRead, gtf.GRange)
            intron.final <- intron.list$type1
            intron.final <- annotateNovel(intron.final, "ri.1", novel, 
                                          gtf.GRange$exonStart, 
                                          gtf.GRange$exonEnd, 
                                          gtf.GRange$splicingLink)

            result.list[["ri.1"]] <- riMessage(intron.final, "Type I", sampleName, sampleID)
            
        }
       
        if ("ri.2" %in% eventType)
        {
            intron.list <- riEvent(middle.reduce, middle.jun.0, nonjun.read, 
                                minReadCounts, 
                                eventType="ri.2",
                                sampleName, sampleID, total.read=readSum$totalRead, gtf.GRange)       
            intron.final <- intron.list$type2
            intron.final <- annotateNovel(intron.final, "ri.2", novel, 
                                          gtf.GRange$exonStart, 
                                          gtf.GRange$exonEnd, 
                                          gtf.GRange$splicingLink)        
            result.list[["ri.2"]] <- riMessage(intron.final, "Type II", sampleName, sampleID)        
        }
    }    
          
  ##########################################################################
  ## 2. exon skipping
  ##########################################################################
  ## 1. use middle.jun.0
  #exonGRange.reduce <- reduce(exonGRange)
  #es.c1 <- countOverlaps(exonGRange.reduce, middle.reduce, type="within")
  
  ## 2. classic: jun read contain only ONE exonreduce
  #exon.1 <- exonGRange[es.c1==1]
  
  ## TODO: include those have more than one exons
  if ("es.1" %in% eventType || "es.2" %in% eventType)
  {

        if ("es.1" %in% eventType) 
        {
            es.list <- esEvent(middle.reduce, middle.jun.0, nonjun.read, 
                         minReadCounts, eventType="es.1", 
                        sampleName, sampleID, 
                         total.read=readSum$totalRead,
                         gtf.GRange)                              
        
            es.final <- es.list$type1
            ## note: for novel, exon start, end switch
            es.final <- annotateNovel(es.final, "es.1", novel, 
                                          gtf.GRange$exonEnd, 
                                          gtf.GRange$exonStart, 
                                          NULL)
            result.list[["es.1"]] <- esMessage(es.final, "Type I", sampleName, sampleID)
            
        }
        
        if ("es.2" %in% eventType)
        {
            es.list <- esEvent(middle.reduce, middle.jun.0, nonjun.read, 
                         minReadCounts, eventType="es.2", 
                         sampleName, sampleID, 
                         total.read=readSum$totalRead, gtf.GRange)          
            es.final <- es.list$type2
            es.final <- annotateNovel(es.final, "es.2", novel, 
                                          gtf.GRange$exonEnd, 
                                          gtf.GRange$exonStart, 
                                          NULL)
            result.list[["es.2"]] <- esMessage(es.final, "Type II", sampleName, sampleID)        
        }
  }
              
  ##########################################################################
  ## 3. alternative donor sites - no reference gtf
  ## - note: need to switch adleft, adright because of strand orientation
  ##########################################################################
  adleftRaw.1 <- asEvent(cis.read, left.jun, right.jun, 
                  sampleName, sampleID, minReadCounts, eventType="adleft.1", novel, 
                  total.read=readSum$totalRead, gtf.GRange)
  
  adleftRaw.2 <- asEvent(cis.read, left.jun, right.jun, 
                  sampleName, sampleID, minReadCounts, eventType="adleft.2", novel, 
                  total.read=readSum$totalRead, gtf.GRange)
 
  adrightRaw.1 <- asEvent(cis.read, left.jun, right.jun, 
                  sampleName, sampleID, minReadCounts, eventType="adright.1", novel, 
                  total.read=readSum$totalRead, gtf.GRange)
  
  adrightRaw.2 <- asEvent(cis.read, left.jun, right.jun, 
                    sampleName, sampleID, minReadCounts, eventType="adright.2", novel, 
                    total.read=readSum$totalRead, gtf.GRange)

  list.1 <- swapStrand(adleftRaw.1, adrightRaw.1)
  adleft.1 <- list.1$adleft
  adright.1 <- list.1$adright
  
  adleft.1 <-  annotateNovel(adleft.1, "adleft.1", novel, 
                             gtf.GRange$exonStart, 
                             gtf.GRange$exonEnd, 
                             gtf.GRange$splicingLink)
                             
  adright.1 <-  annotateNovel(adright.1, "adright.1", novel, 
                             gtf.GRange$exonStart, 
                             gtf.GRange$exonEnd, 
                             gtf.GRange$splicingLink)
  
  list.2 <- swapStrand(adleftRaw.2, adrightRaw.2)
  adleft.2 <- list.2$adleft
  adright.2 <- list.2$adright 

  adleft.2 <-  annotateNovel(adleft.2, "adleft.2",  novel, 
                             gtf.GRange$exonStart, 
                             gtf.GRange$exonEnd, 
                             gtf.GRange$splicingLink)
                             
  adright.2 <-  annotateNovel(adright.2, "adright.2", novel, 
                             gtf.GRange$exonStart, 
                             gtf.GRange$exonEnd, 
                             gtf.GRange$splicingLink)
  
  if (is.null(adleft.1)) print("Alternative Donor...Type I...none...")
  else print("Alternative Donor...Type I...FOUND...")
  result.list[["adleft.1"]] <- adleft.1

  if (is.null(adleft.2)) print("Alternative Donor...Type II...none...")
  else print("Alternative Donor...Type II...FOUND...")  
  result.list[["adleft.2"]] <- adleft.2

  if (is.null(adright.1)) print("Alternative Acceptor...Type I...none...")
  else print("Alternative Acceptor...Type I...FOUND...")     
  result.list[["adright.1"]] <- adright.1

  if (is.null(adright.2)) print("Alternative Acceptor...Type II...none...")
  else print("Alternative Acceptor...Type II...FOUND...")   
  result.list[["adright.2"]] <- adright.2
  
  ## adboth
  adboth.1 <- asEvent(cis.read, left.jun, right.jun,  
                  sampleName, sampleID, minReadCounts, eventType="adboth.1", novel, 
                  total.read=readSum$totalRead, gtf.GRange)
                  
  adboth.1 <-  annotateNovel(adboth.1, "adboth.1", novel, 
                             gtf.GRange$exonStart, 
                             gtf.GRange$exonEnd, 
                             gtf.GRange$splicingLink)
                             
  if (is.null(adboth.1)) print("Alternative Both Sites...Type I...none...")
  else print("Alternative Both Sites...Type I...FOUND...")
  result.list[["adboth.1"]] <- adboth.1
  
  adboth.2 <- asEvent(cis.read, left.jun, right.jun, 
                  sampleName, sampleID, minReadCounts, eventType="adboth.2", novel, 
                  total.read=readSum$totalRead, gtf.GRange)
  adboth.2 <-  annotateNovel(adboth.2, "adboth.2", novel, 
                             gtf.GRange$exonStart, 
                             gtf.GRange$exonEnd, 
                             gtf.GRange$splicingLink)
  if (is.null(adboth.2)) print("Alternative Both Sites...Type II...none...")
  else print("Alternative Both Sites...Type II...FOUND...")
  result.list[["adboth.2"]] <- adboth.2
  
  return(list(type=result.list, summary=readSum))

}

## Main function for gene-level annotation
splicingGene <- function(selectGene,
                      bamFile,
                      select.GRange,
                      eventType=c("ri.1", "ri.2", "es.1", "es.2", 
                                  "adleft.1", "adleft.2", 
                                  "adright.1", "adright.2",
                                  "adboth.1", "adboth.2"),
                      sampleName="sample",
                      sampleID=1, 
                      reportSummary=TRUE,
                      #junctionRead.cutoff=2,
                      #nonJunRead.cutoff=2,
                      minReadCounts = 10, novel="both"
                      )
{
  if (missing(selectGene)) stop("selectGenes Must be chosen!") 
  if (length(selectGene) != 1) stop("Please select only one gene each time!")
  if (any(!(eventType %in% c("ri.1", "ri.2", "es.1", "es.2", 
                                  "adleft.1", "adleft.2", 
                                  "adright.1", "adright.2",
                                  "adboth.1", "adboth.2"))))
      {
          print("eventType does not have this type!")
          return(invisible(NULL))
      }

  if (minReadCounts < 2) stop("Only minReadCounts > =2 is accepted!")

  exonGRange <- select.GRange$exon
  gene.GRange <- select.GRange$gene
  intronGRange <- select.GRange$intron
  exonStart <- select.GRange$exonStart
  exonEnd <- select.GRange$exonEnd
  splicingLink <- select.GRange$splicingLink
  
  #if (!is.null(select.GRange$gene2)) gene.GRange <- c(gene.GRange, select.GRange$gene2)
    
 
  if (length(gene.GRange)==0 || length(exonGRange)==0) stop("select.GRange has 0 gene or exon!")
  if (length(intronGRange)==0) print("select.GRane has 0 intron!")

  ## TODO: delete is.null(), overlap with length()==1 above

      ## get gene range

     geneName <- values(gene.GRange)[, "geneName"]
     gene.ind <- grep(selectGene, geneName)
     if (length(gene.ind)==0) 
     {
        print(paste("There is no gene named as ", selectGene, sep=""))
        return(invisible(NULL))
     }
               
     geneName <- values(exonGRange)[, "geneName"]
     exon.ind <- grep(selectGene, geneName)
     if (length(exon.ind)==0) 
     {
        print(paste("There is no exon named as ", selectGene, sep=""))
        return(invisible(NULL))
     }               
     #geneName <- values(intronGRange)[, "geneName"]
     #intron.ind <- grep(selectGene, geneName)             
            
      gene.GRange <- gene.GRange[gene.ind]
      exonGRange <- exonGRange[exon.ind]
      intronGRange <- setdiff(gene.GRange, reduce(exonGRange))
      
      ## intron annotation 
      #if (length(exonGRange) > 0) elementMetadata(exonGRange)[, "geneName"] <- paste(elementMetadata(exonGRange)[, "geneName"], "exon", 1:length(exonGRange), sep="_")      
      if (length(intronGRange) > 0) 
      {
        elementMetadata(intronGRange)[, "geneName"] <- paste(elementMetadata(gene.GRange)[, "geneName"], "intron", 1:length(intronGRange), sep="_")
        elementMetadata(intronGRange)[, "ID"] <- paste(elementMetadata(gene.GRange)[, "ID"], "intron", 1:length(intronGRange), sep="_")      
      }
  

  result.list <- list()
  result.list[["exon"]] <- exonGRange
  result.list[["intron"]] <- intronGRange
        
      print(paste("GENE - ", selectGene, " are processed...", sep=""))
      mydata <- selectGappedAlignments(bamFile=as.character(bamFile),
                         selectGene, gene.GRange
                         )

      if (length(mydata)==0)  
      {
          print("There is no read data with the select gene!")
          return(NULL)
      }
      
      if (any(!(seqlevels(mydata) %in% seqlevels(gene.GRange))))
      {
          print("some chrNo in bam file are NOT in gtf/gff3 file! The unmatched data will be overlooked...")
          #seqlevels(mydata, force=TRUE) <- seqlevels(gene.GRange)
      }
      if (any(!(seqlevels(gene.GRange) %in% seqlevels(mydata)) ))
      {
          print("some chrNo in gtf/gff3 file are NOT in bam file! The unmatched annotation will be overlooked...")
          seqlevels(gene.GRange, force=TRUE) <- seqlevels(mydata)
          seqlevels(exonGRange, force=TRUE) <- seqlevels(mydata)
          seqlevels(intronGRange, force=TRUE) <- seqlevels(mydata)
      }


      checkedData <- dataCheck(mydata, gene.GRange)
      if (is.null(checkedData)) 
      {
         #print("No junction reads have boundary overlapped! Either read does NOT match gtf file, or the number of reads are NOT enough!")
         return(invisible(NULL))
      }
      left.jun <- checkedData$left.jun
      right.jun <- checkedData$right.jun
      middle.reduce <- checkedData$middle.reduce
      middle.jun.0 <- checkedData$middle.read
      nonjun.read <- checkedData$nonjun.read
      cis.read <- checkedData$cis.read

      readSum <- list(totalRead=0,
                 qualityJunRead=0,
                 qualityNonJunRead=0,
                 sampleLevel=NULL)
                 
      if (reportSummary)
      {                 
        readSum$totalRead <- length(mydata)
        readSum$qualityJunRead  <-  length(cis.read)  
        readSum$qualityNonJunRead  <- length(nonjun.read)  
      }
      
      select.GRange <- list(gene.GRange=gene.GRange, 
                            exonGRange=exonGRange, 
                            intronGRange=intronGRange,
                            splicingLink=splicingLink
                            )
      if ("ri.1" %in% eventType || "ri.2" %in% eventType)
      {         
        if ("ri.1" %in% eventType)
        {                                          
            intron.list <- riEvent(middle.reduce, middle.jun.0, nonjun.read,
                                minReadCounts, 
                                eventType="ri.1",
                                sampleName, sampleID, readSum$totalRead, select.GRange)
            intron.final <- intron.list$type1
            intron.final <- annotateNovel(intron.final, "ri.1", novel, 
                                exonStart, exonEnd, splicingLink)

            result.list[["ri.1"]] <- riMessage(intron.final, "Type I", sampleName, sampleID)

        }
        
        if ("ri.2" %in% eventType)
        {      
        
           intron.list <- riEvent(middle.reduce, middle.jun.0, nonjun.read,
                                minReadCounts,
                                eventType="ri.2", 
                                sampleName, sampleID, readSum$totalRead, select.GRange)                        
            intron.final <- intron.list$type2
            intron.final <- annotateNovel(intron.final, "ri.2", novel, 
                                exonStart, exonEnd, splicingLink)     
            result.list[["ri.2"]] <- riMessage(intron.final, "Type II", sampleName, sampleID)
        }
      }
               
      if ("es.1" %in% eventType || "es.2" %in% eventType)
      {        

        if ("es.1" %in% eventType) 
        {         
            es.list <- esEvent(middle.reduce, middle.jun.0, nonjun.read, 
                          minReadCounts,
                          eventType="es.1", 
                          sampleName, sampleID, readSum$totalRead, select.GRange)                              

            es.final <- es.list$type1
            es.final <- annotateNovel(es.final, "es.1", novel, 
                              exonEnd, exonStart, NULL)            
            result.list[["es.1"]] <- esMessage(es.final, "Type I", sampleName, sampleID)
            
        }
        
        if ("es.2" %in% eventType)
        {
            es.list <- esEvent(middle.reduce, middle.jun.0, nonjun.read, 
                          minReadCounts,
                          eventType="es.2", 
                          sampleName, sampleID, readSum$totalRead, select.GRange)                              
 
            es.final <- es.list$type2
            es.final <- annotateNovel(es.final, "es.2", novel, 
                          exonEnd, exonStart, NULL)
            result.list[["es.2"]] <- esMessage(es.final, "Type II", sampleName, sampleID)        
        }
      }
     
     ## swap strand orientation
     if ("adleft.1" %in% eventType || "adright.1" %in% eventType)
      {      
            adleftRaw.1 <- asEvent(cis.read, left.jun, right.jun,
                                              sampleName, sampleID, minReadCounts,
                                              eventType="adleft.1", novel, readSum$totalRead, 
                                              select.GRange
                                              )
                                              
            adrightRaw.1 <- asEvent(cis.read, left.jun, right.jun,
                sampleName, sampleID, 
                minReadCounts,
                eventType="adright.1" , novel, readSum$totalRead, 
                select.GRange
                )

            list.1 <- swapStrand(adleftRaw.1, adrightRaw.1)
            adleft.1 <- list.1$adleft
            adright.1 <- list.1$adright 
    
            ## annotate novel
            adleft.1 <- annotateNovel(adleft.1, "adleft.1", novel, exonStart, exonEnd, splicingLink)  
            adright.1 <- annotateNovel(adright.1, "adright.1", novel, exonStart, exonEnd, splicingLink)  
                
            result.list[["adleft.1"]] <- adleft.1
            if (is.null(adleft.1)) print("Alternative Donor...Type I...none...")
            else print("Alternative Donor...Type I...FOUND...")

            result.list[["adright.1"]] <- adright.1
            if (is.null(adright.1)) print("Alternative Acceptor...Type I...none...")
            else print("Alternative Acceptor...Type I...FOUND...")
      }
     
     if ("adleft.2" %in% eventType || "adright.2" %in% eventType)
      {
            adleftRaw.2 <- asEvent(cis.read, left.jun, right.jun,
                                              sampleName, sampleID, minReadCounts,
                                              eventType="adleft.2", novel, readSum$totalRead, 
                                              select.GRange
                                              )
            adrightRaw.2 <- asEvent(cis.read, left.jun, right.jun,
                sampleName, sampleID, 
                minReadCounts,
                eventType="adright.2", novel, readSum$totalRead, 
                select.GRange
                )
 
            list.2 <- swapStrand(adleftRaw.2, adrightRaw.2)
            adleft.2 <- list.2$adleft
            adright.2 <- list.2$adright 

            adleft.2 <- annotateNovel(adleft.2, "adleft.2", novel, exonStart, exonEnd, splicingLink)  
            adright.2 <- annotateNovel(adright.2, "adright.2", novel, exonStart, exonEnd, splicingLink)  
               
            result.list[["adleft.2"]] <- adleft.2
            if (is.null(adleft.2)) print("Alternative Donor...Type II...none...")
            else print("Alternative Donor...Type II...FOUND...")

            result.list[["adright.2"]] <- adright.2
            if (is.null(adright.2)) print("Alternative Acceptor...Type II...none...")
            else print("Alternative Acceptor...Type II...FOUND...")
      }
  
      if ("adboth.1" %in% eventType)
      {                
             adboth.1 <- asEvent(cis.read, left.jun, right.jun,
                  sampleName, sampleID, 
                  minReadCounts,
                  eventType="adboth.1", novel, readSum$totalRead, 
                  select.GRange
                  )
            adboth.1 <- annotateNovel(adboth.1, "adboth.1", novel, exonStart, exonEnd, splicingLink)  

            result.list[["adboth.1"]] <- adboth.1
           if (is.null(adboth.1)) print("Alternative Both Sites...Type I...none...")
           else print("Alternative Both Sites...Type I...FOUND...")
      }
      if ("adboth.2" %in% eventType)
      {
             adboth.2 <- asEvent(cis.read, left.jun, right.jun,
                  sampleName, sampleID, 
                  minReadCounts,
                  eventType="adboth.2", novel, readSum$totalRead, 
                  select.GRange
                  )
            adboth.2 <- annotateNovel(adboth.2, "adboth.2", novel, exonStart, exonEnd, splicingLink)  
                              
            result.list[["adboth.2"]] <- adboth.2
           if (is.null(adboth.2)) print("Alternative Both Sites...Type II...none...")
           else print("Alternative Both Sites...Type II...FOUND...")
      }
             
  return(list(type=result.list, summary=readSum))

}



