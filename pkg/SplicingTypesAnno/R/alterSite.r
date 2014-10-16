##################################
##  SplicingTypesAnno
##  - alterSite.r: main algorithms for alternative site,
##              including adleft, adright, adboth
##################################

## convert unknown gene range to a known gene range
## - othersInfo: return also start&end info
convertKnownName <- function(unknown.GRange, known.GRange, otherInfo=FALSE, etype, gene.GRange=NULL)
{

   match.m <- as.matrix(findOverlaps(unknown.GRange, known.GRange)) 
   nondup.ind <- which(!duplicated(match.m[,1]))
   #if (length(nondup.ind) != length(known.GRange)) print("Some unknown are NOT annotated!")
 
   nondup.m <- match.m[nondup.ind,,drop=FALSE]
   known.anno <- elementMetadata(known.GRange)[, "geneName"]
   if(etype=="ri") unknown.anno <- rep("riInsideReducedExon", length(unknown.GRange))
   else if (etype=="es") unknown.anno <- rep("esInsideReducedIntron", length(unknown.GRange))
   else unknown.anno <- rep("UNKNOWN", length(unknown.GRange))
   unknown.anno[nondup.m[,1]] <- known.anno[nondup.m[,2]]
   
   elementMetadata(unknown.GRange)[, "geneName"] <- unknown.anno

   ## annotate geneName if not overlap with known
   if (!is.null(gene.GRange))
   {
      nonnovel.ind <- unique(match.m[,1]) 
      
      if (length(nonnovel.ind)>0)
      {
          novel.GRange <- unknown.GRange[-nonnovel.ind]
          match.g <- as.matrix(findOverlaps(novel.GRange, gene.GRange)) 
          g.ind <- which(!duplicated(match.g[,1]))
 
          g.m <- match.g[g.ind,,drop=FALSE]
          novel.anno <- elementMetadata(novel.GRange)[, "geneName"]
          gene.anno <- elementMetadata(gene.GRange)[, "geneName"]
          novel.anno[g.m[,1]] <- paste(novel.anno[g.m[,1]], gene.anno[g.m[,2]], sep=";")
          
          unknown.anno.1 <- elementMetadata(unknown.GRange)[, "geneName"]
          unknown.anno.1[-nonnovel.ind] <- novel.anno
          elementMetadata(unknown.GRange)[, "geneName"] <- unknown.anno.1
                   
      }
      
   }
   
   if (!otherInfo) return(unknown.GRange)
   else
   {
      known.start <- start(known.GRange)
      unknown.start <- rep(0, length(unknown.GRange))

      unknown.start[nondup.m[,1]] <- known.start[nondup.m[,2]]
      
      known.end <- end(known.GRange)
      unknown.end <- rep(0, length(unknown.GRange))

      unknown.end[nondup.m[,1]] <- known.end[nondup.m[,2]]
   }
   return(list(readGRange=unknown.GRange, 
               other=list(start=unknown.start, end=unknown.end
                          )))
}

## Main function: break junction reads into left and right reads separately
breakJunRead <- function(myread)
{         
    ## simplify the question, only gap 1 allowed!
    myread <- myread[ngap(myread)==1]

    mycigar <- cigar(myread)
    mycigar.tmp <- paste("top=", mycigar, sep="")   
    mycigar <- paste(mycigar.tmp , "=end", sep="")
   
    result <- unlist(strsplit(mycigar, "N"))

    endstr <- result[grep("=end", result)]
    startstr <- result[grep("top=", result)]

    if (length(startstr)==0 || length(endstr)==0) stop("something wrong with junction cigar!")
    
    start.top <- unlist(strsplit(startstr, "M"))
    start.tmp <- start.top[grep("top=", start.top)]
    
    start.done <- gsub("top=", "", start.tmp)
    start.ind <- grep("[a-zA-Z]", start.done) 
    
    if (length(start.ind)!=0) start.done <- as.numeric(start.done[-start.ind])
    else start.done <- as.numeric(start.done)
    
    end.done <- gsub("M=end", "", endstr)
    end.ind <- grep("[a-zA-Z]", end.done)

    na.ind <- c(start.ind, end.ind)
    
    if (length(na.ind)>=1) 
    {
        start.done <- start.done[-na.ind]
        end.done <- end.done[-na.ind]
    }

    if (length(na.ind)!=0) final.read <- myread[-na.ind]
    else final.read <- myread 

      left.jun <- GRanges(seqnames=seqnames(final.read),
                            strand=strand(final.read),
                            ranges=IRanges(start=start(final.read),
                              end=start(final.read) + as.numeric(start.done) - 1)
                              )

     right.jun <- GRanges(seqnames=seqnames(final.read),
                            strand=strand(final.read),
                            ranges=IRanges(start=end(final.read) - as.numeric(end.done) + 1,
                              end=end(final.read) )
                              )
     middle.jun <- GRanges(seqnames=seqnames(final.read),
                            strand=strand(final.read),
                            ranges=IRanges(start=start(final.read) + as.numeric(start.done),
                              end=end(final.read) - as.numeric(end.done) )
                              )
     
     
     if(length(left.jun) != length(final.read)) print("breaking left read FAILs!")
     if(length(right.jun) != length(final.read)) print("breaking right read FAILs!")
                                
     return(list(left=left.jun, right=right.jun, middle=middle.jun, read=final.read))
}


## return all edge-duplicated reads with >= cutoff
copyDfCheck <- function(read, startOrEnd, cutoff)
{
      if (startOrEnd=="start") left.end <- start(read)      
      else left.end <- end(read)
      
      bound.str <- paste(chr=as.vector(seqnames(read)), pos=left.end, strand=as.vector(strand(read)), sep="_")

      left.dup.s <- (table(bound.str) >= cutoff )
      left.dup.true <- left.dup.s[left.dup.s]
      left.end.select <- unlist(attr(left.dup.true, "dimnames"))

      final.ind <- which(bound.str %in% left.end.select)

      return(final.ind)
}

## return all edge-duplicated reads with >= 2
copyStrCheck <- function(right.jun)
{
      right.str <- paste(chr=as.vector(seqnames(right.jun)), pos=start(right.jun), strand=as.vector(strand(right.jun)), sep="_")

      right.ind <- which(duplicated(right.str))
      return(right.ind)
}

## return unique index without duplicate for both sides
copyTwoBothInd <- function(left.jun.2, right.jun.2)
{
      copy.str <- paste(chr=as.vector(seqnames(left.jun.2)),
                            pos=end(left.jun.2),
                            strand=as.vector(strand(left.jun.2)),
                            pos=start(right.jun.2),
                            sep="_")
       copy.ind <- which(!duplicated(copy.str))
       return(copy.ind)
}

## adleft, adright
## Main goal: return unique read with
#  - a. left matched
#  - b. right matched
#  - c. unique no match for two sides: useful for both alternative
adSingleSide <- function(cis.read, left.jun, right.jun, minReadCounts, adType="rightAlternative")
{
  ## duplicate read index
  if (adType=="rightAlternative")
  {
      left.ind <- copyDfCheck(left.jun, "end", minReadCounts * 2)
      right.ind <- copyStrCheck(right.jun)
  }
  else
  {
      right.ind <- copyDfCheck(right.jun, "start", minReadCounts * 2)
      left.ind <- copyStrCheck(left.jun)
  }

  done.ind <- intersect(left.ind, right.ind)

  ## reads - pick unique one combined left and right side
  ## Main role: reduce the no of counts! too many for memory
  ds.read <- cis.read[done.ind]
  ds.tmp.2 <- breakJunRead(ds.read)
  left.jun.2 <- ds.tmp.2$left
  right.jun.2 <- ds.tmp.2$right
  ds.read <- ds.tmp.2$read

  final.ind <-   copyTwoBothInd(left.jun.2, right.jun.2)
  final.read <- ds.read[final.ind]
  left.final <- left.jun.2[final.ind]
  right.final <- right.jun.2[final.ind]

   ## pick one
   ## Note: if one side matches, automatically the other side will NOT match
   ds.left.c <- countOverlaps(left.final, left.final, type="end")
   ds.right.c <- countOverlaps(right.final, right.final, type="start")

  if (adType=="rightAlternative")
  {
      candidate.ind <- which( ds.left.c > 1)
  }
  else
  {
      candidate.ind <- which( ds.right.c > 1)
  }

   finalleft.junread <- left.final[candidate.ind]
   finalright.junread <- right.final[candidate.ind]
   read <- final.read[candidate.ind]
  
   ## limit range to mutilple reads across one/few introns. But each side overlap same locus
   ## NOTE: the locus = exon + following intron

   return(list(left=finalleft.junread, right=finalright.junread, read=read))
}

## adboth
## Main function: return unique read
adbothSide <- function(cis.read, left.jun, right.jun, minReadCounts)
{
      
  ## duplicate index
  left.ind <- copyStrCheck(left.jun)
  right.ind <- copyStrCheck(right.jun)

  done.ind <- intersect(left.ind, right.ind)

  ## reads - pick unique one combined left and right side
  ## Main role: reduce the no of counts! too many for memory
  ds.read <- cis.read[done.ind]
  ds.tmp.2 <- breakJunRead(ds.read)
  left.jun.2 <- ds.tmp.2$left
  right.jun.2 <- ds.tmp.2$right
  ds.read <- ds.tmp.2$read

  ## unique index - only leave one for further, save memory
  final.ind <-   copyTwoBothInd(left.jun.2, right.jun.2)
  final.read <- ds.read[final.ind]
  left.final <- left.jun.2[final.ind]
  right.final <- right.jun.2[final.ind]

   ## differ from single side

   # take off single side
   ds.left.c <- countOverlaps(left.final, left.final, type="end")
   ds.right.c <- countOverlaps(right.final, right.final, type="start")

   candidate.ind <- intersect(which(ds.left.c < 2), which(ds.right.c < 2))
   ## also the reads must overlap at leas in some region. certainly NOT either side.
   ds.overlap <- countOverlaps(final.read, final.read)
   overlap.ind <- which(ds.overlap > 1)
   
   candidate.ind <- intersect(candidate.ind, overlap.ind)
   if (length(candidate.ind)==0) return(invisible(NULL))
   left.done <- left.final[candidate.ind]
   right.done <- right.final[candidate.ind]
   read <- final.read[candidate.ind]

   return(list(left=left.done, right=right.done, read=read))
}

# select type: adleft.1, adleft.2, adright.1, adright.2
## input: unique read representative. it has multiple copies 
## output: unique jun read
typeSelection <- function(left.donor, eventType, intronGRange)
{
    if (is.null(left.donor)) return(invisible(NULL))
    left.jun.type <- left.donor$left
    right.jun.type <- left.donor$right
    read <- left.donor$read

    tmp.c <- countOverlaps(granges(read), intronGRange)
    
    if (eventType %in% c("adleft.1", "adright.1", "adboth.1"))
    {
        # keep those case: one in type I, the other in type II
        if (eventType =="adleft.1")
        {
            tmp.right.jun <- right.jun.type[tmp.c < 2]
            tmax <- as.matrix(findOverlaps(tmp.right.jun, right.jun.type, type="start"))
            tind <- unique(tmax[,2])
            
            final.read <- read[tind] # hard code for the type I
            final.left.jun <- left.jun.type[tind]
            final.right.jun <- right.jun.type[tind]            
            
        }
        else if (eventType =="adright.1")
        {              
            tmp.left.jun <- left.jun.type[tmp.c < 2]
            tmax <- as.matrix(findOverlaps(tmp.left.jun, left.jun.type, type="end"))
            tind <- unique(tmax[,2])
            
            final.read <- read[tind] # hard code for the type I
            final.left.jun <- left.jun.type[tind]
            final.right.jun <- right.jun.type[tind]             
        }
        else
        {   ## adbothSide function has deleted all start/end matching cases                     
            final.read <- read[tmp.c < 2] # hard code for the type I
            final.left.jun <- left.jun.type[tmp.c < 2]
            final.right.jun <- right.jun.type[tmp.c < 2]    
        }

                
    }
    else if (eventType %in% c("adleft.2", "adright.2", "adboth.2"))
    {               
        # NOT keep the special case: one in type I, the other in type II
        # select seeds
        seed.read <- read[tmp.c >= 2]
        seed.left.jun <- left.jun.type[tmp.c >= 2]
        seed.right.jun <- right.jun.type[tmp.c >= 2]

        # select all jun matching seeds
        if (eventType =="adleft.2")
        {
            seed.c <- countOverlaps(seed.right.jun, right.jun.type, type="start")
            
            ## > 1 take off those cases: type II match type I in start/end
            final.left.jun <- seed.left.jun[seed.c > 1]
            final.right.jun <- seed.right.jun[seed.c > 1]
            final.read <- seed.read[seed.c > 1]
        }
        else if (eventType == "adright.2")
        {                
            seed.c <- countOverlaps(seed.left.jun, left.jun.type, type="end")
            final.left.jun <- seed.left.jun[seed.c > 1]
            final.right.jun <- seed.right.jun[seed.c > 1]
            final.read <- seed.read[seed.c > 1]
        }
        else 
        {
            final.left.jun <- seed.left.jun
            final.right.jun <- seed.right.jun
            final.read <- seed.read
        }
       
    }
    else stop("Alternative selection does NOT have this type!")
    
    return(list(left=final.left.jun, right=final.right.jun, read=final.read))
} 

## adleft, adright
## Main function: for single side, calcRatio
getStartEnd <- function(pos)
{
              tmp <- unique(pos)
              if (length(tmp)==1) return(tmp)
              else return(paste(tmp, collapse="_"))
}

## adboth
## for both side, calcBothRatio
getString <- function(string)
{
    tmp <- paste(string, collapse="_")
    return(tmp)
}

## adleft, adright
## Main function: quantify adleft, adright
calcRatio <- function(left.jun, right.jun, right.donor, adType="rightAlternative", 
                      minReadCounts, gRange, novel)
{
  if (missing(gRange)) 
  {
      gene.GRange <- splicingENV$gene.GRange
      exonGRange <- splicingENV$exonGRange
      intronGRange <- splicingENV$intronGRange
      splicingLink <- splicingENV$splicingLink
  }
  else
  {
      gene.GRange <- gRange$gene.GRange
      exonGRange <- gRange$exonGRange
      intronGRange <- gRange$intronGRange
      splicingLink <- gRange$splicingLink
  }
               
   right.jun.cal <- right.donor$right
   left.jun.cal <- right.donor$left
   
   if (length(right.jun.cal)==0 || length(left.jun.cal)==0)
   {
      #print(paste("no value in alternative donor in ", adType, sep=""))
      return(NULL)
   }
      ## read too many, so use the selected THE unique read to count all reads matching it EXACTLY
      ## purpose: 1. counts reads for each unique case
      ##          2. anno unique case with known anno        
      mergeJun <- function(left.jun.cal, right.jun.cal, left.jun, right.jun, mtype)
      {            
          # count overlap reads in "right"
          all.str <- paste(chr=as.vector(seqnames(left.jun)),
                            pos=end(left.jun),
                            strand=as.vector(strand(left.jun)),
                        #chrR=as.vector(seqnames(right.jun.2)),
                            pos=start(right.jun),
                            sep="_")
          jun.str <- paste(chr=as.vector(seqnames(left.jun.cal)),
                            pos=end(left.jun.cal),
                            strand=as.vector(strand(left.jun.cal)),
                        #chrR=as.vector(seqnames(right.jun.2)),
                            pos=start(right.jun.cal),
                            sep="_")

          uniqueCount <- table(all.str)
          junCount <- uniqueCount[names(uniqueCount) %in% jun.str]
          
          ## df.1: string from all reads
          df.1 <- data.frame(str=names(junCount), count=junCount)
          ## df.0, string from target unique read
          df.0 <- data.frame(ID=c(1:length(left.jun.cal)), str=jun.str)

          df.10 <- merge(df.1, df.0, by="str")

          # overlap between read versus reduce(exon)
          all.reduce <- sort(c(reduce(exonGRange), reduce(intronGRange)))

          if (mtype=="end")
          {
              matrix.1 <- as.matrix(findOverlaps(left.jun.cal, all.reduce))
          }
          else
          {
              matrix.1 <- as.matrix(findOverlaps(right.jun.cal, all.reduce)) 
          }
          df.2 <- data.frame(matrix.1)
          colnames(df.2) <- c("ID", "reduceID")

          ## right counts
          df.done <- merge(df.10, df.2, by="ID")
          df.done <- df.done[, c("ID", "count", "reduceID")]
          return(df.done)
      }

       if (adType=="rightAlternative")
       {    
          
          # left sum
          junleft.df <- mergeJun(left.jun.cal, right.jun.cal, left.jun, right.jun, "end")
          colnames(junleft.df) <- c("ID", "count", "reduceIDLeft")
          if (nrow(junleft.df)==0) return(NULL)

          # find overlap in unique reads
          umatrix <- as.matrix(findOverlaps(left.jun.cal, left.jun.cal, type="end") )
          neg.ind <- which(umatrix[,2]-umatrix[,1]<0)
          if (length(neg.ind)>0)
          {    
            not.include <- umatrix[neg.ind, , drop=FALSE][,1]
            not.ind <- which(umatrix[,1] %in% not.include)
            umatrix <- umatrix[-not.ind,]
          }
          
          j.ind <- match(umatrix[,2], junleft.df$ID)
          #tmp.table <- data.frame(subject=umatrix[,1], junleft.df[j.ind,])
          start.id <- umatrix[,1]
          end.id <- umatrix[,2]
          tmp.table <- data.frame(subject=umatrix[,1], junleft.df[j.ind,], 
                                  junLeftEnd=end(left.jun.cal)[start.id],
                                  junRightStart=start(right.jun.cal)[end.id]
                                  )
          ## filter cufoff reads - step 1
          filter.ind <- which( tmp.table$count >= minReadCounts)
          if (length(filter.ind)>0) tmp.table <- tmp.table[filter.ind,]
          else return(NULL)    

          count.sum <- tapply(tmp.table$count, factor(tmp.table$subject), sum)
          count.max <- tapply(tmp.table$count, factor(tmp.table$subject), max)
          
          #########################################
          ## TODO: start/end details
          ## Details later for IGV visualization
          count.start <- tapply(tmp.table$junLeftEnd, factor(tmp.table$subject), getStartEnd)
          count.scount <- tapply(tmp.table$count, factor(tmp.table$subject), getStartEnd)
          count.end <- tapply(tmp.table$junRightStart, factor(tmp.table$subject), getStartEnd)

          ## END #######################################
                    
          df.sum <- data.frame(subject=names(count.sum), countSum=count.sum)
          df.max <- data.frame(subject=names(count.max), countMax=count.max)
          
          #df.table <- merge(tmp.table, df.sum, by="subject")
          #df.table <- merge(df.table, df.max, by="subject")

          df.tmp <- data.frame(subject=names(count.sum),
                                 countSum=count.sum,
                                 countMax=count.max,
                                 junLeftEndCollection=as.character(count.start),
                                 junRightStartCollection=paste(count.end, count.scount, sep=",")
                                 , stringsAsFactors = FALSE
                                 )
          df.table <- merge(tmp.table, df.tmp, by="subject")
          u.ind <- unique(tapply(df.table$ID, df.table$subject, min))

          ## NOTE: df.table has all the read mapping information in that location
          #tmp2.table <- df.table[, c("subject", "reduceIDRight", "countSum", "countMax")]
          final.table.full <- df.table
          final.table <- unique(final.table.full[final.table.full$ID %in% u.ind,])

          ## final.table.full: include all cases.
          ## final.table: collapse those same events to 1 event number
          final.table.full$ratio <- round(final.table.full$count/(final.table.full$countSum),4)
          final.table$ratio <- round(final.table$count/final.table$countSum ,4)
      }
       else
       {  ## 5' alter
                  
          # right sum
          junright.df <- mergeJun(left.jun.cal, right.jun.cal, left.jun, right.jun, "start")
          colnames(junright.df) <- c("ID", "count", "reduceIDRight")
          if (nrow(junright.df)==0) return(NULL)
          
          # find overlap in unique reads
          umatrix <- as.matrix(findOverlaps(right.jun.cal, right.jun.cal, type="start") )
          neg.ind <- which(umatrix[,2]-umatrix[,1]<0)
          if (length(neg.ind)>0)
          {
            not.include <- umatrix[neg.ind,  , drop=FALSE][,1]
            not.ind <- which(umatrix[,1] %in% not.include)
            umatrix <- umatrix[-not.ind,]
          }

          j.ind <- match(umatrix[,2], junright.df$ID)
          #tmp.table <- data.frame(subject=umatrix[,1], junleft.df[j.ind,])
          end.id <- umatrix[,1]
          start.id <- umatrix[,2]
          tmp.table <- data.frame(subject=umatrix[,1], junright.df[j.ind,], 
                                  junLeftEnd=end(left.jun.cal)[start.id],
                                  junRightStart=start(right.jun.cal)[end.id]
                                  )
          ## filter cufoff reads - step 1
          filter.ind <- which( tmp.table$count >= minReadCounts)
          if (length(filter.ind)>0) tmp.table <- tmp.table[filter.ind,]
          else return(NULL)          

          ##
          count.sum <- tapply(tmp.table$count, factor(tmp.table$subject), sum)
          count.max <- tapply(tmp.table$count, factor(tmp.table$subject), max)
          
          #########################################
          ## TODO: start/end details
          ## Details later for IGV visualization
          count.start <- tapply(tmp.table$junLeftEnd, factor(tmp.table$subject), getStartEnd)
          count.scount <- tapply(tmp.table$count, factor(tmp.table$subject), getStartEnd)
          count.end <- tapply(tmp.table$junRightStart, factor(tmp.table$subject), getStartEnd)

          ## END #######################################
          
          df.sum <- data.frame(subject=names(count.sum), countSum=count.sum)
          df.max <- data.frame(subject=names(count.max), countMax=count.max)
         
          df.tmp <- data.frame(subject=names(count.sum),
                                 countSum=count.sum,
                                 countMax=count.max,
                                 junLeftEndCollection=paste(count.start, count.scount, sep=","),
                                 junRightStartCollection=as.character(count.end)
                                 , stringsAsFactors = FALSE
                                 )
          df.table <- merge(tmp.table, df.tmp, by="subject")
          u.ind <- unique(tapply(df.table$ID, df.table$subject, min))
          
          ## NOTE: df.table has all the read mapping information in that location
          #tmp2.table <- df.table[, c("subject", "reduceIDRight", "countSum", "countMax")]          
          final.table.full <- df.table
          final.table <- unique(final.table.full[final.table.full$ID %in% u.ind,])

          ## final.table.full: include all cases.
          ## final.table: collapse those same events to 1 event number
          final.table.full$ratio <- round(final.table.full$count/(final.table.full$countSum),4)
          final.table$ratio <- round(final.table$count/final.table$countSum,4)
       }
                                      

       ## filter cufoff reads - step 2
       filter.ind <- which( (final.table.full$countSum - final.table.full$countMax) != 0 )
       if (length(filter.ind)>0) final.table.full <- final.table.full[filter.ind,]
       else return(NULL)
            
       return(unique(final.table.full))
       
}

## adboth
## Main function: quantify adboth
calcRatioBoth <- function(left.jun, right.jun, both.donor, minReadCounts, gRange, novel)
{
  if (missing(gRange)) 
  {
      gene.GRange <- splicingENV$gene.GRange
      exonGRange <- splicingENV$exonGRange
      intronGRange <- splicingENV$intronGRange
      splicingLink <- splicingENV$splicingLink
  }
  else
  {
      gene.GRange <- gRange$gene.GRange
      exonGRange <- gRange$exonGRange
      intronGRange <- gRange$intronGRange
      splicingLink <- gRange$splicingLink
  }
           
   right.jun.cal <- both.donor$right
   left.jun.cal <- both.donor$left

   if (length( right.jun.cal)!= 0 && length(left.jun.cal)!=0 )
   {

      ## calculate counts for each unique jun read
      mergeJun <- function(left.jun.cal, right.jun.cal, left.jun, right.jun)
      {                 
          # count overlap reads in "right"
          all.str <- paste(chr=as.vector(seqnames(left.jun)),
                            pos=end(left.jun),
                            strand=as.vector(strand(left.jun)),
                            pos=start(right.jun),
                            sep="_")
          jun.str <- paste(chr=as.vector(seqnames(left.jun.cal)),
                            pos=end(left.jun.cal),
                            strand=as.vector(strand(left.jun.cal)),
                            pos=start(right.jun.cal),
                            sep="_")                           

          uniqueCount <- table(all.str)
          junCount <- uniqueCount[names(uniqueCount) %in% jun.str]
          df.1 <- data.frame(str=names(junCount), count=junCount)
          df.0 <- data.frame(ID=c(1:length(left.jun.cal)), 
                        str=jun.str, 
                        junLeftEnd= end(left.jun.cal),
                        junRightStart=start(right.jun.cal)
                        )

          df.10 <- merge(df.1, df.0, by="str")

          # overlap between read versus reduce(exon)
          all.reduce <- sort(c(reduce(exonGRange), reduce(intronGRange)))

          matrix.left <- data.frame(as.matrix(findOverlaps(left.jun.cal, all.reduce)))
          matrix.left <- matrix.left[!duplicated(matrix.left[,1]),]
          colnames(matrix.left) <- c("ID", "reduceIDLeft")
          matrix.right <- data.frame(as.matrix(findOverlaps(right.jun.cal, all.reduce)) )
          colnames(matrix.right) <- c("ID", "reduceIDRight")
          
          df.2 <- merge(matrix.left, matrix.right, by="ID")
          
          df.done <- merge(df.10, df.2, by="ID")

          return(df.done)
      }

          # left sum
      jun.df <- mergeJun(left.jun.cal, right.jun.cal, left.jun, right.jun)
      jun.df$bindstr <- paste(jun.df$reduceIDLeft, jun.df$reduceIDRight, sep="_")
      
      ## check both donor side exist
      dup.str <- jun.df$bindstr[duplicated(jun.df$bindstr)]
      jun.final <- jun.df[jun.df$bindstr %in% dup.str,]
      if (nrow(jun.final)==0) return(NULL)
      
      jun.final$finalstr <- paste(jun.final$str, jun.final$count, sep=",")
      
      ## filter cufoff reads - step 1
      filter.ind <- which( jun.final$count >= minReadCounts)
      if (length(filter.ind)>0) jun.final <- jun.final[filter.ind,]
      else return(NULL)          
                

      count.max <- tapply(jun.final$count, factor(jun.final$bindstr), max)
      count.sum <- tapply(jun.final$count, factor(jun.final$bindstr), sum)
      count.string <- tapply(jun.final$finalstr, factor(jun.final$bindstr), getString)

      count.df <- data.frame(bindstr=names(count.max),
                             countMax=count.max,
                             countSum=count.sum,
                             junLeftRightCollection=count.string)     
      final.df <- merge(count.df, jun.final, by="bindstr")
    
      final.df$ratio <- round(final.df$count/final.df$countSum,4)
      #final.table <- final.df[, c("reduceIDLeft", "reduceIDRight", "countMax", "countSum", "ratio")]
       final.table <- final.df
       filter.ind <- which( (final.table$countSum - final.table$countMax) != 0)
       if (length(filter.ind)>0) final.table <- final.table[filter.ind,]
       else return(NULL)

       return(unique(final.table))
  }
  else return(NULL)
}

## final output format of adleft, adright, adboth
writeADEvent <- function(ratioTable, eventType, all, all.reduce, sampleID=sampleID, sampleName=sampleName, total.read=total.read)
{
      if (is.null(ratioTable)) 
      {
          #print(paste(eventType, " is NOT exist!", sep=""))
          return(NULL)
      }

      ## useful for ad? rpkm normalized counts
      #if (total.read!=0) 
      #{
          #scaleFac <-(10^6)/total.read
          #ratioTable$countSum.norm <- ratioTable$countSum * scaleFac
          #ratioTable$countMax.norm <- ratioTable$countMax * scaleFac
      #}
                  
      if (eventType %in% c("adboth.1", "adboth.2"))
      {
                 
        left.ind <- ratioTable$reduceIDLeft
        right.ind <- ratioTable$reduceIDRight

        left.g <- all.reduce[left.ind]
        right.g <- all.reduce[right.ind]
        
        left.g <- convertKnownName(left.g, all, etype="all")
        right.g <- convertKnownName(right.g, all, etype="all")

        ## add "realJunLocus" to be consistent for ri, es, and ad
        ratioTable$realJunLocus <- paste(paste(as.vector(seqnames(left.g)),
                                             paste(ratioTable$junLeftEnd,ratioTable$junRightStart,sep="-"), sep=":"),
                                             sep="_")
                                             
        ratio.ind <- match(c("count","junLeftEnd","junRightStart",
                                       "countSum","countMax",
                                       "junLeftRightCollection",
                                       "realJunLocus",
                                       "ratio",
                                       "novelSplicingLink",
                                       "novelLeft", "novelRight"
                                       ),
                            colnames(ratioTable))
        ratio.ind <- ratio.ind[!is.na(ratio.ind)]
        targetTable <- ratioTable[,ratio.ind]
                                    
        ad.final <- data.frame(mergeID=paste(paste(as.vector(seqnames(left.g)),
                                             paste(ratioTable$junLeftEnd,ratioTable$junRightStart,sep="-"), sep=":"),
                                             as.vector(strand(left.g)),
                                             sep="_"),
                               geneName=as.character(unlist(values(left.g)[, "geneName"])),
                               chr=as.vector(seqnames(left.g)),
                               locus=paste(as.vector(seqnames(left.g)),
                                       paste(ratioTable$junLeftEnd,ratioTable$junRightStart,sep="-"), sep=":"),
                               #strand= as.vector(strand(left.g)),
                               #leftExon=unlist(values(left.g)),
                               #rightExon=unlist(values(right.g)),
                                #max=ratioTable$countMax,
                                #others=ratioTable$countSum - ratioTable$countMax,
                                #ratio=ratioTable$ratio
                                targetTable,
                                stringsAsFactors = FALSE 
                                )
                           
        if (nrow(ad.final)==0) return(NULL)
        ad.final$ratio <- format(ad.final$ratio, scientific = FALSE)
        
        ## for bed file - start
          adTable <- ad.final[, c("chr", "junLeftEnd", "junRightStart")]
          adTable$ratio <- paste(ad.final$ratio, "=", 
                                 ad.final$junLeftRightCollection, 
                                         "_", sampleName, "_", eventType, sep="")
   
          ad.final$note2 <- adTable$ratio
          eventType <- gsub("\\.", "_", eventType)
          writeBed(adTable, "alternative",
                         paste(sampleName, "_", eventType, ".bed", sep="")) 
       ## for bed file - end
                         
       colnames(ad.final) <-  c("mergeID", "geneName", "chr", "locus",
                               paste(sampleName, sampleID, colnames(ad.final)[-c(1:4)], sep=".")
                               )
                                                                                                                 

      }
      else if  (eventType %in% c("adleft.1", "adleft.2") )
      {            
        right.ind <- ratioTable$reduceIDRight
        right.g <- all.reduce[right.ind]
        right.g <- convertKnownName(right.g, all, etype="all")

        ## add "realJunLocus" to be consistent for ri, es, and ad
        ratioTable$realJunLocus <- paste(paste(as.vector(seqnames(right.g)),
                                             paste(ratioTable$junLeftEnd,ratioTable$junRightStart,sep="-"), sep=":"),
                                             sep="_")

        ## TOCheck: merge by right exon name! instead of gene Name
        ad.final <- data.frame(#geneName=as.character(unlist(values(right.g))),
                               #rightExon=unlist(values(right.g)),
                               mergeID=paste(paste(as.vector(seqnames(right.g)),
                                             paste(ratioTable$junLeftEnd,ratioTable$junRightStart,sep="-"), sep=":"),
                                             as.vector(strand(right.g)),
                                             sep="_"),
                               geneName= unlist(values(right.g)[, "geneName"]),
                               chr=as.vector(seqnames(right.g)),
                               locus=paste(as.vector(seqnames(right.g)),
                                       paste(ratioTable$junLeftEnd,ratioTable$junRightStart,sep="-"), sep=":"),
                               #strand= as.vector(strand(right.g)),
                                #max=ratioTable$countMax,
                                #others=ratioTable$countSum - ratioTable$countMax,
                                #ratio=ratioTable$ratio
                                ratioTable[, -c(1:2,4)],
                                stringsAsFactors = FALSE 
                                )
  
        if (nrow(ad.final)==0) return(NULL)
        ad.final$ratio <- format(ad.final$ratio, scientific = FALSE)
        
        ## for bed file-start
          adTable <- ad.final[, c("chr", "junLeftEnd", "junRightStart")]
          adTable$ratio <- paste(ad.final$ratio, "=",
                                         ad.final$junLeftEndCollection,
                                         "_", sampleName, "_", eventType, sep="")
          
          ad.final$note2 <- adTable$ratio
          
          eventType <- gsub("\\.", "_", eventType) 

          ad.strand <- as.vector(strand(right.g))          
          adTable.p <- adTable[ad.strand=="+",]
          adTable.n <- adTable[ad.strand=="-",]
          if (eventType=="adleft_1")
          {              
              writeBed(adTable.p, "alternative",
                         paste(sampleName, "_", "donor_1", ".bed", sep=""))
              writeBed(adTable.n, "alternative",
                         paste(sampleName, "_", "acceptor_1", ".bed", sep=""))                           
          }
          else
          {
              writeBed(adTable.p, "alternative",
                         paste(sampleName, "_", "donor_2", ".bed", sep=""))
              writeBed(adTable.n, "alternative",
                         paste(sampleName, "_", "acceptor_2", ".bed", sep=""))          
          }
         ## for bed file-end
                         
          colnames(ad.final) <-  c("mergeID", "geneName", "chr", "locus",
                                paste(sampleName, sampleID, colnames(ad.final)[-c(1:4)], sep=".")                                                                   
                                )
      }
      else
      {
               
        left.ind <- ratioTable$reduceIDLeft
        left.g <- all.reduce[left.ind]
        left.g <- convertKnownName(left.g, all, etype="all")

        ## add "realJunLocus" to be consistent for ri, es, and ad
        ratioTable$realJunLocus <- paste(paste(as.vector(seqnames(left.g)),
                                             paste(ratioTable$junLeftEnd,ratioTable$junRightStart,sep="-"), sep=":"),
                                             sep="_")
        ad.final <- data.frame(#geneName=as.character(unlist(values(left.g))),
                               #leftExon=unlist(values(left.g)),
                               mergeID= paste(paste(as.vector(seqnames(left.g)),
                                              paste(ratioTable$junLeftEnd,ratioTable$junRightStart,sep="-"), sep=":"),
                                              as.vector(strand(left.g)),
                                              sep="_"),
                               geneName= unlist(values(left.g)[, "geneName"]),
                               chr = as.vector(seqnames(left.g)), 
                               locus=paste(as.vector(seqnames(left.g)),
                                       paste(ratioTable$junLeftEnd,ratioTable$junRightStart,sep="-"), sep=":"),
                               #strand= as.vector(strand(left.g)),
                                #max=ratioTable$countMax,
                                #others=ratioTable$countSum - ratioTable$countMax,
                                #ratio=ratioTable$ratio
                                ratioTable[, -c(1:2,4)],
                                stringsAsFactors = FALSE 
                                )  
        if (nrow(ad.final)==0) return(NULL)
        ad.final$ratio <- format(ad.final$ratio, scientific = FALSE)
        
        ## for bed file - start
          adTable <- ad.final[, c("chr", "junLeftEnd", "junRightStart")]
          adTable$ratio <- paste(ad.final$ratio, "=",
                                         ad.final$junRightStartCollection,
                                         "_", sampleName, "_", eventType, sep="")
          ad.final$note2 <- adTable$ratio
          eventType <- gsub("\\.", "_", eventType) # change strand .
          
          ad.strand <- as.vector(strand(left.g))
          adTable.p <- adTable[ad.strand=="+",]
          adTable.n <- adTable[ad.strand=="-",]
          if (eventType=="adright_1")
          {              
              writeBed(adTable.p, "alternative",
                         paste(sampleName, "_", "acceptor_1", ".bed", sep=""))
              writeBed(adTable.n, "alternative",
                         paste(sampleName, "_", "donor_1", ".bed", sep=""))                           
          }
          else
          {
              writeBed(adTable.p, "alternative",
                         paste(sampleName, "_", "acceptor_2", ".bed", sep=""))
              writeBed(adTable.n, "alternative",
                         paste(sampleName, "_", "donor_2", ".bed", sep=""))          
          }
                 
        ## for bed file - end
        
       colnames(ad.final) <-  c("mergeID", "geneName", "chr", "locus",
                                paste(sampleName, sampleID, colnames(ad.final)[-c(1:4)], sep=".")                                                                   
                                )
                                   
      }
      #print("alternative donor...processed...")
             
      return(ad.final)

}
