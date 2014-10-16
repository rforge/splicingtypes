##################################
##  SplicingTypesAnno
##  - algorithmFunction.r: heuristic details
##    a. details for RI; b. details for ES; c. simple link to alternative site
##################################

## select genes for analysis
convertParamStr <- function(selectGRange)
{
      chr <- as.vector(seqnames(selectGRange))
      strand <- as.vector(strand(selectGRange))
      start <- start(selectGRange)
      end <- end(selectGRange)

      ## which
      which.str <- paste("RangesList(", chr, "=IRanges(", start, ",", end,"))", sep="")

      ## flag
      if (strand=="+") flag.str <- "scanBamFlag(isMinusStrand=FALSE)"
      else flag.str <- "scanBamFlag(isMinusStrand=TRUE)"

      param.str <- paste("ScanBamParam(flag=", flag.str, ",which=", which.str ,")", sep="")
      return(param.str)
}

## extract selected data for some genes
selectGappedAlignments <- function(bamFile, selectGenes, gene.GRange)
{
    ## not needed in splicingEvent software. already have global geneGRange
    if (length(bamFile)==0) stop("bamFile has length 0!")

    geneName <- values(gene.GRange)[, "geneName"]

    for (i in 1:length(selectGenes))
    {
        gene.ind <- grep(selectGenes[i], geneName)
                          
        # get data for that gene
        selectGRange <- gene.GRange[gene.ind]
        if (length(selectGRange)==0)
        {
            print(paste(selectGenes[i], " does not match those in gtf file!", sep=""))
            return(NULL)
        }
        param.str <- convertParamStr(selectGRange)
        data.str <- paste("mydata <- readGAlignments(bamFile, param=", param.str, ")", sep="")

        eval(parse(text=data.str))
        if (i==1) all.data <- mydata
        else all.data <- c(all.data, mydata)
    }

    return(all.data)
}


## data double string + gap reads select
dataCheck <- function(mydata, gene.GRange)
{

  ## TODO: remove duplicates??

  #####################
  ## junction read filtering
  #####################
  # include trans+cis jun read, within read
  gene.read <- subsetByOverlaps(mydata, gene.GRange)
  if (length(gene.read)==0)
  {
      tmp.anno <- elementMetadata(gene.GRange)[, "geneName"]
      print(paste("No reads are available at the ", unique(tmp.anno)[1], sep=""))
      return(NULL) 
  }

  ## bug-fix: no inseration/deletion/..
  cigar.read <- cigar(gene.read)
  samOtherChars <- c("D","I","S","H", "P", "X", "=")
  useless.ind <- NULL
  sapply(1:length(samOtherChars), function(i)
        {
            tmp.ind <- grep(samOtherChars[i], cigar.read)
            useless.ind <<- c(useless.ind, tmp.ind)
            return(invisible(NULL))
        })
  if (length(useless.ind)!=0)
  {
      gene.read <- gene.read[-useless.ind]
      cigar.read <- cigar(gene.read)
  }
  
  skip.ind <- grep("N", cigar.read)
  if (length(skip.ind)==0) 
  {
      print("No junction reads are available in the bam file. Please check alignment!")
      return(NULL)
  }
   
  jun.read <- gene.read[skip.ind]
  nonjun.read <- gene.read[-skip.ind]
  
  ## jun read inside each genes, i.e., only cis-jun
  cis.read.raw <- subsetByOverlaps(jun.read, gene.GRange, type="within")

  ## make sure cis.read has at least two copys at the boundary of jun
  ds.tmp <- breakJunRead(cis.read.raw)
  cis.read <- ds.tmp$read

  left.jun.0 <- ds.tmp$left
  right.jun.0 <- ds.tmp$right
  middle.jun.0 <- ds.tmp$middle

  copy.str <- paste(chr=as.vector(seqnames(left.jun.0)),
                            pos=end(left.jun.0),
                            strand=as.vector(strand(left.jun.0)), 
                            pos=start(right.jun.0),
                            sep="_")
  dup.str <- copy.str[duplicated(copy.str)]

  copy.ind <- which(copy.str %in% dup.str)
  if (length(copy.ind)==0) 
  {
      print("No junction reads have exactly matching boundary! The number of reads may be NOT enough!")
      return(NULL)
  }

  cis.read <- cis.read[copy.ind]

  ## TODO: need the following two??
  left.jun <- left.jun.0[copy.ind]
  right.jun <- right.jun.0[copy.ind]
  middle.jun.0 <- middle.jun.0[copy.ind]
  middle.reduce <- reduce(middle.jun.0)
  
  return(list(left.jun=left.jun, right.jun=right.jun,
              middle.reduce=middle.reduce, middle.read=middle.jun.0,
              nonjun.read=nonjun.read, cis.read=cis.read))
              
}

## generate bed file for visualization
writeBed <- function(datatable, eventName, fileName)
{
  	
   	## only allow 4 colnames, if string, can use "append"
     colnames(datatable) <- c("track", paste("name=", eventName, sep=""),
   	                          "type=bed", "visibility=full")

     ## transfer to 0-based for output bed format
     datatable[,2] <- as.numeric(datatable[,2]) - 1 
     datatable[,3] <- as.numeric(datatable[,3]) - 1
     if (! file.exists("bedfiles") ) dir.create("bedfiles")
     fileName <- file.path("bedfiles", fileName)
     if (!file.exists(fileName))
        write.table(datatable, file=fileName, 
                  sep="\t", quote=FALSE, row.names=FALSE
                  )      
     else
     {
        write.table(datatable, file=fileName, 
                  sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE
                  )        
     }                              
   	
}

## Main algorithms for detecting RI: intron retention
riEvent <- function(middle.reduce, middle.jun.0, nonjun.read, 
                    minReadCounts=10, eventType=c("ri.1", "ri.2"),
                    sampleName, sampleID, total.read=0, gRange)
{                          
  if (missing(gRange)) 
  {
      gene.GRange <- splicingENV$gene.GRange
      exonGRange <- splicingENV$exonGRange
      intronGRange <- splicingENV$intronGRange
  }
  else
  {
      gene.GRange <- gRange$gene.GRange
      exonGRange <- gRange$exonGRange
      intronGRange <- gRange$intronGRange
  }
  ## NOTE: for small intron, some may extend outside of exon, hard to tell 
  ## e.g. Mlxipl
  #ir.nonjun.read <- subsetByOverlaps(nonjun.read, middle.reduce, type="within")
  ## Final decision, include these small exon extension
  ## e.g. Tubgcp6
  ir.nonjun.read <- subsetByOverlaps(nonjun.read, middle.reduce)           
  
  middle.c <- countOverlaps(middle.reduce, ir.nonjun.read)
  
  if (sum(middle.c >= minReadCounts)==0)
  {
      #print("No intron retition events!")
      return(NULL)
  }

      middle.candidate <- middle.reduce[middle.c >= minReadCounts]

      ## 2. should overlap with intron (NO assemble..based on current annotation)
      middle.c2 <- countOverlaps(middle.candidate, intronGRange)

      ## 3. only consider classic one with one intron overlap!!
      middle.final.1 <- middle.candidate[middle.c2 <= 1]
      middle.final.2 <- middle.candidate[middle.c2 > 1]
      
      # quantify the intron retention                 
      riCalculate <- function(middle.final, ir.nonjun.read, middle.jun.0, 
                         intronGRange, exonGRange, gene.GRange, minReadCounts, 
                         RIcoverage=TRUE, total.read=0
                         )
      {
          if (length(middle.final)==0)
          {
             #print("No intron retition events!")
             return(NULL)
          }
                                   
            ############## algorithm 1
            ## calculate read total length within each intron
            tmp.nonjun.reduce <- reduce(granges(ir.nonjun.read))
            tmp.m <- as.matrix(findOverlaps(middle.final, tmp.nonjun.reduce))
      
            tmp.nonjun.width <- width(tmp.nonjun.reduce)
            tmp.nonjun.df <- data.frame(ID=1:length(tmp.nonjun.reduce), width=tmp.nonjun.width)
      
            tmp.ind <- match(tmp.m[,2], tmp.nonjun.df$ID)
            tmp.df <- data.frame(tmp.m, nonwidth=tmp.nonjun.df[tmp.ind,]$width)
      
            tmp.readWidth <- tapply(tmp.df$nonwidth, tmp.df$queryHits, sum)

            tmp.ratio <- tmp.readWidth / width(middle.final)[as.numeric(names(tmp.readWidth))]

            RIrate <- rep(0, length(middle.final))
            RIrate[as.numeric(names(tmp.readWidth))] <- round(tmp.ratio,2)
            RIrate[RIrate > 1] <- 1
                        
            ## calculate jun link width, real intron width
            ir.nonjun.c <- countOverlaps(middle.final, ir.nonjun.read)
            ir.jun.c <- countOverlaps(middle.final, middle.jun.0)
            
            ir.width <- width(middle.final) ## not anno Intron length, but real intron length, reduce length
            ir.start <- start(middle.final) 
            ir.end <- end(middle.final) 

            ## convert to meaningful exon/intron
            tmp.list <- convertKnownName(middle.final, intronGRange, TRUE, "ri", gene.GRange)
            middle.final <- tmp.list$readGRange
            intron.start <- tmp.list$other$start
            intron.end <- tmp.list$other$end
               
            ## check start/end ?=0 , come from the previous function: convertKnownName
            zero.ind <- NULL
            if (sum(intron.start==0)>0)
            {
                  zero1.ind <- which(intron.start==0)
                  intron.start[zero1.ind] <- ir.start[zero1.ind]
                  zero.ind <- zero1.ind
                  
            }

            if (sum(intron.end==0)>0)
            {
                  zero2.ind <- which(intron.end==0)
                  intron.end[zero2.ind] <- ir.end[zero2.ind]
                  zero.ind <- unique(c(zero.ind, zero2.ind))
               
            }
                                   
            note <- rep("", length(ir.start))
            if (length(zero.ind)!=0) note[zero.ind] <- paste(note[zero.ind], "riInsideReducedExon", sep="")
                    
            ## calcuate intron width
            read.intron.m <- as.matrix(findOverlaps(middle.final, intronGRange))
            intron.width <- rep(0, length(middle.final)) 
                
            width.df <- data.frame(read.intron.m, width=width(intronGRange[read.intron.m[,2]]))                
            width.cal <- tapply(width.df$width, width.df$queryHits, sum)
            
            #width.ind <- match(c(1:length(middle.final)), names(width.cal))
            #intron.width[width.ind] <- width.cal[width.ind]
            
            intron.width[as.numeric(names(width.cal))] <- width.cal            
            intron.zeroind <- which(intron.width==0)
            intron.width[intron.zeroind] <- ir.width[intron.zeroind]

            ## width normalization factor, use intron.width(intron) or ir.width(read)??  
            bam.readLen <- width(ir.nonjun.read)[1] ##TODO: input
            #nonjun.iso.c <- round(ir.nonjun.c * bam.readLen /ir.width,4)
                        
            ## ir.width: intron width inferred from reads
            if (total.read!=0) 
            {                       
              scaleFac.nonjun <-(10^6)/total.read/ir.width
              scaleFac.jun <-(10^6)/total.read/(bam.readLen*2)
              nonjun.iso.norm <- round(ir.nonjun.c*scaleFac.nonjun, 4)
              jun.iso.norm <- round(ir.jun.c*scaleFac.jun, 4)                          
              score.all.ri <- nonjun.iso.norm/(nonjun.iso.norm + jun.iso.norm)
              
              intron.final <- data.frame(mergeID=paste(paste(as.vector(seqnames(middle.final)),
                                                       paste(ir.start,ir.end,sep="-"), sep=":"),
                                                       as.vector(strand(middle.final)),
                                                       sep="_"),      
                                      geneName=unlist(values(middle.final)[, "geneName"]),
                                      chr= as.vector(seqnames(middle.final)),   # chr need for bed file
                                      locus=paste(as.vector(seqnames(middle.final)),
                                            paste(intron.start,intron.end,sep="-"), sep=":"), 
                                      #strand=as.vector(strand(middle.final)),              
                                      nonjun=ir.nonjun.c, jun=ir.jun.c,
                                      nonjun.norm=nonjun.iso.norm,
                                      jun.norm=jun.iso.norm,
                                      #ratio=ir.ratio.norm,
                                      ratio=as.numeric(format(round(score.all.ri,4), scientific = FALSE)),
                                      #junLeftBorder=ir.start, junRightBorder=ir.end,
                                      realIntronStart=ir.start, realIntronEnd=ir.end,
                                      realIntronLocus=as.character(paste(as.vector(seqnames(middle.final)),
                                            paste(ir.start,ir.end,sep="-"), sep=":")), 
                                      intronLocus=paste(as.vector(seqnames(middle.final)),
                                            paste(intron.start,intron.end,sep="-"), sep=":"),
                                      note=note
                                      , RIpercent=RIrate, stringsAsFactors = FALSE                                      
                                      #, score=score.all.ri
                                      )
              }
              else
              {
                  scaleFac.nonjun <-(10^6)/ir.width
                  scaleFac.jun <-(10^6)/(bam.readLen*2)
                  nonjun.iso.norm <- round(ir.nonjun.c*scaleFac.nonjun, 4)
                  jun.iso.norm <- round(ir.jun.c*scaleFac.jun, 4)                                            
                  score.all.ri <- nonjun.iso.norm/(nonjun.iso.norm + jun.iso.norm)            

                  intron.final <- data.frame(mergeID=paste(paste(as.vector(seqnames(middle.final)),
                                                       paste(ir.start,ir.end,sep="-"), sep=":"),
                                                       as.vector(strand(middle.final)),
                                                       sep="_"),
                                      geneName=unlist(values(middle.final)[, "geneName"]),
                                      chr= as.vector(seqnames(middle.final)),   # chr need for bed file
                                      locus=paste(as.vector(seqnames(middle.final)),
                                            paste(intron.start,intron.end,sep="-"), sep=":"),
                                      #strand=as.vector(strand(middle.final)),              
                                      nonjun=ir.nonjun.c, jun=ir.jun.c,
                                      nonjun.norm=nonjun.iso.norm,
                                      jun.norm=jun.iso.norm,                                      
                                      #ratio=ir.ratio.norm,
                                      ratio=as.numeric(format(round(score.all.ri,4), scientific = FALSE)),
                                      #junLeftBorder=ir.start, junRightBorder=ir.end,
                                      realIntronStart=ir.start, realIntronEnd=ir.end,
                                      realIntronLocus=as.character(paste(as.vector(seqnames(middle.final)),
                                            paste(ir.start,ir.end,sep="-"), sep=":")), 
                                      intronLocus=paste(as.vector(seqnames(middle.final)),
                                            paste(intron.start,intron.end,sep="-"), sep=":"),
                                      note=note
                                      , RIpercent=RIrate, stringsAsFactors = FALSE
                                      #, RIcoverage=round(as.numeric(ir.width)/as.numeric(intron.width), 4)                                      
                                      #, score=score.all.ri
                                      )              
              }
   
            filter.ind <- which((intron.final$nonjun >= minReadCounts) & (intron.final$jun >= minReadCounts))
            if (length(filter.ind)>0) 
            {
                #filter.ind <- c(1:length(intron.final))
                intron.final <- intron.final[filter.ind,]
                datatable <- intron.final[, c("chr", "realIntronStart", "realIntronEnd")]

                  datatable$ratio <- paste(intron.final$ratio, "=",
                                         nonjun.iso.norm[filter.ind], "/(",
                                         nonjun.iso.norm[filter.ind], "+", jun.iso.norm[filter.ind], ")", 
                                         #intron.final$nonjun, "/(", intron.final$jun, 
                                         #"*", intron.final$junRightBorder-intron.final$junLeftBorder+1,
                                         #"+", intron.final$nonjun, ")", 
                                         #"_norm_", 
                                         "_", 
                                         sampleName, "_", eventType, ".", 
                                         intron.final$note, sep="")
                        
                                         
                intron.final$note2 <- datatable$ratio # same for intron.final
            
                eventType <- gsub("\\.", "_", eventType)
                writeBed(datatable, "intronRetention",                                
                                      paste(sampleName, "_", eventType, ".bed", sep="")
                                      )         
            }
            else return(NULL)
            
           ## annotate novel           
           #fstart <- intron.final$realIntronStart-1
           #fend <- intron.final$realIntronEnd+1
           #intron.final <- annotateNovel(intron.final, fstart, fend, exonGRange, etype="ri")
            
           return(intron.final)
     }

     if ( ("ri.1" %in% eventType) && ("ri.2" %in% eventType) )
     {             
      intron.1 <- riCalculate(middle.final.1, ir.nonjun.read, middle.jun.0, intronGRange, exonGRange, gene.GRange, minReadCounts,  RIcoverage=TRUE, total.read=total.read)
      intron.2 <- riCalculate(middle.final.2, ir.nonjun.read, middle.jun.0, intronGRange, exonGRange, gene.GRange, minReadCounts, RIcoverage=FALSE, total.read=total.read)
      intron.done <- list(type1=intron.1, type2=intron.2)
   
     }
     else if ("ri.2" %in% eventType) 
     {
      intron.2 <- riCalculate(middle.final.2, ir.nonjun.read, middle.jun.0, intronGRange, exonGRange, gene.GRange, minReadCounts,  total.read=total.read)
      intron.done <- list(type2=intron.2)
      
     }
     else
     {
      intron.1 <- riCalculate(middle.final.1, ir.nonjun.read, middle.jun.0, intronGRange, exonGRange, gene.GRange, minReadCounts, total.read=total.read)
      intron.done <- list(type1=intron.1)
     }

     return(intron.done) 
      

}


## Main algorithms for detecting ES: exon skipping
esEvent <- function(middle.reduce, middle.jun.0, nonjun.read, 
                    minReadCounts=10,
                    eventType=c("es.1", "es.2"), 
                    sampleName, sampleID, total.read=0, gRange 
                    )
{
  if (missing(gRange)) 
  {
      gene.GRange <- splicingENV$gene.GRange
      exonGRange <- splicingENV$exonGRange
      intronGRange <- splicingENV$intronGRange
  }
  else
  {
      gene.GRange <- gRange$gene.GRange
      exonGRange <- gRange$exonGRange
      intronGRange <- gRange$intronGRange
  }  
  
  
  esAlgorithm <- function(middle.reduce.2, middle.jun.0, minReadCounts)  
  {   
      middle.unique <- unique(middle.jun.0)
      t.c <- countOverlaps(middle.reduce.2, middle.unique)
      middle.reduce.3 <- middle.reduce.2[t.c>=3] # 3: itself, left jun, right jun
      middle.unique.2 <- subsetByOverlaps(middle.unique, middle.reduce.3, type="within")  

      t.c.2 <- countOverlaps(middle.unique.2, middle.unique.2)
      target.unique <- middle.unique.2[t.c.2>=3]
      
      t.c.3 <- countOverlaps(middle.unique.2, target.unique, type="equal")
      target.ind <- which(t.c.3!=0)
      nontarget.unique <- middle.unique.2[-target.ind]
      if (length(nontarget.unique)==0) return(NULL)
      
      m.start <- as.matrix(findOverlaps(target.unique, nontarget.unique, type="start"))
      m.end <- as.matrix(findOverlaps(target.unique, nontarget.unique, type="end"))
      if (nrow(m.start)==0 || nrow(m.end)==0) return(NULL)
      
      non.df <- data.frame(subjectHits=c(1:length(nontarget.unique)),
                           start=start(nontarget.unique),
                           end=end(nontarget.unique))
      
      m.start.df <- merge(m.start, non.df, by="subjectHits")
      m.end.df <- merge(m.end, non.df, by="subjectHits")
      
      t.start <- tapply(m.start.df$end, m.start.df$queryHits, min)
      t.end <- tapply(m.end.df$start, m.end.df$queryHits, max)

      t.start.df <- data.frame(ID=as.numeric(names(t.start)), start=t.start)
      t.end.df <- data.frame(ID=as.numeric(names(t.end)), end=t.end)

      t.done <- merge(t.start.df, t.end.df, by="ID")      
      t.done <- t.done[t.done$end > t.done$start,]
            
      ## filter with min reads
      final.middle.0 <- target.unique[as.numeric(t.done$ID)]
      
      final.c <- countOverlaps(final.middle.0, middle.jun.0, type="equal")
      final.ind <- which(final.c>=minReadCounts)
      
      final.middle <- final.middle.0[final.ind]
      final.df <- t.done[final.ind,]
 
      final.exon <- GRanges(seqnames=seqnames(final.middle),
                          strand=strand(final.middle),
                          ranges=IRanges(start=final.df$start+1, end=final.df$end-1
                          ))
      
      if (length(final.middle)==0) return(NULL)
      else
      {
          return(list(final.middle=final.middle, final.exon=final.exon))
      }
  }

  ## type split
  es.canList <-   esAlgorithm(middle.reduce, middle.jun.0, minReadCounts) 
  if (is.null(es.canList)) return(NULL)

  middle.1 <- es.canList$final.middle
  es.1 <- es.canList$final.exon 
  middle.c <- countOverlaps(es.1, exonGRange)
                         
  es.type1 <- es.1[middle.c<=1]
  es.type2 <- es.1[middle.c>1]
  
  middle.type1 <- middle.1[middle.c<=1]
  middle.type2 <- middle.1[middle.c>1]

  esCalculate <- function(middle.type1, es.type1, nonjun.read, middle.jun.0, exonGRange, gene.GRange,  
                    minReadCounts, eventType, total.read=0)  
  {      

      ## count exon
      es.c2 <- countOverlaps(es.type1, nonjun.read)
      es.candidate <- es.type1[es.c2 >= minReadCounts]
      middle.candidate <- middle.type1[es.c2 >= minReadCounts]
    
      if (length(es.candidate)==0)
      {
          #print("No Exon Skipping events.")
          return(NULL)
      }

      ############## algorithm 1

      es.nonjun.c <- countOverlaps(es.candidate, nonjun.read)

      # conver to jun middle one
      es.jun.c <- countOverlaps(middle.candidate, middle.jun.0, type="within")

      ##TODO
      es.width <- width(es.candidate) # width of exon inferred from reads
      es.start <- start(es.candidate)
      es.end <- end(es.candidate) 
            
      ## calculate
      bam.readLen <- width(nonjun.read)[1] ##TODO: input
      
      #es.ratio.raw <- es.jun.c/(es.nonjun.c + es.jun.c)
      #es.ratio.norm <- es.jun.c/(es.nonjun.c/es.width + es.jun.c)

      tmp.list <- convertKnownName(es.candidate, exonGRange, TRUE, "es", gene.GRange)
      es.candidate <- tmp.list$readGRange
      exon.start <- tmp.list$other$start
      exon.end <- tmp.list$other$end 
                      

      ## norm read or not
      if (total.read!=0) 
      {
          scaleFac.nonjun <-(10^6)/total.read/es.width
          scaleFac.jun <-(10^6)/total.read/(bam.readLen*2)
              
          nonjun.iso.norm <- round(es.nonjun.c * scaleFac.nonjun, 4)
          jun.iso.norm <- round(es.jun.c * scaleFac.jun, 4)
          
          score.all.es <- jun.iso.norm/(jun.iso.norm + nonjun.iso.norm)
          #scaleFac <-(10^6)/total.read/es.width
              
          es.final <- data.frame(mergeID=paste(paste(as.vector(seqnames(es.candidate)),
                                                 paste(es.start,es.end,sep="-"), sep=":"),
                                                 as.vector(strand(es.candidate)),
                                                 sep="_"),
                                  geneName=unlist(values(es.candidate)[, "geneName"]),
                                  chr=as.vector(seqnames(es.candidate)),  # chr need for bed file
                                  locus=paste(as.vector(seqnames(es.candidate)),
                                            paste(exon.start,exon.end,sep="-"), sep=":"), 
                                  #strand=as.vector(strand(es.candidate)),          
                                  nonjun=es.nonjun.c, jun=es.jun.c,
                                  nonjun.norm=nonjun.iso.norm,
                                  jun.norm=jun.iso.norm, 
                                  #ratio=es.ratio.norm,
                                  ratio= as.numeric(format(round(score.all.es,4), scientific = FALSE)),
                                  #strand=as.vector(strand(es.candidate)),
                                  #junLeftBorder=start(middle.candidate)-1, junRightBorder=end(middle.candidate)+1,
                                  realExonStart=es.start, realExonEnd=es.end,
                                  realExonLocus=as.character(paste(as.vector(seqnames(es.candidate)),
                                            paste(es.start,es.end,sep="-"), sep=":")), 
                                  realJunLocus=as.character(paste(as.vector(seqnames(middle.candidate)),
                                            paste(start(middle.candidate)-1, end(middle.candidate)+1,sep="-"), 
                                            sep=":")),
                                  exonLocus=paste(as.vector(seqnames(es.candidate)),
                                            paste(exon.start,exon.end,sep="-"), sep=":")                               
                                  , stringsAsFactors = FALSE
                                  #, score=score.all.es
                                  )
      }
      else
      {
          scaleFac.nonjun <-(10^6)/es.width
          scaleFac.jun <-(10^6)/(bam.readLen*2)
              
          nonjun.iso.norm <- round(es.nonjun.c * scaleFac.nonjun, 4)
          jun.iso.norm <- round(es.jun.c * scaleFac.jun, 4)          
          score.all.es <- jun.iso.norm/(jun.iso.norm + nonjun.iso.norm)
          
          es.final <- data.frame(mergeID=paste(paste(as.vector(seqnames(es.candidate)),
                                                 paste(es.start,es.end,sep="-"), sep=":"),
                                                 as.vector(strand(es.candidate)),
                                                 sep="_"),
                                  geneName=unlist(values(es.candidate)[, "geneName"]),
                                  chr=as.vector(seqnames(es.candidate)), # chr need for bed file
                                  locus=paste(as.vector(seqnames(es.candidate)),
                                            paste(exon.start,exon.end,sep="-"), sep=":"),
                                  #strand=as.vector(strand(es.candidate)),          
                                  nonjun=es.nonjun.c, jun=es.jun.c,
                                  nonjun.norm=nonjun.iso.norm,
                                  jun.norm=jun.iso.norm,                                   
                                  #ratio=es.ratio.norm,
                                  ratio=as.numeric(format(round(score.all.es,4), scientific = FALSE)),
                                  #strand=as.vector(strand(es.candidate)),
                                  #junLeftBorder=start(middle.candidate)-1, junRightBorder=end(middle.candidate)+1,
                                  realExonStart=es.start, realExonEnd=es.end, 
                                  realExonLocus=as.character(paste(as.vector(seqnames(es.candidate)),
                                            paste(es.start,es.end,sep="-"), sep=":")), 
                                  realJunLocus=as.character(paste(as.vector(seqnames(middle.candidate)),
                                            paste(start(middle.candidate)-1, end(middle.candidate)+1,sep="-"), 
                                            sep=":")),
                                  exonLocus=paste(as.vector(seqnames(es.candidate)),
                                            paste(exon.start,exon.end,sep="-"), sep=":")                                                        
                                  , stringsAsFactors = FALSE
                                  #, score=score.all.es
                                  )      
      }

      ## TODO filter
      es.filter.ind <- which((es.final$nonjun >= minReadCounts) & (es.final$jun >= minReadCounts))
      if (length(es.filter.ind)>0) 
      {
          #es.filter.ind <- c(1:length(es.final))
          es.final <- es.final[es.filter.ind,]
          
          dataTable <- es.final[, c("chr", "realExonStart", "realExonEnd")]
          dataTable$ratio <- paste(es.final$ratio, "=", 
                                         jun.iso.norm[es.filter.ind], "/(", 
                                         jun.iso.norm[es.filter.ind], 
                                         "+", nonjun.iso.norm[es.filter.ind], ")",
                                         "_", sampleName, "_", eventType, sep="")
          
          es.final$note2 <-  dataTable$ratio
          eventType <- gsub("\\.", "_", eventType)
          writeBed(dataTable, "exonSkipping", 
                         paste(sampleName, "_", eventType, ".bed", sep=""))       
      }
      else return(NULL)
 
      ## annotate novel           
      #fstart <- es.final$realExonStart
      #fend <- es.final$realExonEnd
      #es.final <- annotateNovel(es.final, fstart, fend, exonGRange, etype="es")
     
     return(es.final)
 }
 
     if ( ("es.1" %in% eventType) && ("es.2" %in% eventType) )
     {
      exon.1 <-   esCalculate(middle.type1, es.type1, nonjun.read, middle.jun.0, exonGRange, gene.GRange, 
                    minReadCounts, eventType, total.read=total.read)
      exon.2 <-   esCalculate(middle.type2, es.type2, nonjun.read, middle.jun.0, exonGRange, gene.GRange, 
                    minReadCounts, eventType, total.read=total.read)
      exon.done <- list(type1=exon.1, type2=exon.2)

     }
     else if ("es.2" %in% eventType) 
     {
      exon.2 <-   esCalculate(middle.type2, es.type2, nonjun.read, middle.jun.0, exonGRange, gene.GRange,  
                    minReadCounts, eventType, total.read=total.read)
      exon.done <- list(type2=exon.2)
      
     }
     else
     {
      exon.1 <-   esCalculate(middle.type1, es.type1, nonjun.read, middle.jun.0, exonGRange, gene.GRange, 
                    minReadCounts, eventType, total.read=total.read)
      exon.done <- list(type1=exon.1)
     
     }
      
     return(exon.done)
}


## simple link to main functions for as: alternative sites
asEvent <- function(cis.read.ad, left.jun, right.jun, 
          sampleName, sampleID, minReadCounts, eventType, novel, 
          total.read=0, gRange)
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

    all.reduce <- sort(c(reduce(exonGRange), reduce(intronGRange)))
    all <- sort(c(exonGRange, intronGRange))

    ## type I and II selection
    ## step 1: only work for type I selection, not for type II!

    ## NOTE: special case: some - typeI, some -type II, left same!

    if (eventType %in% c("adleft.1", "adleft.2"))
    {      
        left.donor <- adSingleSide(cis.read.ad, left.jun, right.jun, minReadCounts, adType="leftAlternative")
        left.donor.2 <- typeSelection(left.donor, eventType, intronGRange)
        ratio.ad.left <- calcRatio(left.jun, right.jun, left.donor.2, "leftAlternative", minReadCounts, gRange, novel)
        asTable <- writeADEvent(ratio.ad.left, eventType, all, all.reduce, sampleID, sampleName, total.read=total.read)
    }
    else if (eventType %in% c("adright.1", "adright.2"))
    {
        right.donor <- adSingleSide(cis.read.ad, left.jun, right.jun, minReadCounts, adType="rightAlternative")
        right.donor.2 <- typeSelection(right.donor, eventType, intronGRange)
        ratio.ad.right <- calcRatio(left.jun, right.jun, right.donor.2, "rightAlternative", minReadCounts, gRange, novel)
        asTable <- writeADEvent(ratio.ad.right, eventType, all, all.reduce, sampleID, sampleName, total.read=total.read)
    }
    else if ( eventType %in% c("adboth.1", "adboth.2"))
    {     
        both.donor <- adbothSide(cis.read.ad, left.jun, right.jun, minReadCounts)
        both.donor <- typeSelection(both.donor, eventType, intronGRange)
        ratio.ad.both <- calcRatioBoth(left.jun, right.jun, both.donor, minReadCounts, gRange, novel)
        asTable <- writeADEvent(ratio.ad.both, eventType, all, all.reduce, sampleID, sampleName, total.read=total.read)
    }
    else stop("eventType is WRONG!")



    return(asTable)
}



