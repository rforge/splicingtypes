##################################
##  SplicingTypesAnno
##  - splicingGene.r: gene-level annotation
##################################

## combine as events for all samples
combineGene <- function(sampleList, splicingEvents=c("ri.1"), selectGene, all=TRUE)
{


    if (!(splicingEvents %in% c("ri.1", "ri.2", "es.1", "es.2", 
          "adleft.1", "adleft.2", "adright.1", "adright.2", "adboth.1", "adboth.2"))) stop("splicingEvents Must be ri.1, ri.2, es.1, es.2, adleft.1, adleft.2, adright.1, adright.2, adboth.1, or adboth.2!")

        iintronTable <- NULL
        for (i in 1:length(sampleList))
        {
            if (is.null(iintronTable)) 
            {
                add.table <- sampleList[[i]][[splicingEvents]]                        
                iintronTable <- add.table
            }
            else if (!is.null(sampleList[[i]][[splicingEvents]]))
            {               

                  add.table <- sampleList[[i]][[splicingEvents]]                               
                  iintronTable <- merge(data.frame(iintronTable),
                                add.table, by="mergeID", all=all)

            }

        }

  return(iintronTable)
}

## count exon_intron reads
splicingCount <- function(selectGene,
                      bamFile = "accepted_hits.sort.bam",
                      select.list, 
                      sampleName="sample",
                      sampleID=1, type="any"
                      )
{
  if (length(selectGene) != 1) stop("Please select one gene!")
  if (!file.exists(bamFile)) stop("bamFile does NOT exist!")

  exonGRange <- select.list$exon
  gene.GRange <- select.list$gene
  intronGRange <- select.list$intron
  
  if (length(gene.GRange)==0 || length(exonGRange)==0) stop("select.list has 0 gene or exon!")
  if (length(intronGRange)==0) print("select.GRane has 0 intron!")

  ## preprocess GRange
     geneName <- values(gene.GRange)[, "geneName"]
     gene.ind <- grep(selectGene, geneName)
               
     geneName <- values(exonGRange)[, "geneName"]
     exon.ind <- grep(selectGene, geneName)
               
     geneName <- values(intronGRange)[, "geneName"]
     intron.ind <- grep(selectGene, geneName)             
            
      gene.GRange <- gene.GRange[gene.ind]
      exonGRange <- exonGRange[exon.ind]
      intronGRange <- setdiff(gene.GRange, reduce(exonGRange))
      
      ## intron annotation 
      if (length(exonGRange) > 0) elementMetadata(exonGRange)[, "geneName"] <- paste(elementMetadata(exonGRange)[, "geneName"], "exon", 1:length(exonGRange), sep="_")      
      if (length(intronGRange) > 0) elementMetadata(intronGRange)[, "geneName"] <- paste(elementMetadata(gene.GRange)[, "geneName"], "intron", 1:length(intronGRange), sep="_")

      ## read in data for selected genes
      print(paste("GENE - ", selectGene, " are processed...", sep=""))
      mydata <- selectGappedAlignments(bamFile=as.character(bamFile),
                         selectGene, gene.GRange
                         )
      if (is.null(mydata))  return(NULL)
      
      if (any(!(seqlevels(mydata) %in% seqlevels(gene.GRange))))
      {
          print("some chrNo in bam file are NOT in gtf/gff3 file! The unmatched data will be overlooked...")
      }
      
      if (any(!(seqlevels(gene.GRange) %in% seqlevels(mydata)) ))
      {
          print("some chrNo in gtf/gff3 file are NOT in bam file! The unmatched annotation will be overlooked...")
          seqlevels(gene.GRange, force=TRUE) <- seqlevels(mydata)
          seqlevels(exonGRange, force=TRUE) <- seqlevels(mydata)
          seqlevels(intronGRange, force=TRUE) <- seqlevels(mydata)
      }
               
      # total counts
      total.count <- countOverlaps(gene.GRange, mydata, type=type)
      names(total.count) <- values(gene.GRange)[, "geneName"]
      exon.count <- countOverlaps(exonGRange, mydata, type=type)
      names(exon.count) <- values(exonGRange)[, "geneName"]
      intron.count <- countOverlaps(intronGRange, mydata, type=type)
      names(intron.count) <- values(intronGRange)[, "geneName"]
      return(list(totalCount=total.count, exonCount=exon.count, intronCount=intron.count))
      
}

# combine counts for multiple samples
combineCount <- function(sample.list)
{
   if (length(sample.list) == 0) stop("No sample available!")
   
   total.count <- sample.list[[1]]$totalCount
   exon.count <- sample.list[[1]]$exonCount
   intron.count <- sample.list[[1]]$intronCount
     
   if (length(sample.list)==1)
   {
       return(list(totalCount=total.count, exonCount=exon.count, intronCount=intron.count))
   }
   
   ## more than 1 sample
   exon.matrix <- matrix(0, nrow=length(exon.count), ncol=length(sample.list))
   exon.matrix[,1] <- exon.count
   rownames(exon.matrix) <- names(exon.count) 
   intron.matrix <- matrix(0, nrow=length(intron.count), ncol=length(sample.list))
   intron.matrix[,1] <- intron.count
   rownames(intron.matrix) <- names(intron.count)
   
   for (i in 2:length(sample.list))
   {
      total.count <- c(total.count, sample.list[[i]]$totalCount) 
      exon.matrix[,i] <- sample.list[[i]]$exonCount
      intron.matrix[,i] <- sample.list[[i]]$intronCount
   }
  
   return(list(totalCount=total.count, exonCount=exon.matrix, intronCount=intron.matrix))
   
}


