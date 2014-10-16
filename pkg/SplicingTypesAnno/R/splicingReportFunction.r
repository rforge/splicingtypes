########################################
##  SplicingTypesAnno
##  - function.r: collection of main functions              
########################################

# clean environment
initEnv <- function()
{
  assign("table.list", list(), splicingENV)
  assign("spData", list(), splicingENV)
  
  assign("gene.GRange", NULL, splicingENV)
  assign("exonGRange", NULL, splicingENV)
  assign("intronGRange", NULL, splicingENV)
  assign("splicingLink", NULL, splicingENV)
  
  assign("dataSum", 
        list(ri.1=NULL, ri.2=NULL, es.1=NULL, es.2=NULL, 
        adleft.1=NULL, adleft.2=NULL, adright.1=NULL, adright.2=NULL, 
        adboth.1=NULL, adboth.2=NULL), 
        splicingENV)
  
  #assign("geneSum", 
  #      list(nonOverlap=NULL, overlap=NULL), 
  #      splicingENV)
  
  assign("readSum", 
        list(totalRead=list(), 
             qualityJunRead=list(),
             qualityNonJunRead=list(),
             sampleLevel=NULL
             ), 
        splicingENV)

}

# clean exisiting bed file to clear confusion
cleanBedFiles <- function()
{                   
    bedpath <- file.path(splicingENV$bed.dir)
    if (file.exists(bedpath)) unlink(bedpath, recursive=TRUE)
    
}

# combine splicingTypes for current sample
#  -rbind NOT merge 
addGene <- function(sampleNow, tmp.list, eventsType)
{
    for ( i in 1:length(eventsType))
    {
        if (is.null(splicingENV$table.list[[sampleNow]][[eventsType[i]]]))  splicingENV$table.list[[sampleNow]][[eventsType[i]]] <- tmp.list[[eventsType[i]]]
        else  splicingENV$table.list[[sampleNow]][[eventsType[i]]] <- rbind(splicingENV$table.list[[sampleNow]][[eventsType[i]]] ,tmp.list[[eventsType[i]]])
    }
    return(invisible(NULL))    
}

# generate exon_intron GRange object for selectGenes
selectGRange <- function(selectGenes)
{

      ## get gene range
      gene.ind <- NULL
      exon.ind <- NULL
      intron.ind <- NULL
      sapply(1:length(selectGenes), function(i)
            {
               geneName <- values(splicingENV$gene.GRange)[, "geneName"]
               gene.ind <<- c(gene.ind, grep(selectGenes[i], geneName))
               
               geneName <- values(splicingENV$exonGRange)[, "geneName"]
               exon.ind <<- c(exon.ind, grep(selectGenes[i], geneName))
               
               geneName <- values(splicingENV$intronGRange)[, "geneName"]
               intron.ind <<- c(intron.ind, grep(selectGenes[i], geneName))             
               return(invisible(NULL))
            })
      splicingENV$gene.GRange <- splicingENV$gene.GRange[gene.ind]
      splicingENV$exonGRange <- splicingENV$exonGRange[exon.ind]
      
      ## TODO: more clear for mutlipe gene matching...
      splicingENV$intronGRange <- splicingENV$intronGRange[exon.ind]
}

# Input data configure
# - bamfile + sampleName required; sampleID optional
inputSampleFile <- function(sampleFile=NULL, sampleList=NULL)
{       
    # read sample data
    if (!missing(sampleFile))
    {
        if (! file.exists(sampleFile) )  stop("Sample worksheet does not exist!")
        sample <- read.delim(sampleFile, header=TRUE)
        #sample <- read.xls(sampleFile, sheet=1, skip=1)    
        ## TODO
        #raw.ind <- grep("RawFiles", colnames(sample))
        #colnames(sample)[raw.ind] <- "BamFiles"
        
        if (!all(c("SampleName", "BamFiles") %in% colnames(sample) )) stop("The input file must have column names: SampleName and BamFiles!")
        
        if (!("SampleID" %in% colnames(sample))) sample$SampleID <- c(1:nrow(sample))
        
    }
    else 
    { 
        if (length(sampleList$SampleName)==0) stop("Please provide SampleName in the InputSData!")
        else SampleName= sampleList$SampleName
        
        if (length(sampleList$BamFiles)==0 ) stop("Please provide BamFiles in the InputSData!")
        else BamFiles=sampleList$BamFiles 
        
        if (length(sampleList$SampleName) != length(sampleList$BamFiles) ) stop("SampleName must have same number as BamFiles!")
        
        SampleID <- c(1:length(sampleList$SampleName))
        
        sample <- data.frame(SampleName=as.character(SampleName), 
                              BamFiles=as.character(BamFiles),
                              SampleID=SampleID)
    }
    
    ## clean data
    
    name.noind <- which(sample$SampleName=="")
    bam.noind <- which(sample$BamFiles=="")
    if (length(name.noind) !=0 || length(bam.noind) !=0)
    {
        if (length(name.noind)==length(bam.noind) & all(name.noind %in% bam.noind))
        {
            sample <- sample[-name.noind,]
        }
        else
        {
            stop("Input file contains empty string! Please clean the input file!")
        }
    }
    
    splicingENV$spData <- list(pData=sample)

}

# generate sample section
generateHTML.Raw <- function(raw.dir = raw.dir, outputRawHTML=outputRawHTML)
{   
   html.raw <- paste(
      "<table><tbody><tr class=\"head\"> <td
      		 colspan=\"2\">
      Sample information</td></tr>
      <tr class=\"line1\"><td colspan=\"2\">
      	The analysis is based on the following information. </td></tr>
     <tr class=\"line2\"><td width=\"50%\">Annotation of the <b>samples</b> for which sequence was attempted</td><td width=\"50%\">&nbsp;<a href=\"", raw.dir, "/", outputRawHTML,"\">html</a></td></tr>
     </td></tr></tbody></table><br><br> "
      , sep="")

    sample <- splicingENV$spData$pData
    #sample <- sample[,-1]

    if (! file.exists(raw.dir) ) dir.create(raw.dir)
		p=openPage(file.path(raw.dir,outputRawHTML), title="Sample information", link.css="http://www.ebi.ac.uk/~gpau/hwriter/hwriter.css")
		hwrite("Annotation of the samples" , p, heading=3,center=TRUE)
		hwrite(sample,p,row.names=TRUE,row.bgcolor='#ddaaff',row.style=list('font-weight:bold'),br=TRUE)
		closePage(p)

		print("Generate sample annotation...Done")
		return(html.raw)
      
}

# generate gene overlapping stuats section (TODO: need it?)
generateHTML.Gene <- function(gene.dir, outputGeneHTML, sampleName, sampleID)
{   
 
  ## summary for the gene overlapping
  gene.c <- countOverlaps(splicingENV$gene.GRange, splicingENV$gene.GRange)
  non.no <- sum(gene.c <= 1)
  overlap.no <- sum(gene.c > 1)
  
  geneName.non <- elementMetadata(splicingENV$gene.GRange[gene.c <= 1])[, "ID"]
  geneName.over <- elementMetadata(splicingENV$gene.GRange[gene.c > 1])[, "ID"]  

    if (! file.exists(gene.dir) ) dir.create(gene.dir)
    if (length(geneName.non) > 0) write.csv(geneName.non, 
                                      file=file.path(paste(gene.dir, "geneNonOverlap.csv", sep="/")) )
    if (length(geneName.over) > 0) write.csv(geneName.over,
                                      file=file.path(paste(gene.dir, "geneOverlap.csv", sep="/")) )                                      

   html.summary <- paste(
      "<table><tbody><tr class=\"head2\"> <td
      		 colspan=\"2\">
      Gene summary from GTF/GFF3 file</td></tr>
      <tr class=\"line1\"><td colspan=\"2\">
      	The GTF/GFF3 has the following genes: </td></tr>
     <tr class=\"line2\"><td width=\"50%\">Genes without overlapping</td><td width=\"50%\">", non.no,
      ",&nbsp;<a href=\"", file.path(paste(gene.dir, "geneNonOverlap.csv", sep="/")),
      "\">csv</a></td></tr>
     <tr class=\"line2\"><td width=\"50%\">Genes with overlapping</td><td width=\"50%\">", overlap.no, 
      ",&nbsp;<a href=\"", file.path(paste(gene.dir, "geneOverlap.csv", sep="/")),
      "\">csv</a></td></tr>
      </tbody></table><br><br> ",
      sep="")

		print("Generate gene section from GTF/GFF3...Done")
		return(html.summary)
      
}

# generate read number summary, including total read, quality read, etc.
generateHTML.Read <- function(read.dir, outputReadHTML, sampleName, sampleID, selectGenes)
{   

    if (! file.exists(read.dir) ) dir.create(read.dir)

   html.summary <- paste(
      "<table><tbody><tr class=\"head2\"> <td
      		 colspan=\"2\">
      Read summary</td></tr>
      <tr class=\"line1\"><td colspan=\"2\">
      	The reads information table: </td></tr>
     <tr class=\"line2\"><td width=\"50%\">reads table</td>
     <td width=\"50%\"><a href=\"", 
     file.path(read.dir,outputReadHTML),
      "\">html</a></td></tr>
      </tbody></table><br><br> ",
      sep="")
          
      ## read table
    myread.table <- data.frame(SampleName=sampleName, SampleID=sampleID)
    
    if (length(splicingENV$readSum$totalRead)==0) 
    {
        return(invisible(NULL))
    }
             
    if (is.null(selectGenes))
    {
        myread.table$totalRead <- splicingENV$readSum$totalRead[[1]]
        myread.table$qualityJunRead <- splicingENV$readSum$qualityJunRead[[1]]
        myread.table$qualityNonJunRead <- splicingENV$readSum$qualityNonJunRead[[1]]

    }
    else
    {           
        for (i in 1:length(selectGenes))
        {
            tmp.df <- data.frame(totalRead=splicingENV$readSum$totalRead[[i]], 
                                 qualityJunRead=splicingENV$readSum$qualityJunRead[[i]],
                                 qualityNonJunRead=splicingENV$readSum$qualityNonJunRead[[i]]
                                 )
            colnames(tmp.df) <- paste(selectGenes[i], colnames(tmp.df), sep=".")
            myread.table <- cbind(myread.table, tmp.df)
                                 
        }
    }		
    p=openPage(file.path(read.dir,outputReadHTML), title="read information", link.css="http://www.ebi.ac.uk/~gpau/hwriter/hwriter.css")
		hwrite("reads summary" , p, heading=3,center=TRUE)
		hwrite(myread.table,p,row.names=TRUE,row.bgcolor='#ddaaff',row.style=list('font-weight:bold'),br=TRUE)
		closePage(p)

		print("Generate reads section...Done")
		return(html.summary)
      
}

# generate splicingType sumary
generateHTML.splicingSummary <- function(summary.dir, outputSummaryHTML, sampleName, sampleID)
{   
   html.summary <- paste(
      "<table><tbody><tr class=\"head2\"> <td
      		 colspan=\"2\">
      Splicing summary</td></tr>
      <tr class=\"line1\"><td colspan=\"2\">
      	The analysis found the following events. </td></tr>
     <tr class=\"line2\"><td width=\"50%\">Summary for splicing events</td><td width=\"50%\">&nbsp;<a href=\"", summary.dir, "/", outputSummaryHTML,"\">html</a></td></tr>
     </td></tr></tbody></table><br><br> "
      , sep="")


    splice.table <- data.frame(SampleName=sampleName, SampleID=sampleID)
    if (!is.null(splicingENV$dataSum$ri.1)) splice.table$intronRetention.Type1 <- splicingENV$dataSum$ri.1
    if (!is.null(splicingENV$dataSum$ri.2)) splice.table$intronRetention.Type2 <- splicingENV$dataSum$ri.2
    if (!is.null(splicingENV$dataSum$es.1)) splice.table$exonSkipping.Type1 <- splicingENV$dataSum$es.1 
    if (!is.null(splicingENV$dataSum$es.2)) splice.table$exonSkipping.Type2 <- splicingENV$dataSum$es.2  
    if (!is.null(splicingENV$dataSum$adleft.1)) splice.table$alternativeDonor.Type1 <- splicingENV$dataSum$adleft.1 
    if (!is.null(splicingENV$dataSum$adleft.2)) splice.table$alternativeDonor.Type2 <- splicingENV$dataSum$adleft.2 
    if (!is.null(splicingENV$dataSum$adright.1)) splice.table$alternativeAcceptor.Type1 <- splicingENV$dataSum$adright.1
    if (!is.null(splicingENV$dataSum$adright.2)) splice.table$alternativeAcceptor.Type2 <- splicingENV$dataSum$adright.2
    if (!is.null(splicingENV$dataSum$adboth.1)) splice.table$alternativeBothSites.Type1 <- splicingENV$dataSum$adboth.1                                     
    if (!is.null(splicingENV$dataSum$adboth.2)) splice.table$alternativeBothSites.Type2 <- splicingENV$dataSum$adboth.2    
    
    if (! file.exists(summary.dir) ) dir.create(summary.dir)
		p=openPage(file.path(summary.dir,outputSummaryHTML), title="Summary information", link.css="http://www.ebi.ac.uk/~gpau/hwriter/hwriter.css")
		hwrite("Summary for splicing events" , p, heading=3,center=TRUE)
		hwrite(splice.table,p,row.names=TRUE,row.bgcolor='#ddaaff',row.style=list('font-weight:bold'),br=TRUE)
		closePage(p)

		print("Generate summary section...Done")
		return(html.summary)
      
}

## combine bed files
generateHTML.Bed <- function(sampleName, eventType, outputDir)
{
   bedpath <- file.path(splicingENV$bed.dir)

   bedfiles <- dir(bedpath, pattern=".bed")
   
   targetfiles <- NULL
   eventType <- gsub("\\.", "_", eventType)
   for (i in 1:length(sampleName))
   {
      targetfiles <- c(targetfiles, 
                       paste(sampleName[i], "_", eventType, ".bed", sep="")
                       )
   }
   bedfiles <- bedfiles[ bedfiles %in% targetfiles ]
   if (length(bedfiles)==0) 
   {
		  print("NO SPLICING EVENTS! No bed section...Done")
		  return(NULL)
   }
   bedfiles <- file.path(splicingENV$bed.dir, bedfiles)   
   bed.table <- NULL
   
   for (i in 1:length(bedfiles))
   {
      tmp <- read.table(bedfiles[i], header=TRUE, sep="\t")
      colnames(tmp) <- c("track", "name=splicingEvent",
   	                          "type=bed", "visibility=full")
      bed.table <- rbind(bed.table, tmp)
      bed.table <- bed.table[order(bed.table[,1], bed.table[,2]),]
   }

    write.table(bed.table, file=file.path(splicingENV$bed.dir, "combined.bed"), sep="\t", quote=FALSE, 
                 row.names=FALSE)

   
   bedstr <- unlist(strsplit(bedfiles, "\\."))
   bedstr <- bedstr[seq(1, length(bedstr), by=2)]
   if (length(bedstr) != length(bedfiles)) stop("The bed names are not regular ones!")
   
   linkstr <- paste("<a href=\"", bedfiles, "\">", bedstr,"</a><br>", collapse=" ")

   html.bed <- paste(
      "<table><tbody><tr class=\"head2\"> <td
      		 colspan=\"2\">
      IGV visualization - bed file</td></tr>
      <tr class=\"line1\"><td colspan=\"2\">
      	The visualization files includes:. </td></tr>
     <tr class=\"line2\"><td width=\"50%\"><a href=\"bedfiles/combined.bed\">Combined bed file</a>",
     "</td><td width=\"50%\">Separate files&nbsp;<br>", linkstr,"</td></tr>
     </td></tr></tbody></table><br><br> "
      , sep="")

		print("Generate bed section...Done")
		return(html.bed)
}


## annotate splicing types at gene-level
generateHTML.annotateGene <- function(intron.dir, 
                            bamFile,
                            select.GRange,
                            eventsType,
                            minReadCounts=minReadCounts,
                            selectGenes=selectGenes, 
                            novel=novel,
                            parallel, cpus                        
                            )
{
    if (! file.exists(intron.dir) ) dir.create(intron.dir)
          
    tmp.list <- list()
    for (i in 1: length(bamFile))
    {
             
      for (j in 1:length(selectGenes))
      {
          oneBamFile <- bamFile[i]
          print(paste("Generate intron section:", oneBamFile,"...START......", sep="") )
                        
          # annotate at gene-level for single gene
          tmp.list <- 
          splicingGene(selectGenes[j],
                        bamFile = bamFile[i],
                        select.GRange,
                        eventsType,
                        sampleName=splicingENV$spData$pData$SampleName[i],
                        sampleID=splicingENV$spData$pData$SampleID[i], 
                        reportSummary=TRUE,
                        minReadCounts = minReadCounts, novel=novel
                        )
            ## combine all genes
            if ( (length(splicingENV$table.list) < i) )
            {
                splicingENV$table.list[[i]] <- tmp.list$type
            }
            else
            {
                addGene(i, tmp.list$type, eventsType)
    
            }        
        
        }
                

      }

            print("Generate intron section...end......")                 
}                            

## annotate splicing types at sample level
#   - parallel computing configure
generateHTML.annotateSample <- function(intron.dir, bamFile=as.character(splicingENV$spData$pData$BamFiles), 
                        gtfFile, eventsType, 
                        minReadCounts=minReadCounts, selectGenes=selectGenes, 
                        novel,
                        parallel=FALSE, cpus=2)
{                        
    two.list <- list()
        # export
    gtf.GRange <- list(gene.GRange=splicingENV$gene.GRange,
                       exonGRange=splicingENV$exonGRange,
                       intronGRange=splicingENV$intronGRange,
                       exonStart=splicingENV$exonStart,
                       exonEnd=splicingENV$exonEnd, 
                       splicingLink=splicingENV$splicingLink    
                       )
    #spList <- list(spData=splicingENV$spData, bed.dir=splicingENV$bed.dir)
                  
    if (! file.exists(intron.dir) ) dir.create(intron.dir)
    if (parallel)
    {          
        require(snowfall)
                                    
        sfInit(parallel=TRUE, cpus=cpus)
        #on.exit(sfStop())
        sfLibrary(SplicingTypesAnno)
              
        spData=splicingENV$spData                                        
        #sfExport(list=c("spData"))
        sfExportAll(except=NULL)

         two.list <- sfLapply(1:length(bamFile), function(i)
                            {

          #library(SplicingTypesAnno)
                  #suppressPackageStartupMessages({
                  #              require(SplicingTypesAnno, quietly = T)
                        #})
                         updateTRUEIntron(as.character(bamFile[i]),
                                   gtf.GRange,
                                   spData$pData$SampleName[i],
                                   spData$pData$SampleID[i],
                                   eventsType,
                                   minReadCounts=minReadCounts,
                                   selectGenes=selectGenes,
                                   novel=novel
                                   )
                            })


        sfStop()
        
        print("all parallel done!")
    }
    else
    {
        for (i in 1: length(bamFile))
        {
            oneBamFile <- bamFile[i]
            print(paste("Generate intron section:", oneBamFile,"...START......", sep="") )
                               
            ## update true intron signal, later need to combine with genIntron to save process time
             two.list[[i]] <-  updateTRUEIntron(bamFile = as.character(oneBamFile),
                                   gtf.GRange, 
                                   splicingENV$spData$pData$SampleName[i],
                                   splicingENV$spData$pData$SampleID[i], eventsType,
                                   minReadCounts=minReadCounts,
                                   selectGenes=selectGenes,
                                   novel=novel
                                   )

            print("pass updateTRUEIntron")

            print("Generate intron section...end......")

         }

         
    }          
    
        splicingENV$table.list[[1]] <- two.list[[1]]$type
        splicingENV$readSum$totalRead[[1]] <- two.list[[1]]$summary$totalRead
        splicingENV$readSum$qualityJunRead[[1]] <- two.list[[1]]$summary$qualityJunRead
        splicingENV$readSum$qualityNonJunRead[[1]] <- two.list[[1]]$summary$qualityNonJunRead
                
        ## i: samples
        ## [[1]]: only one virtual gene (global sample)
        if (length(two.list)>1)
        {
          for (i in 2:length(two.list))
          {
            splicingENV$table.list[[i]] <- two.list[[i]]$type
            splicingENV$readSum$totalRead[[1]] <- c(splicingENV$readSum$totalRead[[1]], 
                                                    two.list[[i]]$summary$totalRead)
            splicingENV$readSum$qualityJunRead[[1]] <- c(splicingENV$readSum$qualityJunRead[[1]],
                                                         two.list[[i]]$summary$qualityJunRead)
            splicingENV$readSum$qualityNonJunRead[[1]] <- c(splicingENV$readSum$qualityNonJunRead[[1]],
                                                            two.list[[i]]$summary$qualityNonJunRead)
                      
          }
        }                               

}

# convert gtf
convertGTFvalueToName <- function(tmp.name)
{
  tmp.name <- as.character(tmp.name)
  tmp <- unlist(strsplit(tmp.name, ";"))
  read.ind <- grep("ID=", tmp)

  if (length(read.ind)==0) read.ind <- grep("Parent=", tmp)
  if (length(read.ind)==0) return(invisible(NULL))
  read.tmp <- tmp[read.ind]

  read.name <- gsub("ID=", "", read.tmp)
  read.name <- gsub("Parent=", "", read.name)
  return(read.name)

}

# combine specific splicing types for multiple samples
generateHTML.intronTable <- function(intron.dir, gtfFile,
                            tableList,
                            splicingEvents="ri.1",
                            feature
                            )
{                    
    print(paste("Generate intron TABLE section...START......", splicingEvents, sep=""))
    html <- NULL
    if (! file.exists(intron.dir) ) dir.create(intron.dir)

    if (!(splicingEvents %in% c("ri.1", "ri.2", "es.1", "es.2", 
        "adleft.1", "adleft.2", "adright.1", "adright.2", "adboth.1", "adboth.2"))) stop("splicingEvents Must be ri.1, ri.2, es.1, es.2, adleft.1, adleft.2, adright.1, adright.2, or adboth.1, adboth.2!")
    keyhtml <- switch(splicingEvents,
            ri.1="Intron retention - type 1",
            ri.2="Intron retention - type 2",
            es.1="Exon skipping - type 1",
            es.2="Exon skipping - type 2",
            adleft.1="Alternative donor_type 1",
            adleft.2="Alternative donor_type 2",
            adright.1="Alternative acceptor_type 1",
            adright.2="Alternative acceptor_type 2",
            adboth.1="Alternative donor_both 1",
            adboth.2="Alternative donor_both 2"
            )
    iintronTable <- NULL
     
    allID <- NULL
    add.col <- 0
    event.col <- 0  ## for ratio table
    if (length(tableList) >= 1)
    {              
        ## combine all geneNames to have a big table        
        for (i in 1:length(tableList))
        {                                              
            if (!is.null(tableList[[i]][[splicingEvents]]))
            {
                
                add.table <- tableList[[i]][[splicingEvents]]
                add.col <- ncol(add.table)## take off chr, locus, name

                 allID <- unique(c(allID, as.character(add.table$mergeID)))
                 event.col <-  event.col + 1
                                     
            }
        }
                  
        if (event.col==0)
        {
            return(invisible(NULL))
        }

        ## ratio table
        all.df <- data.frame(mergeID=as.character(allID),
                             geneName=as.character(rep("geneName", length(allID))), 
                            matrix(0, nrow=length(allID), ncol=event.col),
                            stringsAsFactors = FALSE
                            )
                            
        ## note table
        note2.df <- data.frame(matrix(0, nrow=length(allID), ncol=event.col),
                            stringsAsFactors = FALSE)

        ## comprehensive table            
        all.csv <- data.frame(mergeID=as.character(allID), 
                              #chr=rep("chr", length(allID)),
                              geneName=rep("geneName", length(allID)),
                               data.frame(matrix(0, 
                                nrow=length(allID), 
                                ncol=event.col*(add.col-4))), ## don't need the first two col - share; the third one is not used -take off
                                stringsAsFactors = FALSE
                              )
             
        dataframe.col <- 0                   
        for (i in 1:length(tableList))
        {          
            splicingENV$dataSum[[splicingEvents]] <- c(splicingENV$dataSum[[splicingEvents]],
                                        ifelse(is.null(tableList[[i]][[splicingEvents]]), 0, nrow(tableList[[i]][[splicingEvents]]))
                                                  )
            if (!is.null(tableList[[i]][[splicingEvents]]))
            {              
                add.table <- tableList[[i]][[splicingEvents]]
                match.ind <- match(add.table$mergeID, all.df$mergeID)
                     
                note2.ind <- grep("note2", colnames(add.table))
                note2.df[match.ind, dataframe.col+1] <- add.table[,note2.ind]
                colnames(note2.df)[dataframe.col+1] <- colnames(add.table)[note2.ind]

                ratio.ind <- grep("ratio", colnames(add.table))
                add.table[,ratio.ind] <- as.character(add.table[,ratio.ind])
                # first col is gene name, second col is locus
                all.df[match.ind,2] <- as.character(add.table[,2])
                all.df[match.ind,dataframe.col+3] <- add.table[,ratio.ind]
                colnames(all.df)[dataframe.col+3] <- colnames(add.table)[ratio.ind]  

                col.ind <- c((dataframe.col*(add.col-4)+3):(dataframe.col*(add.col-4)+add.col-2))            
                

                all.csv[match.ind, 2] <- as.character(add.table[, 2])
                all.csv[match.ind, col.ind] <- add.table[, -c(1:4)]

                colnames(all.csv)[col.ind] <- colnames(add.table)[-c(1:4)]  
                
                dataframe.col <- dataframe.col + 1        
            }

        }      
              
        ## TODO
        iintronTable <- cbind(all.df, note2.df)
    }
                 
    if ((is.null(iintronTable)) || nrow(iintronTable) == 0)  
    {
        print(paste("iintronTable has 0 row for ", splicingEvents, sep=""))
        return(invisible(NULL))
    }
      
    # write intronTable
    f.intronTable <- iintronTable

    ## at this point, take off "note", "note" only work for intron retention
    note.ind <- grep(".note$", colnames(all.csv))
    if (length(note.ind)>0) all.csv <- all.csv[, -note.ind]

    ## statistical test? - propro.test
    if (splicingEvents %in% c("ri.3"))
    {
        nonjun.ind <- grep("\\.nonjun$", colnames(all.csv))
        jun.ind <- grep("\\.jun$", colnames(all.csv))
        
        nonjun <- all.csv[, nonjun.ind, drop=FALSE]
        sumjun <- all.csv[,nonjun.ind, drop=FALSE] + all.csv[,jun.ind, drop=FALSE]
        
        ## get rid of all zero because prop.test not accept
        row.sum <- rowSums(sumjun)
        row.max <-  apply(sumjun, 1, max)
        zero.ind <- which(row.sum==row.max)
        
        test.csv <- all.csv[-zero.ind,]
        notest.csv <- all.csv[zero.ind,]
        
        nonjun.2 <- test.csv[, nonjun.ind, drop=FALSE]
        sumjun.2 <- test.csv[,nonjun.ind, drop=FALSE] + test.csv[,jun.ind, drop=FALSE]
        
        pvalue <- rep(1, nrow(nonjun.2))
        sapply(1:nrow(nonjun.2), function(i)
              {
                  tmp <- prop.test(unlist(nonjun.2[i,]), unlist(sumjun.2[i,]))
                  pvalue[i] <<- tmp$p.value
              })
        
        test.csv$pvalue <- pvalue
        notest.csv$pvalue <- rep(1, nrow(notest.csv))
        all.csv <- rbind(test.csv, notest.csv)
    }
    
    if (splicingEvents %in% c("ri.3", "ri.4"))
    {
       nonjun.c <- all.csv[,grep("\\.nonjun$", colnames(all.csv))]
       ratio.c <- all.csv[,grep("\\.ratio$", colnames(all.csv))]
       if (ncol(nonjun.c) != ncol(ratio.c)) stop("The column for ratio is WRONG!")
       
       pvalue <- NULL
       sapply(1: nrow(nonjun.c), function(i)
            {                   
                
                if (any(as.numeric(unlist(ratio.c[i,]))==0)) pvalue <<- c(pvalue, "NA")
                else
                {
                  total <- as.numeric(unlist(nonjun.c[i,]))/as.numeric(unlist(ratio.c[i,]))
                  t <- prop.test(as.numeric(unlist(nonjun.c[i,])), total)
                  #if (inherits(f, "try-error")) browser()
                  pvalue <<- c(pvalue, t$p.value)
                }
                return(invisible(NULL))  
            })
       
       all.csv$pvalue <-  pvalue
       ## FDR adjust
            
             
    } 
    

    write.csv(all.csv, file=paste(intron.dir,"intronTable.csv", sep="/"))
    if (nrow(f.intronTable)>5000)
          sortable.html.table(f.intronTable[1:5000,], "intronTable.html",
                  output.directory=intron.dir)
    else
          sortable.html.table(f.intronTable, "intronTable.html",
                  output.directory=intron.dir)
    html <- paste(html,
          "<table><tbody><tr class=\"head\"> <td
          		 colspan=\"2\">",
          keyhtml, "</td></tr>
          <tr class=\"line1\"><td colspan=\"2\">
          reads description</td></tr>
          <tr class=\"line2\"><td width=\"50%\">Intron </td><td width=\"50%\">&nbsp;<a href=\"",
          file.path(intron.dir, "intronTable.html"), "\">Table</a>,&nbsp;<a href=\"",
          file.path(intron.dir, "intronTable.csv"), "\">csv</a></td>",
          "</tr>", "</td></tr></tbody></table><br><br>"
          , sep="")

    print("Generate intron TABLE section...END......")
    return(html)
}



