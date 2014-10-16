########################################
##  SplicingTypesAnno
##  - pipeline.r: report pipeline
########################################

# sampleList <- list(SampleName=c("a1", "a2"), 
#                    BamFiles=c("/work/xsun1/bam1.bam", "/work/xsun1/bam2.bam"),
#                    SampleID=c(1,2)
#                     )
splicingReport <- function(inputData, 
                        gtfFile, 
                        gtfChrFormat="chr", 
                        gtfGeneLabel="gene_id", 
                        gtfGeneValue="g",
                        gtfTranscriptLabel=NULL, ## "transcript_id"
                        outputDir="html", 
                        eventType=c("ri.1", "ri.2", "es.1", "es.2", 
                                  "adleft.1", "adleft.2", 
                                  "adright.1", "adright.2",
                                  "adboth.1", "adboth.2"),
                         inputFile,
                         minReadCounts=10,                     
                         selectGenes=NULL, 
                         novel="both", # "both", "yes", "no"
                         parallel=FALSE, cpus=2
                        )
{
    ## clean environment
    initEnv()
    
    ## parameter validate 
    if (missing(inputFile) & missing(inputData)) stop("Please provide either inputData or inputFile!")
    if (missing(gtfFile)) stop("Please provide gtf file!")

    if (minReadCounts < 2) stop("Only minReadCounts > =2 is accepted!")
        
    ## read sample file or sample data
    if (!missing(inputFile)) inputSampleFile(sampleFile=inputFile)
    else inputSampleFile(sampleList=inputData)

    ## directory check
    if (! file.exists(outputDir) )  dir.create(outputDir)

    oldpath <- getwd()
    setwd(outputDir)
    on.exit(setwd(oldpath))

    ## initialization
    html.riTable.1 <- NULL
    html.riTable.2 <- NULL
    html.esTable.1 <- NULL
    html.esTable.2 <- NULL
    html.adTable.l1 <- NULL
    html.adTable.l2 <- NULL
    html.adTable.r1 <- NULL
    html.adTable.r2 <- NULL
    html.adTable.both1 <- NULL
    html.adTable.both2 <- NULL

   ## bedfile initialization
    cleanBedFiles()
   
    ## start pipeline
    html.done <- paste(splicingENV$html.header,
          "<div class=\"title\"><img src=\"img/dot1.jpg\">",  
          splicingENV$PROJECT_NAME,   
          "<span  class=\"date_TIT\">",
          Sys.time(),
          "</span></div>", sep="")

   
   ## html sample section
   html.raw <- generateHTML.Raw(raw.dir = splicingENV$raw.dir
                          , outputRawHTML=splicingENV$outputRawHTML)

   ## annotate splicing types and store it in table.list()
    if (is.null(selectGenes))
    {
        # generate exon_intron structure from gtf/gff
        all.list <- translateGTF(gtfFile, gtfChrFormat=gtfChrFormat, 
                                  gtfGeneLabel=gtfGeneLabel, gtfGeneValue=gtfGeneValue,
                                  gtfTranscriptLabel=gtfTranscriptLabel, # memeory issue
                                  geneOverlap="no", ## TODO later
                                  exonBoundary="yes"
                                  )
        splicingENV$gene.GRange <- all.list$gene
        splicingENV$exonGRange <- all.list$exon
        splicingENV$intronGRange <- all.list$intron
        splicingENV$exonStart <- all.list$exonStart
        splicingENV$exonEnd <- all.list$exonEnd
        splicingENV$splicingLink <- all.list$splicingLink
        if (is.null(splicingENV$gene.GRange)) return(invisible(NULL))
        
        ## annotate splicing types
        generateHTML.annotateSample(splicingENV$intron.dir, 
                            as.character(splicingENV$spData$pData$BamFiles),
                            gtfFile, eventType,
                            minReadCounts=minReadCounts,
                            selectGenes=selectGenes, novel=novel,
                            parallel, cpus
                            )
    }
    else                    
    {   
        # generate exon_intron structure from gtf/gff      
        select.GRange <- translateGTF(gtfFile, gtfChrFormat=gtfChrFormat, 
                                  gtfGeneLabel=gtfGeneLabel, gtfGeneValue=gtfGeneValue,
                                  gtfTranscriptLabel=gtfTranscriptLabel, # memeory issue 
                                selectGenes=selectGenes,
                                exonBoundary="yes"
                                )
        if (is.null(select.GRange$gene)) return(invisible(NULL))

        # annotate splicing types
        generateHTML.annotateGene(splicingENV$intron.dir, 
                            as.character(splicingENV$spData$pData$BamFiles),
                            select.GRange,
                            eventType,
                            minReadCounts=minReadCounts,
                            selectGenes=selectGenes, novel=novel,
                            parallel, cpus
                            )        
        
    }


   ## html section for splicing types
   if ("ri.1" %in% eventType)
   {
     html.riTable.1 <- generateHTML.intronTable(splicingENV$RI1.dir, gtfFile, 
                              splicingENV$table.list, 
                              splicingEvents="ri.1",
                              "gene"
                              )
   }

   if ("ri.2" %in% eventType)
   {
     html.riTable.2 <- generateHTML.intronTable(splicingENV$RI2.dir, gtfFile, 
                              splicingENV$table.list, 
                              splicingEvents="ri.2",
                              "gene"
                              )
   }
   
   if ("es.1" %in% eventType)
   {
      html.esTable.1 <- generateHTML.intronTable(splicingENV$ES1.dir,gtfFile, 
                            splicingENV$table.list, 
                            splicingEvents="es.1",
                            "exon" 
                            )
   }

   if ("es.2" %in% eventType)
   {
      html.esTable.2 <- generateHTML.intronTable(splicingENV$ES2.dir,gtfFile, 
                            splicingENV$table.list, 
                            splicingEvents="es.2",
                            "exon" 
                            )

   }

   if ("adleft.1" %in% eventType)
   {
     html.adTable.l1 <- generateHTML.intronTable(splicingENV$ADLEFT1.dir,gtfFile,
                              splicingENV$table.list, 
                              splicingEvents="adleft.1",
                              "exon"
                              )
   }
   
   if ("adleft.2" %in% eventType)
   {
     html.adTable.l2 <- generateHTML.intronTable(splicingENV$ADLEFT2.dir,gtfFile,
                              splicingENV$table.list, 
                              splicingEvents="adleft.2",
                              "exon"
                              )
   }
   
   if ("adright.1" %in% eventType)
   {
      html.adTable.r1 <- generateHTML.intronTable(splicingENV$ADRIGHT1.dir,gtfFile,
                              splicingENV$table.list,
                              splicingEvents="adright.1",
                              "exon"
                              )

  }
  
   if ("adright.2" %in% eventType)
   {
      html.adTable.r2 <- generateHTML.intronTable(splicingENV$ADRIGHT2.dir,gtfFile,
                              splicingENV$table.list,
                              splicingEvents="adright.2",
                              "exon"
                              )

  }

   if ("adboth.1" %in% eventType)
   {
     html.adTable.both1 <- generateHTML.intronTable(splicingENV$ADBOTH1.dir,gtfFile,
                              splicingENV$table.list,
                              splicingEvents="adboth.1",
                              "exon"
                              )
   }
   
   if ("adboth.2" %in% eventType)
   {
     html.adTable.both2 <- generateHTML.intronTable(splicingENV$ADBOTH2.dir,gtfFile,
                              splicingENV$table.list,
                              splicingEvents="adboth.2",
                              "exon"
                              )
   }
  
  ## html summary section
  sampleName <- splicingENV$spData$pData$SampleName
  sampleID <- splicingENV$spData$pData$SampleID
  
    # data summary
   html.datasum <- generateHTML.splicingSummary(splicingENV$summary.dir, splicingENV$outputSummaryHTML, 
            sampleName, sampleID)

   
    # gene summary
   if (is.null(selectGenes)) html.geneSum <- generateHTML.Gene(splicingENV$gene.dir, splicingENV$outputGeneHTML, 
                                     sampleName, sampleID)

    # read summary
   html.readSum <- generateHTML.Read(splicingENV$read.dir, splicingENV$outputReadHTML, 
                                     sampleName, sampleID, selectGenes)

   ## pipeline output: 1. splicing types + summary
   html.done <- paste(html.done , html.raw, sep="")
   #if (is.null(selectGenes)) html.done <- paste(html.done, html.geneSum, sep="")
   html.done <- paste(html.done, html.readSum, sep="")
   
   html.done <- paste(html.done, html.datasum, sep="")
      
   if (!is.null(html.riTable.1)) html.done <- paste(html.done, html.riTable.1, sep="")
   if (!is.null(html.riTable.2)) html.done <- paste(html.done, html.riTable.2, sep="")
   if (!is.null(html.esTable.1)) html.done <- paste(html.done, html.esTable.1, sep="")
   if (!is.null(html.esTable.2)) html.done <- paste(html.done, html.esTable.2, sep="")
   if (!is.null(html.adTable.l1)) html.done <- paste(html.done, html.adTable.l1, sep="")
   if (!is.null(html.adTable.l2)) html.done <- paste(html.done, html.adTable.l2, sep="")
   if (!is.null(html.adTable.r1)) html.done <- paste(html.done, html.adTable.r1, sep="")
   if (!is.null(html.adTable.r2)) html.done <- paste(html.done, html.adTable.r2, sep="")
   if (!is.null(html.adTable.both1)) html.done <- paste(html.done, html.adTable.both1, sep="") 
   if (!is.null(html.adTable.both2)) html.done <- paste(html.done, html.adTable.both2, sep="")     
   
   ## pipeline output: 2. splicing visualization - bed file
   html.bedFile <- generateHTML.Bed(sampleName, eventType, outputDir)
   if (!is.null(html.bedFile)) html.done <- paste(html.done, html.bedFile, sep="")   
   
   ## pipeline end
   html.done <- paste(html.done, "<hr><p style=\"text-align:right\"><i> SplicingTypesAnno @ 2012</i></p><p></p><p></p></body></html>", sep="")

   seq.html <- file(splicingENV$html.index.name, "w")  
   writeLines(html.done, con = seq.html, sep = "\n")
   close(seq.html)
   
   ## final image file copied for html
    img.path <- system.file("img", package="SplicingTypesAnno")
    checkImg <- file.copy(img.path, ".", recursive = TRUE)
    if (!checkImg) print("You don't have permission to copy image files to the output folder!")
      
    
}


## browse the report
splicingShowReport <- function(outputDir="html")
{
    if (!file.exists(outputDir)) stop("The html folder does NOT exist!")
    
    file.location <- file.path(getwd(),outputDir, splicingENV$html.index.name)
    if (!file.exists(file.location)) stop("The report does NOT complete, and there is NO html file!")
    
    browseURL(file.location, browser = getOption("browser"))
    return(invisible(NULL))
}


splicingCleanReport <- function(outputDir="html")
{
    if (!file.exists(outputDir)) 
    { 
        print("The html folder does NOT exist! No need to clean the report")
        return(invisible(NULL))
    }
    else  unlink(outputDir, recursive=TRUE)
    return(invisible(NULL))
}

