########Wormtracks server####
###Amhed Vargas
###amhed.velazquez@kaust.edu.sa
###Server

##Required packages
#install.packages("shiny")
#install.packages("shinythemes")
#install.packages("ggvis")
#install.packages("ggplot2")
#install.packages("DT")
#install.packages("shinyWidgets")

#Load libraries
library(shiny)
library(shinythemes)
library(ggvis)
library(ggplot2)
library(DT)
library(shinyWidgets)

##Add info for database search
#Library
PiLib=read.table("DB/Lib_info_simplified.tsv",sep="\t",header=F, stringsAsFactors=F)
colnames(PiLib)=c("Plate", "Well", "Pool", "Spot","PrimerOne","PrimerTwo","Gene","crRNA", "Type", "Oligo")

#Main DB with names
MainDB=read.table("DB/Main_DB.tsv",sep="\t",header=F, stringsAsFactors=F)
colnames(MainDB)=c("Accesion","ID","Locus","Transcript","Alias","Type")

  shinyServer(function(input, output, session) {
    
    ##Retrieve unique ID for the session
    session_id <- session$token
  
    ##Functions related to IGV browser
    # prints actual tab
    observeEvent(input$panels,{
      sendM(input$panels);
    })
    
    ##function to send check tab status and change visibility of browser
    sendM = function(x){
      session$sendCustomMessage("igvstat-change", x)
    }
    
    ###Control panels functions##########################################################################################
    ##Functions needed to generate links between panels
    observeEvent(input$link_to_tabpanel_genome_browser, {
      newvalue <- "Browse"
      updateTabsetPanel(session, "panels", newvalue)
    })
    
    observeEvent(input$link_to_tabpanel_title, {
      newvalue <- "Query"
      updateTabsetPanel(session, "panels", newvalue)
    })
    
    observeEvent(input$link_to_tabpanel_download, {
      newvalue <- "Downloads"
      updateTabsetPanel(session, "panels", newvalue)
    })
    observeEvent(input$link_to_tabpanel_about, {
      newvalue <- "About"
      updateTabsetPanel(session, "panels", newvalue)
    })
    
    geteinfo = function(x){
      if (is.null(x)) return(NULL)
      gidi=""
      if(x %in% as.character(ChrisPAT$Wormbase.ID)){gidi = x}
      if(x %in% as.character(ChrisPAT$Gene.name)){gidi = rownames(ChrisPAT[which(as.character(ChrisPAT$Gene.name) == x),])}
      if(gidi ==""){return(paste("Gene not found! try with Wormbase ID"))}
      
      gene <- ChrisPAT[as.character(gidi),]
      
      #runjs(paste("igv.browser.search(",gene$Chromosome,":",gene$Gene.start,"-",gene$Gene.end,")", sep=""))
      #paste0("<script>igv.browser.search(",gene$Chromosome,":",gene$Gene.start,"-",gene$Gene.end,")</script>")
      igvloc=paste(gene$Chromosome,":",gene$Gene.start,"-",gene$Gene.end,sep="")
      session$sendCustomMessage("gene-coordinates", igvloc)
    }
    
    #############################################################################################3
    ####Functions for search of gene
    
    ##Observers for action button search
    observeEvent(input$actiongenesearch, {
      output$ErrorMessage <- renderText({})
      wbid=""
      mygene = as.character(input$geneinput)
      
      if(mygene %in% as.character(MainDB$ID)){
        wbid=as.character(MainDB[which(as.character(MainDB$ID)==mygene)[1],1])
      }else{
      if(mygene %in% as.character(MainDB$Locus)){
        wbid=as.character(MainDB[which(as.character(MainDB$Locus)==mygene)[1],1])
      }else{
        
        if(mygene %in% as.character(MainDB$Transcript)){
          wbid=as.character(MainDB[which(as.character(MainDB$Transcript)==mygene)[1],1])
        }
        }
      }
      
      if(wbid == ""){
        output$SelPiTabSummary <- renderUI({
          HTML("<b>Gene not found</b>")
        })
        output$ErrorMessage <- renderText({})
        output$SelPiTab=DT::renderDataTable({})
        output$SimpleFragment <- renderText({})
        output$extraoui <- renderUI({})
      }else{
        tmpL=PiLib[which(as.character(PiLib$Gene) %in% wbid),]
        
        tata=tmpL[order(tmpL$Spot),]
        tmpL=tata[order(tata$Pool),]
        
        
        #output$SelPiTab<- DT::renderDataTable(tmpL)
        
        output$SelPiTab <- DT::renderDataTable({
          Pdata=data.frame(
            Plate = tmpL[,1],
            Well = tmpL[,2],
            Spot = tmpL[,4],
            PrimerFwd = tmpL[,5],
            PrimerRev= tmpL[,6],
            crRNA=tmpL[,8],
            Type=tmpL[,9],
            Oligo = shinyInput(actionButton, as.character(tmpL[,10]), 'button_', label = "Show", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' ),
            stringsAsFactors = FALSE
          )
          
          colnames(Pdata)[4]="Forward primer"
          colnames(Pdata)[5]="Reverse primer"
          
          Pdata
          #},server = FALSE, escape = FALSE, selection = 'none'))
        },server = FALSE, escape = FALSE, selection = 'none')
        
        ##Done
        
        }

    }, ignoreInit = T)
    
    ##Show functions
    ##Handle shiny to add dynamic button
    shinyInput <- function(FUN, seqs, id, ...) {
      inputs <- character(length(seqs))
      for (i in 1:length(seqs)) {
        inputs[i] <- as.character(FUN(paste0(id, seqs[i]), ...))
      }
      inputs
    }
    
    observeEvent(input$select_button, {
      selectedSeq <- as.character(strsplit(input$select_button, "_")[[1]][2])
      
      output$SimpleFragment <- renderText({paste0(selectedSeq)})
      
      })
    
    
    ##Observers for action button browser
    observeEvent(input$actionbrow, {
      gege =input$genebrow
      output$genebrowsearch <- renderUI({
        HTML(
          geteinfo(gege)
        )
      })
    }, ignoreInit = T)
    
    #######################FUnctions related to browser function
    #To do after first gui, basically copy and paste of info displayed after parsing
    
    ##Observe Js value
    observeEvent(input$jsValue, {
      ##Consider uncommenting lines only for parsing
      #cat("\nMessage Received\n")
      #cat(paste(input$jsValue),sep="\n")
      #output$igv_id <- renderText({
        #paste("Selected ID is:", input$jsValue)
       # paste(unlist(strsplit(input$jsValue,";")),sep="\n")
      #})
    
      message=as.character(input$jsValue)
      
      seq=unlist(strsplit(message,";"))[4]

      tmpL=PiLib[which(as.character(PiLib$crRNA) %in% seq),]
      
      tata=tmpL[order(tmpL$Spot),]
      tmpL=tata[order(tata$Pool),]
      
      
      #output$SelPiTab<- DT::renderDataTable(tmpL)
      
      output$SelPiTabBrowser <- DT::renderDataTable({
        Pdata=data.frame(
          Plate = tmpL[,1],
          Well = tmpL[,2],
          Spot = tmpL[,4],
          PrimerFwd = tmpL[,5],
          PrimerRev= tmpL[,6],
          crRNA=tmpL[,8],
          Type=tmpL[,9],
          Oligo = shinyInput(actionButton, as.character(tmpL[,10]), 'button_', label = "Show", onclick = 'Shiny.onInputChange(\"select_button2\",  this.id)' ),
          stringsAsFactors = FALSE
        )
        
        colnames(Pdata)[4]="Forward primer"
        colnames(Pdata)[5]="Reverse primer"
        
        Pdata
        #},server = FALSE, escape = FALSE, selection = 'none'))
      },server = FALSE, escape = FALSE, selection = 'none')
      
      ##Done
      
      })
    
    observeEvent(input$select_button2, {
      selectedSeq <- as.character(strsplit(input$select_button2, "_")[[1]][2])
      
      output$SimpleFragmentBrowser <- renderText({paste0(selectedSeq)})
      
    })
    
    })  

