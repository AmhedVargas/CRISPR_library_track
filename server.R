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
#install.packages("base64enc")

#Load libraries
library(shiny)
library(shinythemes)
library(ggvis)
library(ggplot2)
library(DT)
library(shinyWidgets)
library(base64enc)

##Add info for database search
#Library
PiLib=read.table("DB/Lib_info_simplified.tsv",sep="\t",header=F, stringsAsFactors=F)
colnames(PiLib)=c("Plate", "Well", "Pool", "Spot","PrimerOne","PrimerTwo","Gene","crRNA", "Type", "Oligo")

#Main DB with names
MainDB=read.table("DB/Main_DB.tsv",sep="\t",header=F, stringsAsFactors=F)
colnames(MainDB)=c("Accesion","ID","Locus","Transcript","Alias","Type")

  shinyServer(function(input, output, session) {
    
    ###############Session functions
    ##Retrieve unique ID for the session
    session_id <- session$token
    
    ##Create temporary folder for unique user
    system(paste("mkdir -p WorkingSpace/users/",session_id,sep=""))
    
    ##On leaving remove directory
    session$onSessionEnded(function(){
      system(paste("rm -rf WorkingSpace/users/",session_id,sep=""))
    }
    )
    
    ##Create path to remember
    UserPath = paste("WorkingSpace/users/",session_id,"/",sep="")
    
    ##Functions related to IGV browser
    # prints actual tab
    observeEvent(input$panels,{
      sendM(input$panels);
    })
    
    #######Messages
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
            Oligo = shinyInput(actionButton, paste(as.character(tmpL[,1]),"_",as.character(tmpL[,2]),"_",as.character(tmpL[,4]),"_",as.character(tmpL[,5]),"_",as.character(tmpL[,6]),"_",as.character(tmpL[,8]),"_",as.character(tmpL[,10]),sep=""), 'button_', label = "Download", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' ),
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
    
    
    ##Observers for action button browser
    observeEvent(input$actionbrow, {
      gege =input$genebrow
      output$genebrowsearch <- renderUI({
        HTML(
          geteinfo(gege)
        )
      })
    }, ignoreInit = T)
    
    
    
    ######################Functions to process oligos
    ##Make ape with oligo annotations
    OligoApe = function(sequence, FwdPrimerN, RevPrimerN, crRNASeq, Plate, Well, Spot){
      if (is.null(sequence)){return(NULL)}
      FileLines=c()
      ##Main definitions
      FileLines=append(FileLines,paste("LOCUS",paste(Well,Plate,Spot,crRNASeq,"OligoStructure",sep="_",collapse=""),paste(nchar(sequence),"bp ds-DNA", sep=""),"linear",paste(c(unlist(strsplit(date()," ")))[c(3,2,5)],sep="",collapse="-"),sep="     "))
      FileLines=append(FileLines,paste("DEFINITION",".",sep="     "))
      FileLines=append(FileLines,paste("ACCESSION",".",sep="     "))
      FileLines=append(FileLines,paste("VERSION",".",sep="     "))
      FileLines=append(FileLines,paste("SOURCE",".",sep="     "))
      FileLines=append(FileLines,paste("ORGANISM","C.elegans",sep="     "))
      
      ##Comments
      FileLines=append(FileLines,paste("COMMENT",paste("Plate:",as.character(Plate)),sep="     "))
      FileLines=append(FileLines,paste("COMMENT",paste("Well:",as.character(Well)),sep="     "))
      FileLines=append(FileLines,paste("COMMENT",paste("Spot:",as.character(Spot)),sep="     "))
      FileLines=append(FileLines,paste("COMMENT",paste("crRNA:",as.character(crRNASeq)),sep="     "))
      FileLines=append(FileLines,paste("COMMENT",paste(),sep="     "))
      FileLines=append(FileLines,paste("COMMENT",paste("Note: sequence of homology arms might differ from endogenous sequence as some were modified to prevent CRISPR re-cuting or enzyme digestion"),sep="     "))
      FileLines=append(FileLines,paste("COMMENT","Generated using wormbuilder.org",sep="     "))
      FileLines=append(FileLines,paste("COMMENT","ApEinfo:methylated:1",sep="     "))
      
      ##Features
      #Start
      FileLines=append(FileLines,paste("FEATURES             Location/Qualifiers",sep=""))
      #Constant info
      FileLines=append(FileLines,paste("     primer_bind     ","1..20",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","Universal Primer","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","Universal Primer","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#66ffcb","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#ff2600","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0\"",sep=""))
      FileLines=append(FileLines,paste("     promoter        ","98..167",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","PCeN50-2","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","PCeN50-2","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#0f7ffe","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#346ee0","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0\"",sep=""))
      FileLines=append(FileLines,paste("     misc_feature    ","43..46",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","GG Overhang","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","GG Overhang","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#ffd478","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#ffd478","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
      FileLines=append(FileLines,paste("     misc_feature    ","255..258",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","GG Overhang(1)","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","GG Overhang","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#ffd478","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#ffd478","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
      FileLines=append(FileLines,paste("     misc_feature    ","188..210",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","crRNA scaffold + Term","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","crRNA scaffold + Term","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#fc6fce","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","green","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0\"",sep=""))
      FileLines=append(FileLines,paste("     misc_feature    ","169..187",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","Spacer","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","Spacer","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#7980ff","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","green","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0\"",sep=""))
      FileLines=append(FileLines,paste("     misc_feature    ","211..254",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","Left Homology","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","Left Homology","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#75d5ff","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","green","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0\"",sep=""))
      FileLines=append(FileLines,paste("     misc_feature    ","47..90",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","Right Homology","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","Right Homology","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#75d5ff","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","green","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0\"",sep=""))
      FileLines=append(FileLines,paste("     primer_bind     ","complement(281..300)",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","Universal Primer(1)","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","Universal Primer","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#14c0bd","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#ff7d78","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0\"",sep=""))
      FileLines=append(FileLines,paste("     misc_feature    ","260..265",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","BsaI","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","BsaI","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#fefc78","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","green","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
      FileLines=append(FileLines,paste("     misc_feature    ","36..41",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","BsaI(1)","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","BsaI","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#fefc78","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","green","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
      FileLines=append(FileLines,paste("     misc_feature    ","92..97",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","XbaI","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","XbaI","\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#f1b1b4","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","green","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0\"",sep=""))
      FileLines=append(FileLines,paste("     primer_bind     ","complement(266..280)",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"",RevPrimerN,"\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"",RevPrimerN,"\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#7ffe07","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#fb0106","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0\"",sep=""))
      FileLines=append(FileLines,paste("     primer_bind     ","21..35",sep=""))
      FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"",FwdPrimerN,"\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"",FwdPrimerN,"\"",sep="",collapse=""),sep="     "))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#21fe80","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#fb0106","\"",sep=""))
      FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0\"",sep=""))
      
      ##Origin and sequence
      FileLines=append(FileLines,paste("ORIGIN"))
      
      Compseq=unlist(strsplit(sequence,""))
      partseq=c()
      
      for(i in seq(1,length(Compseq),10)){
        endseq=i+9
        if(length(Compseq)-i < 9){endseq=length(Compseq)}
        partseq=append(partseq,paste(Compseq[i:endseq],collapse=""))
        
      }
      
      i=1
      for(num in seq(1,length(Compseq),60)){
        index=as.character(num)
        spaces=paste(rep(" ",6-nchar(index)),collapse="")
        endseq=i+5
        if((length(partseq)-i) < 5){endseq=length(partseq)}
        FileLines=append(FileLines , paste(spaces,index," ",paste(partseq[i:(endseq)],collapse=" "),sep=""))
        
        i=i+6
      }
      
      FileLines=append(FileLines,paste("//"))
      return(FileLines)
      }
    
    
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
          Oligo = shinyInput(actionButton, paste(as.character(tmpL[,1]),"_",as.character(tmpL[,2]),"_",as.character(tmpL[,4]),"_",as.character(tmpL[,5]),"_",as.character(tmpL[,6]),"_",as.character(tmpL[,8]),"_",as.character(tmpL[,10]),sep=""), 'button_', label = "Download", onclick = 'Shiny.onInputChange(\"select_button2\",  this.id)' ),
          stringsAsFactors = FALSE
        )
        
        colnames(Pdata)[4]="Forward primer"
        colnames(Pdata)[5]="Reverse primer"
        
        Pdata
        #},server = FALSE, escape = FALSE, selection = 'none'))
      },server = FALSE, escape = FALSE, selection = 'none')
      
      ##Done
      
      })
    
    observeEvent(input$select_button, {
      selectedPlate <- as.character(strsplit(input$select_button, "_")[[1]][2])
      selectedWell <- as.character(strsplit(input$select_button, "_")[[1]][3])
      selectedSpot <- as.character(strsplit(input$select_button, "_")[[1]][4])
      selectedFwdP <- as.character(strsplit(input$select_button, "_")[[1]][5])
      selectedRevP <- as.character(strsplit(input$select_button, "_")[[1]][6])
      selectedSeq <- as.character(strsplit(input$select_button, "_")[[1]][7])
      selectedOligo <- as.character(strsplit(input$select_button, "_")[[1]][8])
      
      tmpLines = OligoApe(selectedOligo, selectedFwdP, selectedRevP, selectedSeq, selectedPlate, selectedWell, selectedSpot)
      
      writeLines(tmpLines,paste(UserPath,"oligo.gb",sep=""))
      
      outfile <- file.path(UserPath, "oligo.gb")
      
      b64 <- dataURI(
        file = outfile, 
        mime = "text/plain;charset=US-ASCII"
      )
      session$sendCustomMessage("downloadApe64", b64)
      
    })
    
    
    
    ##Observer of second output table
    observeEvent(input$select_button2, {
      
      selectedPlate <- as.character(strsplit(input$select_button2, "_")[[1]][2])
      selectedWell <- as.character(strsplit(input$select_button2, "_")[[1]][3])
      selectedSpot <- as.character(strsplit(input$select_button2, "_")[[1]][4])
      selectedFwdP <- as.character(strsplit(input$select_button2, "_")[[1]][5])
      selectedRevP <- as.character(strsplit(input$select_button2, "_")[[1]][6])
      selectedSeq <- as.character(strsplit(input$select_button2, "_")[[1]][7])
      selectedOligo <- as.character(strsplit(input$select_button2, "_")[[1]][8])
      
      #output$SimpleFragmentBrowser <- renderText({paste0(selectedSeq)})
      
      tmpLines = OligoApe(selectedOligo, selectedFwdP, selectedRevP, selectedSeq, selectedPlate, selectedWell, selectedSpot)
      
      writeLines(tmpLines,paste(UserPath,"oligo.gb",sep=""))
      
      outfile <- file.path(UserPath, "oligo.gb")

      b64 <- dataURI(
        file = outfile, 
        mime = "text/plain;charset=US-ASCII"
      )
      session$sendCustomMessage("downloadApe64", b64)
      
    })
    
    })  

