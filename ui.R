########crRNA library interface####
###Amhed Vargas
###amhed.velazquez@kaust.edu.sa
###UI

#Load libraries
library(shiny)
library(shinythemes)
library(ggvis)
library(ggplot2)
library(DT)
library(shinyWidgets)
library(base64enc)

# Define User interface
shinyUI(
    fluidPage(
        ##Add same bootstrap theme as css
        tags$head(
            tags$link(rel="stylesheet",type = "text/css", href="bootstrap.min.css")
        ),
        
    ##Custom extra styles: single sliders background and title of navbar as transparent 
    tags$style(type = 'text/css', 
               ".js-irs-none .irs-single, .js-irs-none .irs-bar-edge, .js-irs-none .irs-bar {
                          background: transparent;
                          border-top-color: transparent;
                          border-bottom-color: transparent;
                          border-left-color: transparent;
                          border-right-color: transparent}
               .navbar-default .navbar-brand:hover {color: #ffffff;}
               "),
    #script to renew igv browser everytime a tab is changed. This solves issue regarding Nan positioning when browser is not in main tab.
    tags$script("
                     Shiny.addCustomMessageHandler('igvstat-change', function(panel) {
                     igv.visibilityChange()
                     });
                                 "),
    ##Script to create Ape annotations based on 64 bytes encoding
    tags$script('
                     Shiny.addCustomMessageHandler("downloadApe64", function(params) {
                     var namas = params[0];
                     var b64 = params[1];
                     const a = document.createElement("a");
                     document.body.append(a);
                     //a.download = "oligo.gb";
                     //a.download = b64.mime;
                     a.download = namas;
                     a.href = b64;
                     a.click();
                     a.remove();
                    });
                                 '),
    #Main tab pages
    navbarPage(
        ### Add action links to the tittle of piRNAi app
        title=actionLink("link_to_tabpanel_title", HTML("<b><i>C. elegans</i> CRISPR-Cas9 library</b>")),
        windowTitle="C. elegans CRISPR-Cas9 library",
        id = "panels",
        
        tabPanel("Query",
                 fluidRow(
                     h2("Search for gene targets"),
                     
                         #tabPanel("Multiple search",
                         tabPanel("Gene search",
                                 br(),
                                 textAreaInput("MultipleGeness", placeholder = "Wormbase IDs, transcript or gene names separated by newlines or commas", label = "List of genes", value = "", cols= 100, rows=5, width = "600px"),
                                 selectInput("crRNA_multiplegenes_type", label = HTML("<b>Type (location)</b>"), 
                                             choices = c("All","500 bp promoter", "250 bp promoter", "ATG", "CDS", "Stop"), 
                                             selected = 1),
                                 ##Basic search button
                                 actionButton("actionmultiplesearch", label = "Search"),
                                 hr(),
                                 ##Error in case something fails
                                 verbatimTextOutput("ErrorMessageMultiple"),
                                 ##Basic table
                                 DT::dataTableOutput('crRNAMultipleTab'),
                                 htmlOutput("MultipleExtraoui"),
                                 class = "span7")
                     #)
                     )),
        
        
        ####Genome browser
        tabPanel("Browse",
                 #mainPanel(
                 fluidRow(
                   h3("Browse your favorite location"),
                   includeScript("https://igv.org/web/release/2.10.5/dist/igv.js"),
                   includeHTML("www/igv.html"),
                   hr(),
                   # Listen for gene-coordinates messages
                     tags$script("
                     Shiny.addCustomMessageHandler('gene-coordinates', function(coordinates) {
                     igv.browser.search(coordinates);
                     });
                                 "),
                   htmlOutput("genebrowsearch"),
                  #   ),
                 #hr(),
                 verbatimTextOutput("igv_id"),
                 verbatimTextOutput("ErrorMessageBrowser"),
                 verbatimTextOutput("SimpleFragmentBrowser"),
                 actionButton("userigvlocation", label = "Display guide sequences seen in browser"),
                 actionButton("userigvbasket", label = "Send guide sequences seen in browser to basket"),
                 hr(),
                 ##Basic table
                 DT::dataTableOutput('SelPiTabBrowser')
                 )),
        
        ##Basket tab
        tabPanel("Basket",
                 fluidRow(
                     h3("Download oligos in bulk"),
                     DT::dataTableOutput('BigBasket'),
                     ##Add a conditional panel with rv counter higher than 0 to add download button, alternatively associate that to the value # this option is better so
                     uiOutput("DownloadBasket")
                 )
        ),
        
        ##Download of tracks
        tabPanel("Downloads",
                 fluidRow(
                     h3("Tracks"),
                     
                     HTML("<br>Bed track with protospacers (ce11/WS282)
                              <p align=\"justify\">
                         <a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/tracks/wormTracks/crRNAlib.v4.bed\">Download (45.5 MB)</a><br></p>"),
                     
                     HTML("<br>BigWig track with AlphaFold confidence scores (ce11/WS282)
                              <p align=\"justify\">
                         <a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/tracks/wormTracks/AF_confidence.bigwig\">Download (145 MB)</a><br></p>")
                     )),

        ###About
        tabPanel("About",
                 fluidRow(
                     HTML("<h3><i>C. elegans</i> CRISPR-Cas9 library DB</h3>"),
                     HTML("<p align=\"justify\">
                      This website is generated via custom modified css/html code running in R via the shiny library.
                 <br>All the templates, libraries, and programs used to produce this site are under the MIT and GNU licenses.
                    <br>The R scripts and files to render this site can be found <a href=\"https://github.com/AmhedVargas/CRISPR_library_track\">here</a></p>
                 This website was designed by <a href=\"https://scholar.google.com/citations?user=WokAqF8AAAAJ&hl=en&authuser=1\">Mohammed AlJohani</a> (<a href=\"https://www.linkedin.com/in/mohammedaljohani/\">LinkedIn</a>), <a href=\"https://www.researchgate.net/profile/Amhed_Vargas_Velazquez\">Amhed Missael Vargas Velazquez</a> (<a href=\"https://github.com/AmhedVargas\">Github</a>), and 
                          <a href=\"https://www.kaust.edu.sa/en/study/faculty/christian-jensen\">Christian Froekjaer-Jensen</a> from the Laboratory of Synthetic Genome Biology<br>"),
                 h3("The Laboratory of Synthetic Genome Biology"),
                 HTML("<p align=\"justify\">
                 The Laboratory of Synthetic Genome Biology is located in building 2 - level 3 (Ibn Al-Haytham – Above Spine) at King Abdullah University of Science and Technology (KAUST).
                 <br><i>Contact info</i>:<br>Christian Froekjaer-Jensen, Ph.D. 
                 <br>Assistant Professor of Bioscience
                 <br><a href=\"https://wormbuilder.org/laboratory-of-synthetic-genome-biology/\">Laboratory of Synthetic Genome Biology</a>
                 <br>Email: <a href=\"mailto:cfjensen@kaust.edu.sa\">cfjensen@kaust.edu.sa</a>
                      </p>")
        )
    )),
    hr(),
    HTML("<a href=\"https://www.kaust.edu.sa/en/study/faculty/christian-jensen\">Froekjaer-Jensen Lab @KAUST</a><br>"),
    HTML("<a href=\"http://www.wormbuilder.org/\">Wormbuilder</a><br>"),
    HTML("<a href=\"mailto:cfjensen@kaust.edu.sa\">Contact us!</a>")
    

)
)
