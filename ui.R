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
                     Shiny.addCustomMessageHandler("downloadApe64", function(b64) {
                     const a = document.createElement("a");
                     document.body.append(a);
                     a.download = "oligo.gb";
                     a.href = b64;
                     a.click();
                     a.remove();
                    });
                                 '),
    #Main tab pages
    navbarPage(
        ### Add action links to the tittle of piRNAi app
        title=actionLink("link_to_tabpanel_title", HTML("<b>crRNA library</b>")),
        windowTitle="crRNA library",
        id = "panels",
        
        tabPanel("Query",
                 mainPanel(
                     h2("Search for gene targets in our crRNA library"),
                         tabPanel("Gene search",
                                  br(),
                                  ##Basic search text
                                  textAreaInput("geneinput", label = "Target gene", value = "", resize="none", placeholder= "WormbaseID, transcript or gene name", rows=1),
                                  ##Basic search button
                                  actionButton("actiongenesearch", label = "Search"),
                                  hr(),
                                  ##Error in case something fails
                                  verbatimTextOutput("ErrorMessage"),
                                  verbatimTextOutput("SimpleFragment"),
                                  ##Basic table
                                  DT::dataTableOutput('SelPiTab'),
                                  htmlOutput("SelPiTabSummary"),
                                  uiOutput("extraoui")
                         ))),
        
        
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
                 hr(),
                 verbatimTextOutput("igv_id"),
                 verbatimTextOutput("ErrorMessageBrowser"),
                 verbatimTextOutput("SimpleFragmentBrowser"),
                 ##Basic table
                 DT::dataTableOutput('SelPiTabBrowser'),
                 )),
        ##Download of tracks
        tabPanel("Downloads",
                 mainPanel(
                     h3("Tracks"),
                     
                     HTML("<br>Bed track with crRNA library.
                              
                              <p align=\"justify\">
                         <a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/tracks/wormTracks/crRNAlib.v3.bed\">Download (48 MB)</a><br></p>")
                     )),

        ###About
        tabPanel("About",
                 mainPanel(
                 h3("The app"),
                 HTML("<p align=\"justify\">Wormbuilder tracks works as a repository for genomic coordinates of amenabe regions in our <a href=\"https://syngenbio.kaust.edu.sa/Pages/Home.aspx\">lab</a>...Still in development
                      </p>")
        )
    )),
    hr(),
    HTML("<a href=\"https://syngenbio.kaust.edu.sa/Pages/Home.aspx\">Syntetic genome biology laboratory @KAUST</a><br>"),
    HTML("<a href=\"http://www.wormbuilder.org/\">Wormbuilder</a><br>"),
    HTML("<a href=\"mailto:cfjensen@kaust.edu.sa\">Contact us!</a>")
    

)
)
