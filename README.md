# crRNA_LibraryTrack
Shiny app that allows the exploration and adoption of a large scale CRISPR oligo library for *C.elegans*.

See the app in action [here](https://wormbuilder.dev/crRNALib/). Also, see the [crRNA-DB repository](https://github.com/AmhedVargas/crRNA_library_2021) for more information regarding the creation of crRNA sequences.

## Structure of the repository
The shiny app core lies onto two R scripts, the "server.R" and "ui.R" codes, and text-based databases of crRNA oligo sequences (contained within the "DB" folder). While the ui code handles the user queries, the server code computes them and load fragments from the database via the sequence identifier (see below).

## Dependencies
**R**

While the program has been successfully tested in R 3.6.1 (2019-07-05, and later versions) on a x86_64-pc-linux-gnu platform, the application should work fine in any other version as long as the following libraries are properly installed:

*   install.packages("shiny")
*   install.packages("shinythemes")
*   install.packages("DT")
*   install.packages("ggvis")
*   install.packages("ggplot2")
*   install.packages("shinyWidgets")
*   install.packages("base64enc")

## Deployment and implementation
You can see the app in action by going to [here](https://wormbuilder.dev/crRNALib/).

To run the app locally, make sure you have installed R along with all its dependencies. Follow by cloning this repository:

`git clone https://github.com/AmhedVargas/CRISPR_library_track`

Run the program via command line specifying an open port, e.g., 5100, and open a web browser to access the app.

`R -e "shiny::runApp('CRISPR_library_track',host="0.0.0.0",port=5100)"`

Alternatively, you can run the app in a graphical environment such as Rstudio.

## Usage
**Query tab**
Search for crRNA targets using WormBase IDs, locus names, or transcript ids.

**Browse tab**
We have incorporated the javascript version of IGV to show the location of the crRNAs. Feel free to click on the crRNAs to obtain their information.

**Basket tab**
Retrieve oligo sequences in bulk.

 
## Troubleshoot

Please feel free to [e-mail me](mailto:amhed.velazquez@kaust.edu.sa) for any question, doubt or error in the app.

## Citation
For now there is no available publication

