########Wormtracks server####
###Amhed Vargas
###amhed.velazquez@kaust.edu.sa
###Server

#Load libraries
library(shiny)
library(shinythemes)
library(ggvis)
library(ggplot2)
library(DT)
library(shinyWidgets)
library(base64enc)

##Create functions for Ape files
######################Functions to process oligos
##Make ape with oligo annotations
OligoApe = function(sequence, FwdPrimerN, RevPrimerN, crRNASeq, Plate, Well, Target){
	if (is.null(sequence)){return(NULL)}
	FileLines=c()
	##Main definitions
	FileLines=append(FileLines,paste("LOCUS",paste(Well,Plate,crRNASeq,"OligoStructure",sep="_",collapse=""),paste(nchar(sequence),"bp ds-DNA", sep=""),"linear",paste(c(unlist(strsplit(date()," ")))[c(3,2,5)],sep="",collapse="-"),sep="     "))
	FileLines=append(FileLines,paste("DEFINITION",".",sep="     "))
	FileLines=append(FileLines,paste("ACCESSION",".",sep="     "))
	FileLines=append(FileLines,paste("VERSION",".",sep="     "))
	FileLines=append(FileLines,paste("SOURCE",".",sep="     "))
	FileLines=append(FileLines,paste("ORGANISM","C.elegans",sep="     "))
	
	##Comments
	FileLines=append(FileLines,paste("COMMENT",paste("Plate:",as.character(Plate)),sep="     "))
	FileLines=append(FileLines,paste("COMMENT",paste("Well:",as.character(Well)),sep="     "))
	#FileLines=append(FileLines,paste("COMMENT",paste("Spot:",as.character(Spot)),sep="     "))
	FileLines=append(FileLines,paste("COMMENT",paste("Target:",as.character(Target)),sep="     "))
	FileLines=append(FileLines,paste("COMMENT",paste("Spacer:",as.character(crRNASeq)),sep="     "))
	FileLines=append(FileLines,paste("COMMENT",paste(),sep="     "))
	FileLines=append(FileLines,paste("COMMENT",paste("Note: sequence of homology arms might differ from endogenous sequence as some were modified to prevent CRISPR re-cuting or enzyme digestion"),sep="     "))
	FileLines=append(FileLines,paste("COMMENT",paste("Also, the first nucleotide of the guide sequence has been changed into G to promote its transcription."),sep="     "))
	FileLines=append(FileLines,paste("COMMENT",paste(),sep="     "))
	FileLines=append(FileLines,paste("COMMENT","Generated using wormbuilder.org",sep="     "))
	FileLines=append(FileLines,paste("COMMENT","ApEinfo:methylated:1",sep="     "))
	
	##Features
	#Start
	FileLines=append(FileLines,paste("FEATURES             Location/Qualifiers",sep=""))
	#Constant info
	FileLines=append(FileLines,paste("     primer_bind     ","1..20",sep=""))
	FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","UniF","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","UniF","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#66ffcb","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#ff2600","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0\"",sep=""))
	FileLines=append(FileLines,paste("     promoter        ","98..167",sep=""))
	FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","Pol III promoter","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","Pol III promoter","\"",sep="",collapse=""),sep="     "))
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
	FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","UniR","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","UniR","\"",sep="",collapse=""),sep="     "))
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

RBPOligoApe = function(sequence, GeneFwdPrimer, GeneRevPrimer, OligoFwdPrimer, OligoRevPrimer, crRNASeq, Plate, Well, Target){
	if (is.null(sequence)){return(NULL)}
	FileLines=c()
	##Main definitions
	FileLines=append(FileLines,paste("LOCUS",paste(Well,Plate,crRNASeq,"OligoStructure",sep="_",collapse=""),paste(nchar(sequence),"bp ds-DNA", sep=""),"linear",paste(c(unlist(strsplit(date()," ")))[c(3,2,5)],sep="",collapse="-"),sep="     "))
	FileLines=append(FileLines,paste("DEFINITION",".",sep="     "))
	FileLines=append(FileLines,paste("ACCESSION",".",sep="     "))
	FileLines=append(FileLines,paste("VERSION",".",sep="     "))
	FileLines=append(FileLines,paste("SOURCE",".",sep="     "))
	FileLines=append(FileLines,paste("ORGANISM","C.elegans",sep="     "))
	
	##Comments
	FileLines=append(FileLines,paste("COMMENT",paste("Plate:",as.character(Plate)),sep="     "))
	FileLines=append(FileLines,paste("COMMENT",paste("Well:",as.character(Well)),sep="     "))
	FileLines=append(FileLines,paste("COMMENT",paste("Target:",as.character(Target)),sep="     "))
	FileLines=append(FileLines,paste("COMMENT",paste("Spacer:",as.character(crRNASeq)),sep="     "))
	FileLines=append(FileLines,paste("COMMENT",paste(),sep="     "))
	FileLines=append(FileLines,paste("COMMENT","Generated using wormbuilder.org",sep="     "))
	FileLines=append(FileLines,paste("COMMENT","ApEinfo:methylated:1",sep="     "))
	
	##Features
	#Start
	FileLines=append(FileLines,paste("FEATURES             Location/Qualifiers",sep=""))
	#Constant info
	FileLines=append(FileLines,paste("     primer_bind     ","89..107",sep=""))
	FileLines=append(FileLines,paste("                     ",paste("/PCR_conditions=","\"","primer sequence:TAATACGACTCACTATAGG","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","T7 Promoter (including GG)","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/label=","\"","T7 Promoter (including GG)","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","T7 Promoter (including GG)","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#FFED7A","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#012ca2","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
	FileLines=append(FileLines,paste("     misc_feature    ","complement(266..280)",sep=""))
	FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","Gene-specific Reverse","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/label=","\"",GeneRevPrimer,"\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"",GeneRevPrimer,"\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#ff2600","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#ff2600","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
	FileLines=append(FileLines,paste("     misc_feature    ","complement(251..265)",sep=""))
	FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","Oligo-specific Reverse","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/label=","\"",OligoRevPrimer,"\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"",OligoRevPrimer,"\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#ff9393","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#ff9393","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
	FileLines=append(FileLines,paste("     misc_feature    ","36..50",sep=""))
	FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","Oligo-specific Forward","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/label=","\"",OligoFwdPrimer,"\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"",OligoFwdPrimer,"\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#00f900","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#00f900","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
	FileLines=append(FileLines,paste("     misc_feature    ","21..35",sep=""))
	FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","Gene-specific Forward","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/label=","\"",GeneFwdPrimer,"\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"",GeneFwdPrimer,"\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#d4fb78","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#d4fb78","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
	FileLines=append(FileLines,paste("     primer_bind     ","127..206",sep=""))
	FileLines=append(FileLines,paste("                     ",paste("/PCR_conditions=","\"","primer sequence:GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","sgRNA Scaffold","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/label=","\"","sgRNA Scaffold","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","sgRNA Scaffold","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#F3A7F7","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#c04829","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
	FileLines=append(FileLines,paste("     misc_feature    ","108..126",sep=""))
	FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","spacer","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/label=","\"","spacer","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","spacer","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#0096ff","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#0096ff","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
	FileLines=append(FileLines,paste("     primer_bind     ","206..211",sep=""))
	FileLines=append(FileLines,paste("                     ",paste("/PCR_conditions=","\"","primer sequence:TCTAGA","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","XbaI","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/label=","\"","XbaI","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","XbaI","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#ADB2FD","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#ff0000","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
	FileLines=append(FileLines,paste("     primer_bind     ","1..20",sep=""))
	FileLines=append(FileLines,paste("                     ",paste("/PCR_conditions=","\"","primer sequence:CTCACCGCTCTTGTAGCATG","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","UniversalForward","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/label=","\"","UniversalForward","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","UniversalForward","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#359245","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#db2626","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
	FileLines=append(FileLines,paste("     primer_bind     ","complement(281..300)",sep=""))
	FileLines=append(FileLines,paste("                     ",paste("/PCR_conditions=","\"","primer sequence:GACCGGCAATCTCTTCCTGG","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"","UniversalReverse","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/label=","\"","UniversalReverse","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"","UniversalReverse","\"",sep="",collapse=""),sep="     "))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"","#359245","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"","#db2626","\"",sep=""))
	FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
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

##Add Databases
#Names of primers
priminam=read.delim("DB/PrimerNaming.csv", sep=",",header=T, stringsAsFactors = F)

#Library
PiLib=read.table("DB/Libraries_data.tsv",sep="\t",header=T, stringsAsFactors=F)
#Rename primers based on sequences
rownames(priminam)=priminam$Sequence
PiLib$GenePFwd=priminam[PiLib$GenePFwd,"Name"]
PiLib$GenePRev=priminam[PiLib$GenePRev,"Name"]
PiLib$OligoPFwd=priminam[PiLib$OligoPFwd,"Name"]
PiLib$OligoPRev=priminam[PiLib$OligoPRev,"Name"]

#Main DB with names
MainDB=read.table("DB/Gene_info.tsv",sep="\t",header=F, stringsAsFactors=F)
colnames(MainDB)=c("Accesion","ID","Locus","Transcript","Alias","Type","Name2Use")

wbname=MainDB[,c("ID","Name2Use")]
rownames(wbname)=as.character(wbname[,1])

##crRNA coordinates
coordcrRNA = read.table("DB/Libraries.simplified.bed",sep="\t", header=F, stringsAsFactors=F)
colnames(coordcrRNA) = c("chr","midpos","seq")

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
    
    
    ###Basket values
    rv <- reactiveValues(
      # And here is our main data frame
      basket = data.frame(Plate = c(""), Well = c(""), 
      										OligoPFwd=c(""), OligoPRev=c(""), GenePFwd=c(""), GenePRev=c(""), 
      										Target=c(""),crRNA=c(""), Type=c(""), Oligo=c(""), RNPOligo=c(""), 
      										stringsAsFactors=F), 
      # And our counter
      counter = 0,
      #and igv location
      igvlox = ""
    )
    
    ##Make reactive stuff
    makeReactiveBinding("rv$basket")
    makeReactiveBinding("rv$counter")
    makeReactiveBinding("rv$igvlox")
    
    ###Testsss
    rdf <- reactive(rv$basket)
    rdc <- reactive(rv$counter)
    rdi <- reactive(rv$igvlox)
    
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
    
    #############################################################################################3
    ####Functions for search of gene
    ###Observer for multiple search 
    observeEvent(input$actionmultiplesearch, {
      
      output$ErrorMessageMultiple <- renderText({})
      
      wbidsIN=c()
      wbidsOut=c()
      oriris=c()
      
      ##Genes
      mygenes = as.character(input$MultipleGeness)
      mygenes = gsub(" ", "", mygenes)
      mygenes = gsub("\n", ",", mygenes)
      mygenes = unlist(strsplit(mygenes,","))
      
      ##Type of crRNA
      mytipin = as.character(input$crRNA_multiplegenes_type)
      
      if(length(mygenes) > 0){
        for(mygene in mygenes){
          
          if(mygene != ""){
          if(mygene %in% as.character(MainDB$ID)){
            wbidsIN=append(wbidsIN,as.character(MainDB[which(as.character(MainDB$ID)==mygene)[1],1]))
            oriris=append(oriris,mygene)
          }else{
            if(mygene %in% as.character(MainDB$Locus)){
              wbidsIN=append(wbidsIN,as.character(MainDB[which(as.character(MainDB$Locus)==mygene)[1],1]))
              oriris=append(oriris,mygene)
            }else{
              
              if(mygene %in% as.character(MainDB$Transcript)){
                wbidsIN=append(wbidsIN,as.character(MainDB[which(as.character(MainDB$Transcript)==mygene)[1],1]))
                oriris=append(oriris,mygene)
              }else{
                wbidsOut=append(wbidsOut,mygene)
                }
            }
          }
          
            }
          }
      }
      
      if(length(wbidsIN) == 0){
        output$MultipleExtraoui <- renderUI({
          HTML("<b>Not a single gene was found.</b>")
        })
        output$ErrorMessageMultiple <- renderText({})
        output$crRNAMultipleTab=DT::renderDataTable({})
        output$SimpleFragment <- renderText({})
      }else{
        if(length(wbidsOut)>0){
          output$ErrorMessageMultiple <- renderText({paste0(c("The following genes were not found:",wbidsOut),sep="\n")})
        }
        
        tmpT=PiLib[which(as.character(PiLib$Gene) %in% wbidsIN[1]),]
        #tmpT$Gene = rep(oriris[1],nrow(tmpT))
        
        tmpL = tmpT
        if(length(wbidsIN) > 1){
        for(q in 2:length(wbidsIN)){
          tmpT=PiLib[which(as.character(PiLib$Gene) %in% wbidsIN[q]),]
          #tmpT$Gene = rep(oriris[q],nrow(tmpT))
          tmpL=rbind(tmpL,tmpT)
          }
        
        }
        
        tmpL=unique(tmpL[,-4])
        
        crRNAmatches=c()
        
        ##Now, filter by type
        if(mytipin != "All"){
        mycrRNAtypes= as.character(tmpL$Gene)
        crRNAmatches=grep(mytipin,mycrRNAtypes)
        }else{
          crRNAmatches=1:nrow(tmpL)
          }
        
        ##
        if(length(crRNAmatches) == 0){
          output$ErrorMessageMultiple <- renderText({paste0("Some genes were found but there is no crRNA with the selected type. Try looking for all instead")})
          output$crRNAMultipleTab=DT::renderDataTable({})
          }else{
            tmpL=tmpL[crRNAmatches,]    
        tmpL=tmpL[order(tmpL$Well),]
        tmpL=tmpL[order(tmpL$Plate),]
        
        output$crRNAMultipleTab <- DT::renderDataTable({
          Pdata=data.frame(
          	Plate = tmpL$Plate,
          	Well = tmpL$Well,
          	OligoFwd = tmpL$OligoPFwd,
          	OligoRev= tmpL$OligoPRev,
          	GeneFwd = tmpL$GenePFwd,
          	GeneRev= tmpL$GenePRev,
          	Target= wbname[as.character(tmpL$Gene),2],
          	crRNA=tmpL$crRNA,
          	Type=tmpL$Type,
          	RNPOligo = shinyInput(actionButton, paste(as.character(tmpL[,1]),"_",as.character(tmpL[,2]),"_",as.character(tmpL[,3]),"_",as.character(tmpL[,4]),"_",as.character(tmpL[,5]),"_",as.character(tmpL[,6]),"_",as.character(tmpL[,7]),"_",as.character(wbname[as.character(tmpL[,8]),2]),"_",as.character(tmpL[,9]),"_",as.character(tmpL[,10]),"_",as.character(tmpL[,11]),"_",as.character(tmpL[,12]),"_RNP",sep=""), 'button_', label = "RNP GenBank", onclick = 'Shiny.onInputChange(\"select_buttonRNP\",  this.id.concat(\"_\", Math.random()))'),
						Oligo = shinyInput(actionButton, paste(as.character(tmpL[,1]),"_",as.character(tmpL[,2]),"_",as.character(tmpL[,3]),"_",as.character(tmpL[,4]),"_",as.character(tmpL[,5]),"_",as.character(tmpL[,6]),"_",as.character(tmpL[,7]),"_",as.character(wbname[as.character(tmpL[,8]),2]),"_",as.character(tmpL[,9]),"_",as.character(tmpL[,10]),"_",as.character(tmpL[,11]),"_",as.character(tmpL[,12]),"_Oligo",sep=""), 'button_', label = "DNA GenBank", onclick = 'Shiny.onInputChange(\"select_button\",  this.id.concat(\"_\", Math.random()))' ),
            Save = shinyInput(actionButton, paste(as.character(tmpL[,1]),"separator",as.character(tmpL[,2]),"separator",as.character(tmpL[,3]),"separator",as.character(tmpL[,4]),"separator",as.character(tmpL[,5]),"separator",as.character(tmpL[,6]),"separator",as.character(tmpL[,7]),"separator",as.character(tmpL[,8]),"separator",as.character(tmpL[,9]),"separator",as.character(tmpL[,10]),"separator",as.character(tmpL[,11]),"separator",as.character(tmpL[,12]),sep=""), 'buttonseparator', label = "Add to basket", onclick = 'Shiny.onInputChange(\"add_button2\",  this.id)' ),
            stringsAsFactors = FALSE
          )
          
          colnames(Pdata)[3]="Oligo-specific Fwd primer"
          colnames(Pdata)[4]="Oligo-specific Rev primer"
		  colnames(Pdata)[5]="Gene-specific Fwd primer (RNP only)"
          colnames(Pdata)[6]="Gene-specific Rev primer (RNP only)"
          colnames(Pdata)[8]="Guide sequence"
          colnames(Pdata)[10]="RNP annotated sequence (ApE)"
          colnames(Pdata)[11]="DNA annotated sequence (ApE)"
          Pdata
          #},server = FALSE, escape = FALSE, selection = 'none'))
        },server = FALSE, escape = FALSE, selection = 'none')
          }
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
    
    
    
    #######################FUnctions related to browser function
    #To do after first gui, basically copy and paste of info displayed after parsing
    observeEvent(input$locationIGV, {
      #cat(paste(input$locationIGV),sep="\n")
      #writeLines(as.character(input$locationIGV),paste(UserPath,"lastlocation.txt",sep=""))
      rv$igvlox <- as.character(input$locationIGV)
      })
    
    
    ##Observer for where I am in IGV browser
    observeEvent(input$userigvlocation, {
      #cat(rv$igvlox)
      
      tempis = unlist(strsplit(rv$igvlox,":"))
      uschr = tempis[1]
      
      uspos = unlist(strsplit(tempis[2],"-"))
      
      ussta = as.integer(gsub(pattern="," , replacement="", x=uspos[1]))
      usend = as.integer(gsub(pattern="," , replacement="", x=uspos[2]))
      
      tempos = which(coordcrRNA$chr %in% uschr)
      
      if(length(tempos) > 0){
        tempta = coordcrRNA[tempos,]
        tempos = which((tempta$midpos >= ussta) & (tempta$midpos <= usend))
      }
      
      if(length(tempos) > 0){
        tempseq = tempta[tempos,"seq"]
        
        idsss=which(as.character(PiLib$crRNA) %in% tempseq)
        
        #cat(paste(input$jsValue),sep="\n")
        
        #cat(paste(idsss),sep="\n")
        
        if(length(idsss) > 0){
          
          tmpL=PiLib[idsss,]
          
          tmpL=unique(tmpL[,-4])
          tmpL=tmpL[order(tmpL$Well),]
          tmpL=tmpL[order(tmpL$Plate),]
          
          
          #output$SelPiTab<- DT::renderDataTable(tmpL)
          
          output$SelPiTabBrowser <- DT::renderDataTable({
          	Pdata=data.frame(
          		Plate = tmpL$Plate,
          		Well = tmpL$Well,
          		OligoFwd = tmpL$OligoPFwd,
          		OligoRev= tmpL$OligoPRev,
          		GeneFwd = tmpL$GenePFwd,
          		GeneRev= tmpL$GenePRev,
          		Target= wbname[as.character(tmpL$Gene),2],
          		crRNA=tmpL$crRNA,
          		Type=tmpL$Type,
          		RNPOligo = shinyInput(actionButton, paste(as.character(tmpL[,1]),"_",as.character(tmpL[,2]),"_",as.character(tmpL[,3]),"_",as.character(tmpL[,4]),"_",as.character(tmpL[,5]),"_",as.character(tmpL[,6]),"_",as.character(tmpL[,7]),"_",as.character(wbname[as.character(tmpL[,8]),2]),"_",as.character(tmpL[,9]),"_",as.character(tmpL[,10]),"_",as.character(tmpL[,11]),"_",as.character(tmpL[,12]),"_RNP",sep=""), 'button_', label = "RNP GenBank", onclick = 'Shiny.onInputChange(\"select_buttonRNP\",  this.id.concat(\"_\", Math.random()))'),
          		Oligo = shinyInput(actionButton, paste(as.character(tmpL[,1]),"_",as.character(tmpL[,2]),"_",as.character(tmpL[,3]),"_",as.character(tmpL[,4]),"_",as.character(tmpL[,5]),"_",as.character(tmpL[,6]),"_",as.character(tmpL[,7]),"_",as.character(wbname[as.character(tmpL[,8]),2]),"_",as.character(tmpL[,9]),"_",as.character(tmpL[,10]),"_",as.character(tmpL[,11]),"_",as.character(tmpL[,12]),"_Oligo",sep=""), 'button_', label = "DNA GenBank", onclick = 'Shiny.onInputChange(\"select_button\",  this.id.concat(\"_\", Math.random()))' ),
          		Save = shinyInput(actionButton, paste(as.character(tmpL[,1]),"separator",as.character(tmpL[,2]),"separator",as.character(tmpL[,3]),"separator",as.character(tmpL[,4]),"separator",as.character(tmpL[,5]),"separator",as.character(tmpL[,6]),"separator",as.character(tmpL[,7]),"separator",as.character(tmpL[,8]),"separator",as.character(tmpL[,9]),"separator",as.character(tmpL[,10]),"separator",as.character(tmpL[,11]),"separator",as.character(tmpL[,12]),sep=""), 'buttonseparator', label = "Add to basket", onclick = 'Shiny.onInputChange(\"add_button2\",  this.id)' ),
          		stringsAsFactors = FALSE
          	)
          	
          	colnames(Pdata)[3]="Oligo-specific Fwd primer"
          	colnames(Pdata)[4]="Oligo-specific Rev primer"
          	colnames(Pdata)[5]="Gene-specific Fwd primer (RNP only)"
          	colnames(Pdata)[6]="Gene-specific Rev primer (RNP only)"
          	colnames(Pdata)[8]="Guide sequence"
          	colnames(Pdata)[10]="RNP annotated sequence (ApE)"
          	colnames(Pdata)[11]="DNA annotated sequence (ApE)"
            
            Pdata
          },server = FALSE, escape = FALSE, selection = 'none')
        }
        
      }else{
        ##Render no sequences here
        output$ErrorMessageBrowser <- renderText({paste("There is not a single guide sequence in this location.")})
        output$SelPiTabBrowser <- DT::renderDataTable({})
        }
      
    }, ignoreInit = T)
    
    
    ##Add to browser IGV location
    ##Observer for where I am in IGV browser
    observeEvent(input$userigvbasket, {
      ##Same as before
      tempis = unlist(strsplit(rv$igvlox,":"))
      uschr = tempis[1]
      uspos = unlist(strsplit(tempis[2],"-"))
      ussta = as.integer(gsub(pattern="," , replacement="", x=uspos[1]))
      usend = as.integer(gsub(pattern="," , replacement="", x=uspos[2]))
      tempos = which(coordcrRNA$chr %in% uschr)
      if(length(tempos) > 0){
        tempta = coordcrRNA[tempos,]
        tempos = which((tempta$midpos >= ussta) & (tempta$midpos <= usend))
      }
      if(length(tempos) > 0){
        tempseq = tempta[tempos,"seq"]
        idsss=which(as.character(PiLib$crRNA) %in% tempseq)
        if(length(idsss) > 0){
          tmpL=PiLib[idsss,]
          tmpL=unique(tmpL[,-4])
          tmpL=tmpL[order(tmpL$Well),]
          tmpL=tmpL[order(tmpL$Plate),]

          tmpsequs = cbind(tmpL[,1],tmpL[,2],tmpL[,4],tmpL[,5],tmpL[,6],tmpL[,7],tmpL[,8],tmpL[,9],tmpL[,10],tmpL[,11],tmpL[,12])
          
          colnames(tmpsequs) = colnames(rv$basket)
          
          rv$basket <- rbind(rv$basket,tmpsequs)
          rv$counter <- rv$counter + nrow(tmpsequs)
          
          showNotification(paste(nrow(tmpsequs), "sequences were added to the basket"))
        }

      }else{
        ##Render no sequences here
        output$ErrorMessageBrowser <- renderText({paste("There is not a single guide sequence to add to basket")})
        output$SelPiTabBrowser <- DT::renderDataTable({})
      }
      
    }, ignoreInit = T)
    
    
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
      
      idsss=which(as.character(PiLib$crRNA) %in% seq)
      
      #cat(paste(input$jsValue),sep="\n")
      
      #cat(paste(idsss),sep="\n")
      
      if(length(idsss) > 0){

      tmpL=PiLib[idsss,]
      
      tmpL=unique(tmpL[,-4])
      tmpL=tmpL[order(tmpL$Well),]
      tmpL=tmpL[order(tmpL$Plate),]
      
      
      #output$SelPiTab<- DT::renderDataTable(tmpL)
      
      output$SelPiTabBrowser <- DT::renderDataTable({
      	Pdata=data.frame(
      		Plate = tmpL$Plate,
      		Well = tmpL$Well,
      		OligoFwd = tmpL$OligoPFwd,
      		OligoRev= tmpL$OligoPRev,
      		GeneFwd = tmpL$GenePFwd,
      		GeneRev= tmpL$GenePRev,
      		Target= wbname[as.character(tmpL$Gene),2],
      		crRNA=tmpL$crRNA,
      		Type=tmpL$Type,
      		RNPOligo = shinyInput(actionButton, paste(as.character(tmpL[,1]),"_",as.character(tmpL[,2]),"_",as.character(tmpL[,3]),"_",as.character(tmpL[,4]),"_",as.character(tmpL[,5]),"_",as.character(tmpL[,6]),"_",as.character(tmpL[,7]),"_",as.character(wbname[as.character(tmpL[,8]),2]),"_",as.character(tmpL[,9]),"_",as.character(tmpL[,10]),"_",as.character(tmpL[,11]),"_",as.character(tmpL[,12]),"_RNP",sep=""), 'button_', label = "RNP GenBank", onclick = 'Shiny.onInputChange(\"select_buttonRNP\",  this.id.concat(\"_\", Math.random()))'),
      		Oligo = shinyInput(actionButton, paste(as.character(tmpL[,1]),"_",as.character(tmpL[,2]),"_",as.character(tmpL[,3]),"_",as.character(tmpL[,4]),"_",as.character(tmpL[,5]),"_",as.character(tmpL[,6]),"_",as.character(tmpL[,7]),"_",as.character(wbname[as.character(tmpL[,8]),2]),"_",as.character(tmpL[,9]),"_",as.character(tmpL[,10]),"_",as.character(tmpL[,11]),"_",as.character(tmpL[,12]),"_Oligo",sep=""), 'button_', label = "DNA GenBank", onclick = 'Shiny.onInputChange(\"select_button\",  this.id.concat(\"_\", Math.random()))' ),
      		Save = shinyInput(actionButton, paste(as.character(tmpL[,1]),"separator",as.character(tmpL[,2]),"separator",as.character(tmpL[,3]),"separator",as.character(tmpL[,4]),"separator",as.character(tmpL[,5]),"separator",as.character(tmpL[,6]),"separator",as.character(tmpL[,7]),"separator",as.character(tmpL[,8]),"separator",as.character(tmpL[,9]),"separator",as.character(tmpL[,10]),"separator",as.character(tmpL[,11]),"separator",as.character(tmpL[,12]),sep=""), 'buttonseparator', label = "Add to basket", onclick = 'Shiny.onInputChange(\"add_button2\",  this.id)' ),
      		stringsAsFactors = FALSE
      	)
      	
      	colnames(Pdata)[3]="Oligo-specific Fwd primer"
      	colnames(Pdata)[4]="Oligo-specific Rev primer"
      	colnames(Pdata)[5]="Gene-specific Fwd primer (RNP only)"
      	colnames(Pdata)[6]="Gene-specific Rev primer (RNP only)"
      	colnames(Pdata)[8]="Guide sequence"
      	colnames(Pdata)[10]="RNP annotated sequence (ApE)"
      	colnames(Pdata)[11]="DNA annotated sequence (ApE)"
        
        Pdata
      },server = FALSE, escape = FALSE, selection = 'none')
      }
      ##Done
      
      })
    
    
    ##Observer of second output table
    observeEvent(input$select_button, {
    	##plus one as there's an extra dash
    	selectedPlate <- as.character(strsplit(input$select_button, "_")[[1]][2])
    	selectedWell <- as.character(strsplit(input$select_button, "_")[[1]][3])
    	#selectedSpot <- as.character(strsplit(input$select_buttonRNP, "_")[[1]][4])
    	selectedOFwdP <- as.character(strsplit(input$select_button, "_")[[1]][5])
    	selectedORevP <- as.character(strsplit(input$select_button, "_")[[1]][6])
    	selectedGFwdP <- as.character(strsplit(input$select_button, "_")[[1]][7])
    	selectedGRevP <- as.character(strsplit(input$select_button, "_")[[1]][8])
    	selectedTarget <- as.character(strsplit(input$select_button, "_")[[1]][9])
    	selectedSeq <- as.character(strsplit(input$select_button, "_")[[1]][10])
    	selectedType <- as.character(strsplit(input$select_button, "_")[[1]][11])
    	selectedOligo <- as.character(strsplit(input$select_button, "_")[[1]][12])
    	selectedRNPOligo <- as.character(strsplit(input$select_button, "_")[[1]][13])
      
      #output$SimpleFragmentBrowser <- renderText({paste0(selectedSeq)})
      
      tmpLines = OligoApe(selectedOligo, selectedOFwdP, selectedORevP, selectedSeq, selectedPlate, selectedWell, selectedTarget)
      
      writeLines(tmpLines,paste(UserPath,"oligo.gb",sep=""))
      
      outfile <- file.path(UserPath, "oligo.gb")

      b64 <- dataURI(
        file = outfile, 
        mime = "text/plain;charset=US-ASCII"
      )
      
      fnamas=paste(paste("DNA",selectedPlate, selectedWell,selectedOFwdP, selectedORevP, selectedTarget, selectedType, selectedSeq, collapse="", sep="_"),".gb",sep="")
      
      params=c(fnamas,b64)
      session$sendCustomMessage("downloadApe64", params)
      
    })
    
	##Observer of second output table
    observeEvent(input$select_buttonRNP, {
    	##plus one as there's an extra dash
    	selectedPlate <- as.character(strsplit(input$select_buttonRNP, "_")[[1]][2])
    	selectedWell <- as.character(strsplit(input$select_buttonRNP, "_")[[1]][3])
    	#selectedSpot <- as.character(strsplit(input$select_buttonRNP, "_")[[1]][4])
    	selectedOFwdP <- as.character(strsplit(input$select_buttonRNP, "_")[[1]][5])
    	selectedORevP <- as.character(strsplit(input$select_buttonRNP, "_")[[1]][6])
    	selectedGFwdP <- as.character(strsplit(input$select_buttonRNP, "_")[[1]][7])
    	selectedGRevP <- as.character(strsplit(input$select_buttonRNP, "_")[[1]][8])
    	selectedTarget <- as.character(strsplit(input$select_buttonRNP, "_")[[1]][9])
    	selectedSeq <- as.character(strsplit(input$select_buttonRNP, "_")[[1]][10])
    	selectedType <- as.character(strsplit(input$select_buttonRNP, "_")[[1]][11])
    	selectedOligo <- as.character(strsplit(input$select_buttonRNP, "_")[[1]][12])
    	selectedRNPOligo <- as.character(strsplit(input$select_buttonRNP, "_")[[1]][13])
    	
    	
    	tmpLines = RBPOligoApe(selectedRNPOligo, selectedGFwdP, selectedGRevP, selectedOFwdP, selectedORevP, selectedSeq, selectedPlate, selectedWell, selectedTarget)
    	
    	writeLines(tmpLines,paste(UserPath,"oligo.gb",sep=""))
    	
    	outfile <- file.path(UserPath, "oligo.gb")
    	
    	b64 <- dataURI(
    		file = outfile, 
    		mime = "text/plain;charset=US-ASCII"
    	)
    	
    	fnamas=paste(paste("RNP",selectedPlate, selectedWell, selectedGFwdP, selectedGRevP, selectedOFwdP, selectedORevP, selectedTarget, selectedType, selectedSeq, collapse="", sep="_"),".gb",sep="")
    	
    	params=c(fnamas,b64)
    	session$sendCustomMessage("downloadApe64", params)
    	
    })
    
    
    ##Add to basket on genome browser
    observeEvent(input$add_button2, {
    	##plus one as there's an extra dash
    	selectedPlate <- as.character(strsplit(input$add_button2, "separator")[[1]][2])
    	selectedWell <- as.character(strsplit(input$add_button2, "separator")[[1]][3])
    	#selectedSpot <- as.character(strsplit(input$select_buttonRNP, "_")[[1]][4])
    	selectedOFwdP <- as.character(strsplit(input$add_button2, "separator")[[1]][5])
    	selectedORevP <- as.character(strsplit(input$add_button2, "separator")[[1]][6])
    	selectedGFwdP <- as.character(strsplit(input$add_button2, "separator")[[1]][7])
    	selectedGRevP <- as.character(strsplit(input$add_button2, "separator")[[1]][8])
    	selectedTarget <- as.character(strsplit(input$add_button2, "separator")[[1]][9])
    	selectedSeq <- as.character(strsplit(input$add_button2, "separator")[[1]][10])
    	selectedType <- as.character(strsplit(input$add_button2, "separator")[[1]][11])
    	selectedOligo <- as.character(strsplit(input$add_button2, "separator")[[1]][12])
    	selectedRNPOligo <- as.character(strsplit(input$add_button2, "separator")[[1]][13])
    	
      rv$basket <- rbind(rv$basket,c(selectedPlate,selectedWell,
      															 selectedOFwdP,selectedORevP, selectedGFwdP,selectedGRevP,
      															 selectedTarget, selectedSeq,selectedType,
      															 selectedOligo,selectedRNPOligo))
      rv$counter <- rv$counter + 1
      
      showNotification(paste(selectedSeq, "was added to the basket"))
    })
    
    ##Add to basket on genome browser
    observeEvent(input$remove_button, {
      
      removeID <- as.integer(strsplit(input$remove_button, "_")[[1]][3])
      
      idss = c(1:nrow(rv$basket))[-removeID]
      rv$basket <- rv$basket[idss,]
      rv$counter <- rv$counter - 1
      
    })
    
    observeEvent(rv$counter,{
      
      if(rv$counter == 0){
        output$BigBasket <- DT::renderDataTable({})
        output$DownloadBasket <- renderUI({})
        }else{
          output$DownloadBasket <- renderUI({
            tagList(
          downloadButton('DownBasketTable', 'Download table'),
          downloadButton('DownBasketRNP', 'Download RNP annotated genbank files in bulk'),
          downloadButton('DownBasketApe', 'Download DNA annotated genbank files in bulk'),
          actionButton("Cleanbasket", label = "Empty basket")
            )
          })
          
      ##Add buttons to remove and download
      output$BigBasket <- DT::renderDataTable({
        dtt=rv$basket
        dtt=dtt[-1,]
        
        Pdata=data.frame(
          Plate = dtt[,1],
          Well = dtt[,2],
          OligoPFwd = dtt[,3],
          OligoPRev= dtt[,4],
          GenePFwd = dtt[,5],
          GenePRev= dtt[,6],
          Target= wbname[as.character(dtt[,7]),2],
          crRNA=dtt[,8],
          Type=dtt[,9],
          RNPOligo = shinyInput(actionButton, paste(as.character(dtt[,1]),"_",as.character(dtt[,2]),"_",as.character(dtt[,2]),"_",as.character(dtt[,3]),"_",as.character(dtt[,4]),"_",as.character(dtt[,5]),"_",as.character(dtt[,6]),"_",as.character(wbname[as.character(dtt[,7]),2]),"_",as.character(dtt[,8]),"_",as.character(dtt[,9]),"_",as.character(dtt[,10]),"_",as.character(dtt[,11]),"_RNPBk",sep=""), 'button_', label = "RNP GenBank", onclick = 'Shiny.onInputChange(\"select_buttonRNP\",  this.id.concat(\"_\", Math.random()))'),
          Oligo = shinyInput(actionButton, paste(as.character(dtt[,1]),"_",as.character(dtt[,2]),"_",as.character(dtt[,2]),"_",as.character(dtt[,3]),"_",as.character(dtt[,4]),"_",as.character(dtt[,5]),"_",as.character(dtt[,6]),"_",as.character(wbname[as.character(dtt[,7]),2]),"_",as.character(dtt[,8]),"_",as.character(dtt[,9]),"_",as.character(dtt[,10]),"_",as.character(dtt[,11]),"_OligoBk",sep=""), 'button_', label = "DNA GenBank", onclick = 'Shiny.onInputChange(\"select_button\",  this.id.concat(\"_\", Math.random()))' ),
          Remove = shinyInput(actionButton, paste("Basket_",as.character((1:nrow(dtt))+1),sep=""), 'button_', label = "Remove from basket", onclick = 'Shiny.onInputChange(\"remove_button\",  this.id.concat(\"_\", Math.random()))' ),
          stringsAsFactors = FALSE
        )
           
        colnames(Pdata)[3]="Oligo-specific Fwd primer"
        colnames(Pdata)[4]="Oligo-specific Rev primer"
        colnames(Pdata)[5]="Gene-specific Fwd primer (RNP only)"
        colnames(Pdata)[6]="Gene-specific Rev primer (RNP only)"
        colnames(Pdata)[8]="Guide sequence"
        colnames(Pdata)[10]="RNP annotated sequence (ApE)"
        colnames(Pdata)[11]="DNA annotated sequence (ApE)"
        
        
        rownames(Pdata)=1:nrow(Pdata)
        Pdata
        #},server = FALSE, escape = FALSE, selection = 'none'))
      },server = FALSE, escape = FALSE, selection = 'none')
      
        }
      
      })
    
    ##Download table
    output$DownBasketTable <- downloadHandler(
      filename = function() {
        paste("Basket", "tsv", sep=".")
      },
      content = function(file) {
        dtt=rv$basket
        dtt=dtt[-1,]
        ##Change names
        dtt[,7]=as.character(wbname[as.character(dtt[,7]),2])
        colnames(dtt)[8] = "Guide sequence"
        
        fname = paste(UserPath,"Basket.tsv",sep="")
        write.table(x=dtt,fname,row.names=F,sep="\t", quote = F)
        
        file.copy(fname, file)
      }
    )
    
    ###Download zip files
    output$DownBasketRNP <- downloadHandler(
    	filename = function() {
    		paste("RNP-Basket", "zip", sep=".")
    	},
    	content = function(file) {
    		dtt=rv$basket
    		dtt=dtt[-1,]
    		
    		fs <- c()
    		fname = paste(UserPath,"Basket.zip",sep="")
    		for(i in 1:nrow(dtt)){
    			tmpLines = RBPOligoApe(dtt[i,11], dtt[i,5], dtt[i,6], dtt[i,3], dtt[i,4], dtt[i,8], dtt[i,1],dtt[i,2],wbname[as.character(dtt[i,7]),2])
    			##Add type so files are not collapsed
    			ppath=paste(UserPath,paste("RNP",dtt[i,1], dtt[i,2],dtt[i,5], dtt[i,6], dtt[i,3], dtt[i,4], wbname[as.character(dtt[i,7]),2], dtt[i,9], dtt[i,8], sep="_"),".gb",sep="")
    			writeLines(tmpLines,ppath)
    			fs <- c(fs, ppath)
    		}
    		
    		zip(zipfile=fname, files=fs, flags= "-r9Xj")
    		
    		file.copy(fname, file)
    	},
    	contentType = "application/zip"
    )
    
    output$DownBasketApe <- downloadHandler(
      filename = function() {
        paste("DNA-Basket", "zip", sep=".")
      },
      content = function(file) {
        dtt=rv$basket
        dtt=dtt[-1,]

        fs <- c()
        fname = paste(UserPath,"Basket.zip",sep="")
        for(i in 1:nrow(dtt)){
        	tmpLines = OligoApe(dtt[i,10], dtt[i,3], dtt[i,4], dtt[i,8], dtt[i,1],dtt[i,2],wbname[as.character(dtt[i,7]),2])          ##Add type so files are not collapsed
          ppath=paste(UserPath,paste("DNA",dtt[i,1], dtt[i,2],dtt[i,3], dtt[i,4], wbname[as.character(dtt[i,7]),2], dtt[i,9], dtt[i,8], sep="_"),".gb",sep="")
          writeLines(tmpLines,ppath)
          fs <- c(fs, ppath)
          }
        
        zip(zipfile=fname, files=fs, flags= "-r9Xj")
        
        file.copy(fname, file)
      },
      contentType = "application/zip"
    )
    
    ##Clean basket
    ##Add to basket on genome browser
    observeEvent(input$Cleanbasket, {
      rv$basket <- data.frame(Plate = c(""), Well = c(""), 
      												OligoPFwd=c(""), OligoPRev=c(""), GenePFwd=c(""), GenePRev=c(""), 
      												Target=c(""),crRNA=c(""), Type=c(""), Oligo=c(""), RNPOligo=c(""), 
      												stringsAsFactors=F)
      rv$counter <- 0
    })
    
    
    })  

