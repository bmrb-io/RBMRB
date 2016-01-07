
devtools::use_package("httr")
devtools::use_package("reshape2")
#'Downloads the chemical shift data from BMRB for a given BMRB entry/list of BMRB entries
#'
#'@param BMRBidlist ==> sinlge BMRB ID / list of BMRB IDs in csv format
#'@return all available chemical shift data in R data frame
#'@export fetchBMRB
#'@examples
#'df<-fetchBMRB('15060')
#'df<-fetchBMRB('15060,15070,8898,99')

fetchBMRB<-function(BMRBidlist){
  rawdata<-httr::GET('http://manta.bmrb.wisc.edu/api/a/chemshifts.php',
                     query=list(idlist=BMRBidlist))
  c<-httr::content(rawdata,'text')
  if (nchar(c)>5){
    d<-gsub("\\]","",gsub("\\[","",c))
    d2<-gsub("\n\"Entry_ID\",\"id\",\"Entity_ID\",\"Comp_index_ID\",\"Comp_ID\",\"Atom_ID\",\"Atom_type\",\"Val\",\"Val_err\",\"Ambiguity_code\",\"Assigned_chem_shift_list_ID\"","",d)
    t<-read.table(textConnection(d2),sep="\n")
    outdata<-reshape2::colsplit(t$V1,",",names=c('x','BMRB_ID','Entry_ID','Entity_ID','Comp_index_ID','Comp_ID','Atom_ID','Atom_type','Chemical_shift','err','Ambiguity_code','Assigned_chem_shift_list_ID'))
    outdata$x<-NULL
  }
  else{
    cat("Invalid BMRB ID")
    outdata<-NA
  }
  return (outdata)
}


#'Converts chemical shift data frame into H1-N15 HSQC data frame
#'
#'@param csdf ==> chemical shift data frame from fetchBMRB
#'@return 1H-N15 chemical shift list on the same row combined using comp index ID and bmrb ID
#'@export N15HSQC
#'@examples
#'df<-fetchBMRB('15060')
#'hsqc<-N15HSQC(df)

N15HSQC<-function(csdf){
  shiftH<-subset(csdf,Atom_ID=="H")
  names(shiftH)[names(shiftH)=="Chemical_shift"]<-"H"
  shiftN<-subset(csdf,Atom_ID=="N")
  names(shiftN)[names(shiftN)=="Chemical_shift"]<-"N"
  shiftHN<-merge(shiftH,shiftN,by=c('BMRB_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  outdat<-shiftHN[,c("BMRB_ID","Comp_index_ID","Assigned_chem_shift_list_ID","Comp_ID.x","Comp_ID.y","H","N")]
  names(outdat)[names(outdat)=="Comp_ID.x"]<-"Comp_ID_H"
  names(outdat)[names(outdat)=="Comp_ID.y"]<-"Comp_ID_N"
  return(outdat)
}

#'Reads the full chemical shift csv frile from BMRB ftp site. Not recommended unless you have high speed internet connection
#'
#'@param No parameter is needed
#'@return list of all atom chemical shifts for all BMRB entries as a R data frame
#'@export fetchallBMRB
#'@examples
#'df<-fetchallBMRB()

fetchallBMRB<-function(){
  outdat<-read.csv('http://www.bmrb.wisc.edu/ftp/pub/bmrb/relational_tables/nmr-star3.1/Atom_chem_shift.csv', header = T)
  return(outdat)
}


