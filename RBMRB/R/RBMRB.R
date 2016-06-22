
bmrb_api<<-"http://webapi.bmrb.wisc.edu/v0.4/jsonrpc"

#'Downloads chemical shift data from BMRB for a given BMRB entry/list of BMRB entries
#'
#'@param BMRBidlist ==> sinlge BMRB ID / list of BMRB IDs in csv format
#'@return all available chemical shift data in R data frame
#'@export fetchBMRB
#'@examples
#'df<-fetchBMRB('15060')
#'df<-fetchBMRB('15060,15070,8898,99')

fetchBMRB<-function(BMRBidlist){
  query=rjson::toJSON(list(method='loop',jsonrpc='2.0',params=list(ids=BMRBidlist,keys=list('_Atom_chem_shift')),id=1))
  rawdata<-httr::POST(bmrb_api,encode='json',body=query)
  c<-rjson::fromJSON(httr::content(rawdata,'text'))
  if (length(c$result)!=0){
  for (x in c$result){
    for (y in x$`_Atom_chem_shift`){
      csdata<-data.table::as.data.table(y$data)
      cstags<-as.data.frame(data.table::as.data.table(y$tags))$V1
      if (exists('cs_data')){
        cs_data<-rbind(cs_data,as.data.frame(data.table::data.table(t(csdata))))
      }else{
        cs_data<-as.data.frame(data.table::data.table(t(csdata)))

      }
    }
  }
  colnames(cs_data)<-cstags

  }
  else{
    warning('Entry not found')
  }
#   if (nchar(c)>5){
#     d<-gsub("\\]","",gsub("\\[","",c))
#     d2<-gsub("\n\"Entry_ID\",\"id\",\"Entity_ID\",\"Comp_index_ID\",\"Comp_ID\",\"Atom_ID\",\"Atom_type\",\"Val\",\"Val_err\",\"Ambiguity_code\",\"Assigned_chem_shift_list_ID\"","",d)
#     t<-read.table(textConnection(d2),sep="\n")
#     outdata<-reshape2::colsplit(t$V1,",",names=c('x','BMRB_ID','Entry_ID','Entity_ID','Comp_index_ID','Comp_ID','Atom_ID','Atom_type','Chemical_shift','err','Ambiguity_code','Assigned_chem_shift_list_ID'))
#     outdata$x<-NULL
#   }
#   else{
#     cat("Invalid BMRB ID")
#     outdata<-NA
#   }
  return (cs_data)
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
  names(shiftH)[names(shiftH)=="Val"]<-"H"
  shiftN<-subset(csdf,Atom_ID=="N")
  names(shiftN)[names(shiftN)=="Val"]<-"N"
  shiftHN<-merge(shiftH,shiftN,by=c('Entry_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  outdat<-shiftHN[,c("Entry_ID","Comp_index_ID","Assigned_chem_shift_list_ID","Comp_ID.x","Comp_ID.y","H","N")]
  names(outdat)[names(outdat)=="Comp_ID.x"]<-"Comp_ID_H"
  names(outdat)[names(outdat)=="Comp_ID.y"]<-"Comp_ID_N"
  return(outdat)
}

#'Reads the full chemical shift csv frile from BMRB ftp site. Not recommended unless you have high speed internet connection
#'
#'@param Specify the csv file downloaded from http://www.bmrb.wisc.edu/ftp/pub/bmrb/relational_tables/nmr-star3.1/Atom_chem_shift.csv  with full path;May leave it empty to fetch the file directly from internet
#'@return list of all atom chemical shifts for all BMRB entries as a R data frame
#'@export fetchallBMRB
#'@examples
#'df<-fetchallBMRB()
#'df<-fetchallBMRB('~/Downloads/Atom_chem_shift.csv')
fetchallBMRB<-function(csvpath='http://www.bmrb.wisc.edu/ftp/pub/bmrb/relational_tables/nmr-star3.1/Atom_chem_shift.csv'){
  outdat<-read.csv(csvpath, header = T)
  return(outdat)
}





