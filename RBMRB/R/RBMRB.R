#'Downloads chemical shift data from BMRB for a given BMRB entry/list of BMRB entries
#'
#'@param BMRBidlist ==> sinlge BMRB ID / list of BMRB IDs in csv format
#'@return all available chemical shift data in R data frame
#'@export fetchBMRB
#'@examples
#'df<-fetchBMRB('15060')
#'df<-fetchBMRB('15060,15070,8898,99')
fetchBMRB<-function(BMRBidlist){
  bmrb_apiurl_json<<-"http://webapi.bmrb.wisc.edu/v0.4/jsonrpc"
  query=rjson::toJSON(list(method='loop',jsonrpc='2.0',params=list(ids=BMRBidlist,keys=list('_Atom_chem_shift')),id=1))
  rawdata<-httr::POST(bmrb_apiurl_json,encode='json',body=query)
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
  cs_data$Val<-as.numeric(cs_data$Val)
  cs_data$Val_err<-as.numeric(cs_data$Val_err)

  }
  else{
    warning('Entry not found')
    cs_data<-NA
  }
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
#'@param atom ==> atom name like CA,CB2
#'@param db ==> macromolecules, metabolomics
#'@return list of all atom chemical shifts for all BMRB entries as a R data frame
#'@export fetchCSList
#'@examples
#'df<-fetchCSList('CB2','macromolecules')
#'df<-fetchCSList('C1','metabolomics')
fetchCSList<-function(atom,db){
  bmrb_api<<-"http://webapi.bmrb.wisc.edu/"
  raw_data<-httr::GET(bmrb_api,path=paste0("/v0.4/rest/chemical_shifts/",atom,"/",db))
  dat<-httr::content(raw_data,'parsed')
  if (length(dat$data)==0){
    warning('Atom or db wrong')
    dat_frame<-NA
  }else{
  dat_tab<-data.table::as.data.table(dat$data)
  dat_frame<-as.data.frame(data.table::data.table(t(dat_tab)))
  for (name in names(dat_frame)){dat_frame[[name]]<-unlist(as.character(dat_frame[[name]]))}
  dat_tags<-as.data.frame(data.table::data.table(t(data.table::as.data.table(dat$columns))))
  for (i in 1:length(dat_tags$V1)){dat_tags$V1[i]<-strsplit(dat_tags$V1[i],"[.]")[[1]][2]}
  colnames(dat_frame)<-dat_tags$V1
  dat_frame$Val<-as.numeric(dat_frame$Val)
  dat_frame$Val_err<-as.numeric(dat_frame$Val_err)
  }
  return(dat_frame)
}


