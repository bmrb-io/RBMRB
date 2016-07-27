#'NMR Chemical shifts data in R data frame
#'
#'Downloads NMR chemical shift data from BMRB(www.bmrb.wisc.edu) for a given Entry ID or list of Entry IDs
#'@param BMRBidlist ==> sinlge BMRB ID / list of BMRB IDs in csv format exampe:c(17074,17076,17077)
#'@return all available chemical shift data in R data frame
#'@export fetch_entry_chemical_shifts
#'@examples
#'df<-fetch_entry_chemical_shifts(15060)
#'df<-fetch_entry_chemical_shifts(c(17074,17076,17077))
fetch_entry_chemical_shifts<-function(BMRBidlist){
  bmrb_apiurl_json<-"http://webapi.bmrb.wisc.edu/v0.4/jsonrpc"
  query=rjson::toJSON(list(method='loop',jsonrpc='2.0',params=list(ids=BMRBidlist,keys=list('_Atom_chem_shift')),id=1))
  rawdata<-httr::POST(bmrb_apiurl_json,encode='json',body=query)
  c<-rjson::fromJSON(httr::content(rawdata,'text'))
  if (length(c$result)!=0){
  for (x in c$result){
    for (y in x$`_Atom_chem_shift`){
      if (is.null(y)){
        warnings('No Chemical shift list found')}
      else{
        csdata<-data.table::as.data.table(y$data)
        cstags<-as.data.frame(data.table::as.data.table(y$tags))$V1
        if (exists('cs_data')){
          cs_data<-rbind(cs_data,as.data.frame(data.table::data.table(t(csdata))))}
        else{
          cs_data<-as.data.frame(data.table::data.table(t(csdata)))}
        }
      }
    }
    if (exists('cs_data')){
      colnames(cs_data)<-cstags
      cs_data$Val<-as.numeric(cs_data$Val)
      cs_data$Val_err<-as.numeric(cs_data$Val_err)}
    else{
      warning('No data')
      cs_data<-NA}
    }
  else{
    warning('Entry not found')
    cs_data<-NA
  }
  return (cs_data)
}


#'Reformats chemical shift dataframe into H1-N15 HSQC dataframe
#'
#'Reformats the output dataframe from fetch_entry_chemical_shifts() into a simple dataframe contains only H1 and N15 chemical shifts, which could be used to simulate H1-N15 HSQC spectra
#'@param csdf ==> chemical shift data frame from fetch_entry_chemical_shift
#'@return 1H-N15 chemical shift list on the same row combined using comp index ID and bmrb ID
#'@export convert_cs_to_n15hsqc
#'@examples
#'df<-fetch_entry_chemical_shifts(15060)
#'hsqc<-convert_cs_to_n15hsqc(df)
convert_cs_to_n15hsqc<-function(csdf){
  if (all(is.na(csdf))){
    warning('No data')
    outdat<-NA
  }else{
  shiftH<-subset(csdf,Atom_ID=="H")
  names(shiftH)[names(shiftH)=="Val"]<-"H"
  shiftN<-subset(csdf,Atom_ID=="N")
  names(shiftN)[names(shiftN)=="Val"]<-"N"
  shiftHN<-merge(shiftH,shiftN,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  outdat<-shiftHN[,c("Entry_ID","Comp_index_ID","Entity_ID","Assigned_chem_shift_list_ID","Comp_ID.x","Comp_ID.y","H","N")]
  names(outdat)[names(outdat)=="Comp_ID.x"]<-"Comp_ID_H"
  names(outdat)[names(outdat)=="Comp_ID.y"]<-"Comp_ID_N"
  }
  return(outdat)
}

#'NMR Chemical shifts list
#'
#'Downloads the full list of chemical shifts from BMRB(www.bmrb.wisc.edu) macromolecular/metabolomics database
#'@param atom ==> atom name like CA,CB2
#'@param db ==> macromolecules, metabolomics
#'@return list of all atom chemical shifts for all BMRB entries as a R data frame
#'@export fetch_atom_chemical_shifts
#'@examples
#'df<-fetch_atom_chemical_shifts('CB2','macromolecules')
#'df<-fetch_atom_chemical_shifts('C1','metabolomics')
fetch_atom_chemical_shifts<-function(atom,db='macromolecules'){
  bmrb_api<-"http://webapi.bmrb.wisc.edu/"
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

#'Simulate H1-N15 HSQC spectra
#'
#'Simulates H1-N15 HSQC spectra directly from BMRB(www.bmrb.wisc.edu) database. Default plot type will be 'scatter'.Peaks from different spectra(entries) can be connected based on residue numbers by specifying plot type as 'line'
#'@param idlist ==> list of bmrb ids c(17074,17076,17077)
#'@param type ==> scatter/line default=scatter
#'@return plot object
#'@export simulate_n15hsqc
#'@examples
#'plot_hsqc<-simulate_n15hsqc(c(17074,17076,17077))
#'plot_hsqc<-simulate_n15hsqc(18857,'line')
simulate_n15hsqc<-function(idlist,type='scatter'){
  cs_data<-fetch_entry_chemical_shifts(idlist)
  hsqc_data<-convert_cs_to_n15hsqc(cs_data)
  if (all(is.na(hsqc_data))){
    warning('No data')
    plt2<-NA}
  else{
  hsqc_data$Info=NA
  hsqc_data$Info=paste(hsqc_data$Comp_index_ID,hsqc_data$Entity_ID,hsqc_data$Comp_ID_H,hsqc_data$Assigned_chem_shift_list_ID,sep=",")
  if (type=='scatter'){
    plt<-ggplot2::ggplot(hsqc_data)+
    ggplot2::geom_point(ggplot2::aes(x=H,y=N,color=Entry_ID,label=Info))#+
    #ggplot2::geom_line(ggplot2::aes(x=H,y=N,color=Comp_index_ID,label=key)) #+ theme(legend.position="none")

  } else {
    plt<-ggplot2::ggplot(hsqc_data)+
    ggplot2::geom_line(ggplot2::aes(x=H,y=N,group=Comp_index_ID,label=Info))+ #+ theme(legend.position="none")
    ggplot2::geom_point(ggplot2::aes(x=H,y=N,color=Entry_ID,label=Info))
  }
  plt2<-plotly::plotly_build(plt)
  plt2$layout$annotations=F
  #plt2$layout$showlegend=F
  plt2$layout$xaxis$autorange = "reversed"
  plt2$layout$yaxis$autorange = "reversed"
  }
  return(plt2)
}

