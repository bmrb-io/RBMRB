#'Imports NMR Chemical Shift data from BMRB into R data frame
#'
#'Downloads NMR chemical shift data from BMRB database for a given Entry ID or list of Entry IDs
#'@param BMRBidlist sinlge BMRB ID (or) list of BMRB IDs in csv format
#'For metabolomics entries entry id should have 'bmse' prefix example: c('bmse000034','bmse000035','bmse000036')
#'@return R data frame that contains relevant NMR chemical shifts from BMRB database
#'@export fetch_entry_chemical_shifts
#'@examples
#'df<-fetch_entry_chemical_shifts(15060)
#'# Downloads NMR chemical shifts of a single entry from BMRB
#'df<-fetch_entry_chemical_shifts(c(17074,17076,17077))
#'# Downloads NMR chemical shifts of multiple entries from BMRB
#'df<-fetch_entry_chemical_shifts(c('bmse000034','bmse000035','bmse000036'))
#'# Downloads data from BMRB metabolomics database
#'@seealso \code{\link{fetch_atom_chemical_shifts}}
fetch_entry_chemical_shifts<-function(BMRBidlist){
  bmrb_apiurl_json<-"http://webapi.bmrb.wisc.edu/v1/jsonrpc"
  query=rjson::toJSON(list(method='loop',jsonrpc='2.0',params=list(ids=BMRBidlist,keys=list('_Atom_chem_shift')),id=1))
  rawdata<-httr::POST(bmrb_apiurl_json,encode='json',body=query)
  c<-rjson::fromJSON(httr::content(rawdata,'text',encoding = 'UTF-8'))
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
      cs_data$Val<-suppressWarnings(as.numeric(cs_data$Val))
      cs_data$Val_err<-suppressWarnings(as.numeric(cs_data$Val_err))}
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


#'Reformats chemical shift dataframe for easy plotting
#'
#'Reformats the output dataframe from \link{fetch_entry_chemical_shifts} into a simple dataframe that contains algorithmically combined proton and nitrogen chemical shifts in two columns.
#'This will be helpful to plot 1H-15N HSQC(Hetronuclear Single Quantum Coherence) spectrum.
#'@param csdf Chemical shift data frame from \link{fetch_entry_chemical_shifts}
#'@return R data frame that contains proton and nitrogen chemical shifts in two columns for each residue
#'@export convert_cs_to_n15hsqc
#'@examples
#'df<-fetch_entry_chemical_shifts(15060)
#'#Downloads the chemical shift data from BMRB
#'hsqc<-convert_cs_to_n15hsqc(df)
#'#Reformats for easy plotting
#'@seealso \code{\link{convert_cs_to_c13hsqc}} and \code{\link{convert_cs_to_tocsy}}
convert_cs_to_n15hsqc<-function(csdf){
  with(csdf,{
  if (all(is.na(csdf))){
    warning('No data')
    outdat<-NA
  }else{

  #all amide proton nitrogen
  shiftH<-subset(csdf,Atom_ID=="H")
  names(shiftH)[names(shiftH)=="Val"]<-"H"
  shiftN<-subset(csdf,Atom_ID=="N")
  names(shiftN)[names(shiftN)=="Val"]<-"N"
  shiftHNt<-merge(shiftH,shiftN,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))

  #ARG side chain
  shiftHH11<-subset(csdf,Atom_ID=="HH11")
  names(shiftHH11)[names(shiftHH11)=="Val"]<-"H"
  shiftHH12<-subset(csdf,Atom_ID=="HH12")
  names(shiftHH12)[names(shiftHH12)=="Val"]<-"H"
  shiftHH21<-subset(csdf,Atom_ID=="HH21")
  names(shiftHH21)[names(shiftHH21)=="Val"]<-"H"
  shiftHH22<-subset(csdf,Atom_ID=="HH22")
  names(shiftHH22)[names(shiftHH22)=="Val"]<-"H"
  shiftHE<-subset(csdf,Atom_ID=="HE")
  names(shiftHE)[names(shiftHE)=="Val"]<-"H"
  shiftNE<-subset(csdf,Atom_ID=="NE")
  names(shiftNE)[names(shiftNE)=="Val"]<-"N"
  shiftNH1<-subset(csdf,Atom_ID=="NH1")
  names(shiftNH1)[names(shiftNH1)=="Val"]<-"N"
  shiftNH2<-subset(csdf,Atom_ID=="NH2")
  names(shiftNH2)[names(shiftNH2)=="Val"]<-"N"
  shiftNEHE=merge(shiftHE,shiftNE,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNEHE)))){shiftNEHE$Comp_index_ID=paste(shiftNEHE$Comp_index_ID,"HE")}
  shiftNHH11=merge(shiftHH11,shiftNH1,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNHH11)))){shiftNHH11$Comp_index_ID=paste(shiftNHH11$Comp_index_ID,"HH11")}
  shiftNHH12=merge(shiftHH12,shiftNH1,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNHH12)))){shiftNHH12$Comp_index_ID=paste(shiftNHH12$Comp_index_ID,"HH12")}
  shiftNHH21=merge(shiftHH21,shiftNH2,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNHH21)))){shiftNHH21$Comp_index_ID=paste(shiftNHH21$Comp_index_ID,"HH21")}
  shiftNHH22=merge(shiftHH22,shiftNH2,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNHH22)))){shiftNHH22$Comp_index_ID=paste(shiftNHH22$Comp_index_ID,"HH22")}


  #GLN sidechain
  shiftHE21<-subset(csdf,Atom_ID=="HE21")
  names(shiftHE21)[names(shiftHE21)=="Val"]<-"H"
  shiftHE22<-subset(csdf,Atom_ID=="HE22")
  names(shiftHE22)[names(shiftHE22)=="Val"]<-"H"
  shiftNE<-subset(csdf,Atom_ID=="NE2")
  names(shiftNE)[names(shiftNE)=="Val"]<-"N"
  shiftNEHE21=merge(shiftHE21,shiftNE,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNEHE21)))){shiftNEHE21$Comp_index_ID=paste(shiftNEHE21$Comp_index_ID,"HE21")}
  shiftNEHE22=merge(shiftHE22,shiftNE,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNEHE22)))){shiftNEHE22$Comp_index_ID=paste(shiftNEHE22$Comp_index_ID,"HE22")}

  #ASN sidechain
  shiftHD21<-subset(csdf,Atom_ID=="HD21")
  names(shiftHD21)[names(shiftHD21)=="Val"]<-"H"
  shiftHD22<-subset(csdf,Atom_ID=="HD22")
  names(shiftHD22)[names(shiftHD22)=="Val"]<-"H"
  shiftND<-subset(csdf,Atom_ID=="ND2")
  names(shiftND)[names(shiftND)=="Val"]<-"N"
  shiftNDHD21=merge(shiftHD21,shiftND,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNDHD21)))){shiftNDHD21$Comp_index_ID=paste(shiftNDHD21$Comp_index_ID,"HD21")}
  shiftNDHD22=merge(shiftHD22,shiftND,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNDHD22)))){shiftNDHD22$Comp_index_ID=paste(shiftNDHD22$Comp_index_ID,"HD22")}

  #HIS sidechain
  shiftHD1<-subset(csdf,Atom_ID=="HD1")
  names(shiftHD1)[names(shiftHD1)=="Val"]<-"H"
  shiftND1<-subset(csdf,Atom_ID=="ND1")
  names(shiftND1)[names(shiftND1)=="Val"]<-"N"
  shiftND1HD1=merge(shiftHD1,shiftND1,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftND1HD1)))){shiftND1HD1$Comp_index_ID=paste(shiftND1HD1$Comp_index_ID,"HD1")}
  shiftHE2<-subset(csdf,Atom_ID=="HE2")
  names(shiftHE2)[names(shiftHE2)=="Val"]<-"H"
  shiftNE2<-subset(csdf,Atom_ID=="NE2")
  names(shiftNE2)[names(shiftNE2)=="Val"]<-"N"
  shiftNE2HE2=merge(shiftHE2,shiftNE2,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNE2HE2)))){shiftNE2HE2$Comp_index_ID=paste(shiftNE2HE2$Comp_index_ID,"HE2")}

  #TRP sidechain
  shiftHE1<-subset(csdf,Atom_ID=="HE1")
  names(shiftHE1)[names(shiftHE1)=="Val"]<-"H"
  shiftNE1<-subset(csdf,Atom_ID=="NE1")
  names(shiftNE1)[names(shiftNE1)=="Val"]<-"N"
  shiftNE1HE1=merge(shiftHE1,shiftNE1,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNE1HE1)))){shiftNE1HE1$Comp_index_ID=paste(shiftNE1HE1$Comp_index_ID,"HE1")}


  #LYS sidechain
  shiftHZ<-subset(csdf,Atom_ID=="HZ")
  names(shiftHZ)[names(shiftHZ)=="Val"]<-"H"
  shiftHZ1<-subset(csdf,Atom_ID=="HZ1")
  names(shiftHZ1)[names(shiftHZ1)=="Val"]<-"H"
  shiftHZ2<-subset(csdf,Atom_ID=="HZ2")
  names(shiftHZ2)[names(shiftHZ2)=="Val"]<-"H"
  shiftHZ3<-subset(csdf,Atom_ID=="HZ3")
  names(shiftHZ3)[names(shiftHZ3)=="Val"]<-"H"
  shiftNZ<-subset(csdf,Atom_ID=="NZ")
  names(shiftNZ)[names(shiftNZ)=="Val"]<-"N"
  shiftNHZ=merge(shiftHZ,shiftNZ,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNHZ)))){shiftNHZ$Comp_index_ID=paste(shiftNHZ$Comp_index_ID,"HZ")}
  shiftNHZ1=merge(shiftHZ1,shiftNZ,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNHZ1)))){shiftNHZ1$Comp_index_ID=paste(shiftNHZ1$Comp_index_ID,"HZ1")}
  shiftNHZ2=merge(shiftHZ2,shiftNZ,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNHZ2)))){shiftNHZ2$Comp_index_ID=paste(shiftNHZ2$Comp_index_ID,"HZ2")}
  shiftNHZ3=merge(shiftHZ3,shiftNZ,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  if (!(all(is.na(shiftNHZ3)))){shiftNHZ3$Comp_index_ID=paste(shiftNHZ3$Comp_index_ID,"HZ3")}

  shiftHN<-rbind(rbind(rbind(rbind(rbind(rbind(rbind(rbind(rbind(rbind(rbind(rbind(rbind(rbind(rbind(rbind(shiftNEHE21,shiftNEHE22),shiftNDHD21),shiftNDHD22),shiftHNt),shiftND1HD1),shiftNE2HE2),shiftNEHE),shiftNHH11),shiftNHH12),shiftNHH21),shiftNHH22),shiftNHZ),shiftNHZ1),shiftNHZ2),shiftNHZ3),shiftNE1HE1)
  outdat<-shiftHN[,c("Entry_ID","Comp_index_ID","Entity_ID","Assigned_chem_shift_list_ID","Comp_ID.x","Comp_ID.y","H","N")]
  names(outdat)[names(outdat)=="Comp_ID.x"]<-"Comp_ID_H"
  names(outdat)[names(outdat)=="Comp_ID.y"]<-"Comp_ID_N"
  }
  return(outdat)
  })
}

#'Reformats chemical shift dataframe for easy plotting
#'
#'Reformats the output dataframe from \link{fetch_entry_chemical_shifts} into a simple dataframe that contains algorithmically combined proton shifts in two columns.
#'This will be helpful to plot TOCSY(TOtal Correlation SpectroscpY) spectrum
#'@param csdf chemical shift data frame from \link{fetch_entry_chemical_shifts}
#'@return R data frame that contains all possible combinations of proton chemical shifts in two columns
#'@export convert_cs_to_tocsy
#'@examples
#'df<-fetch_entry_chemical_shifts(15060)
#'# Downloads data from BMRB
#'tocsy<-convert_cs_to_tocsy(df)
#'# Reformats for easy plotting
#'@seealso \code{\link{convert_cs_to_c13hsqc}} and \code{\link{convert_cs_to_n15hsqc}}
convert_cs_to_tocsy<-function(csdf){
  with(csdf,{
  if (all(is.na(csdf))){
    warning('No data')
    outdat<-NA
  }else{
    cs_h<-subset(csdf,Atom_type=="H")
    outdat<-merge(cs_h,cs_h,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  }
  return(outdat)
  })
}

#'Reformats chemical shift dataframe for easy plotting
#'
#'Reformats the output dataframe from \link{fetch_entry_chemical_shifts} into a simple dataframe that contains proton and carbon chemical shifts in tow columns.
#'This will be helpful to plot 1H-13C HSQC(Hetronuclear Single Quantum Coherence) spectrum
#'@param csdf chemical shift data frame from \link{fetch_entry_chemical_shifts}
#'@return R data frame that contains proton and carbon chemical shifts in two columns for each residue
#'@export convert_cs_to_c13hsqc
#'@examples
#'df<-fetch_entry_chemical_shifts(15060)
#'# Downloads data from BMRB
#'hsqc<-convert_cs_to_c13hsqc(df)
#'# Reformats for easy plotting
#'@seealso \code{\link{convert_cs_to_n15hsqc}} and \code{\link{convert_cs_to_tocsy}}
convert_cs_to_c13hsqc<-function(csdf){
  with(csdf,{
  if (all(is.na(csdf))){
    warning('No data')
    outdat<-NA
  }else{
    shiftHA<-subset(csdf,Atom_ID=="HA")
    #names(shiftHA)[names(shiftHA)=="Val"]<-"HA"
    shiftHA2<-subset(csdf,Atom_ID=="HA2")
    #names(shiftHA2)[names(shiftHA2)=="Val"]<-"HA2"
    shiftHA3<-subset(csdf,Atom_ID=="HA3")
    #names(shiftHA3)[names(shiftHA3)=="Val"]<-"HA3"

    shiftHB<-subset(csdf,Atom_ID=="HB")
    #names(shiftHB)[names(shiftHB)=="Val"]<-"HB"
    shiftHB1<-subset(csdf,Atom_ID=="HB1")
    #names(shiftHB1)[names(shiftHB1)=="Val"]<-"HB1"
    shiftHB2<-subset(csdf,Atom_ID=="HB2")
    #names(shiftHB2)[names(shiftHB2)=="Val"]<-"HB2"
    shiftHB3<-subset(csdf,Atom_ID=="HB3")
    #names(shiftHB3)[names(shiftHB3)=="Val"]<-"HB3"

    shiftHG<-subset(csdf,Atom_ID=="HG")
    #names(shiftHG)[names(shiftHG)=="Val"]<-"HG"
    shiftHG2<-subset(csdf,Atom_ID=="HG2")
    #names(shiftHG2)[names(shiftHG2)=="Val"]<-"HG2"
    shiftHG3<-subset(csdf,Atom_ID=="HG3")
    #names(shiftHG3)[names(shiftHG3)=="Val"]<-"HG3"
    shiftHG11<-subset(csdf,Atom_ID=="HG11")
    #names(shiftHG11)[names(shiftHG11)=="Val"]<-"HG11"
    shiftHG12<-subset(csdf,Atom_ID=="HG12")
    #names(shiftHG12)[names(shiftHG12)=="Val"]<-"HG12"
    shiftHG13<-subset(csdf,Atom_ID=="HG13")
    #names(shiftHG13)[names(shiftHG13)=="Val"]<-"HG13"
    shiftHG21<-subset(csdf,Atom_ID=="HG21")
    #names(shiftHG21)[names(shiftHG21)=="Val"]<-"HG21"
    shiftHG22<-subset(csdf,Atom_ID=="HG22")
    #names(shiftHG22)[names(shiftHG22)=="Val"]<-"HG22"
    shiftHG23<-subset(csdf,Atom_ID=="HG23")
    #names(shiftHG23)[names(shiftHG23)=="Val"]<-"HG23"

    shiftHD1<-subset(csdf,Atom_ID=="HD1")
    #names(shiftHD1)[names(shiftHD1)=="Val"]<-"HD1"
    shiftHD2<-subset(csdf,Atom_ID=="HD2")
    #names(shiftHD2)[names(shiftHD2)=="Val"]<-"HD2"
    shiftHD3<-subset(csdf,Atom_ID=="HD3")
    #names(shiftHD3)[names(shiftHD3)=="Val"]<-"HD3"
    shiftHD11<-subset(csdf,Atom_ID=="HD11")
    #names(shiftHD11)[names(shiftHD11)=="Val"]<-"HD11"
    shiftHD12<-subset(csdf,Atom_ID=="HD12")
    #names(shiftHD12)[names(shiftHD12)=="Val"]<-"HD12"
    shiftHD13<-subset(csdf,Atom_ID=="HD13")
    #names(shiftHD13)[names(shiftHD13)=="Val"]<-"HD13"
    shiftHD21<-subset(csdf,Atom_ID=="HD21")
    #names(shiftHD21)[names(shiftHD21)=="Val"]<-"HD21"
    shiftHD22<-subset(csdf,Atom_ID=="HD22")
    #names(shiftHD22)[names(shiftHD22)=="Val"]<-"HD22"
    shiftHD23<-subset(csdf,Atom_ID=="HD23")
    #names(shiftHD23)[names(shiftHD23)=="Val"]<-"HD23"

    shiftHE21<-subset(csdf,Atom_ID=="HE21")
    #names(shiftHE21)[names(shiftHE21)=="Val"]<-"HE21"
    shiftHE22<-subset(csdf,Atom_ID=="HE22")
    #names(shiftHE22)[names(shiftHE22)=="Val"]<-"HE22"
    shiftHE1<-subset(csdf,Atom_ID=="HE1")
    #names(shiftHE1)[names(shiftHE1)=="Val"]<-"HE1"
    shiftHE2<-subset(csdf,Atom_ID=="HE2")
    #names(shiftHE2)[names(shiftHE2)=="Val"]<-"HE2"
    shiftHE3<-subset(csdf,Atom_ID=="HE3")
    #names(shiftHE3)[names(shiftHE3)=="Val"]<-"HE3"

    shiftHZ<-subset(csdf,Atom_ID=="HZ")
    #names(shiftHZ)[names(shiftHZ)=="Val"]<-"HZ"
    shiftHZ2<-subset(csdf,Atom_ID=="HZ2")
    #names(shiftHZ2)[names(shiftHZ2)=="Val"]<-"HZ2"
    shiftHZ3<-subset(csdf,Atom_ID=="HZ3")
    #names(shiftHZ3)[names(shiftHZ3)=="Val"]<-"HZ3"

    shiftCA<-subset(csdf,Atom_ID=="CA")
    #names(shiftCA)[names(shiftCA)=="Val"]<-"CA"

    shiftCB<-subset(csdf,Atom_ID=="CB")
    #names(shiftCB)[names(shiftCB)=="Val"]<-"CB"

    shiftCG<-subset(csdf,Atom_ID=="CG")
    #names(shiftCG)[names(shiftCG)=="Val"]<-"CG"
    shiftCG1<-subset(csdf,Atom_ID=="CG1")
    #names(shiftCG1)[names(shiftCG1)=="Val"]<-"CG1"
    shiftCG2<-subset(csdf,Atom_ID=="CG2")
    #names(shiftCG2)[names(shiftCG2)=="Val"]<-"CG2"

    shiftCD<-subset(csdf,Atom_ID=="CD")
    #names(shiftCD)[names(shiftCD)=="Val"]<-"CD"
    shiftCD1<-subset(csdf,Atom_ID=="CD1")
    #names(shiftCD1)[names(shiftCD1)=="Val"]<-"CD1"
    shiftCD2<-subset(csdf,Atom_ID=="CD2")
    #names(shiftCD2)[names(shiftCD2)=="Val"]<-"CD2"

    shiftCE<-subset(csdf,Atom_ID=="CE")
    #names(shiftCE)[names(shiftCE)=="Val"]<-"CE"
    shiftCE1<-subset(csdf,Atom_ID=="CE1")
    #names(shiftCE1)[names(shiftCE1)=="Val"]<-"CE1"
    shiftCE2<-subset(csdf,Atom_ID=="CE2")
    #names(shiftCE2)[names(shiftCE2)=="Val"]<-"CE2"
    shiftCE3<-subset(csdf,Atom_ID=="CE3")
    #names(shiftCE3)[names(shiftCE3)=="Val"]<-"CE3"

    shiftCZ<-subset(csdf,Atom_ID=="CZ")
    #names(shiftCZ)[names(shiftCZ)=="Val"]<-"CZ"
    shiftCZ2<-subset(csdf,Atom_ID=="CZ2")
    #names(shiftCZ2)[names(shiftCZ2)=="Val"]<-"CZ2"
    shiftCZ3<-subset(csdf,Atom_ID=="CZ3")
    #names(shiftCZ3)[names(shiftCZ3)=="Val"]<-"CZ3"

    shiftCAHA<-merge(shiftCA,shiftHA,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCAHA)))){shiftCAHA$Comp_index_ID=paste(shiftCAHA$Comp_index_ID,"HA")}
    shiftCAHA2<-merge(shiftCA,shiftHA2,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCAHA2)))){shiftCAHA2$Comp_index_ID=paste(shiftCAHA2$Comp_index_ID,"HA2")}
    shiftCAHA3<-merge(shiftCA,shiftHA3,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCAHA3)))){shiftCAHA3$Comp_index_ID=paste(shiftCAHA3$Comp_index_ID,"HA3")}
    shiftCAH<-rbind(rbind(shiftCAHA,shiftCAHA2),shiftCAHA3)
    #shiftCAH$type="CA"

    shiftCBHB<-merge(shiftCB,shiftHB,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCBHB)))){shiftCBHB$Comp_index_ID=paste(shiftCBHB$Comp_index_ID,"HB")}
    shiftCBHB1<-merge(shiftCB,shiftHB1,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCBHB1)))){shiftCBHB1$Comp_index_ID=paste(shiftCBHB1$Comp_index_ID,"HB1")}
    shiftCBHB2<-merge(shiftCB,shiftHB2,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCBHB2)))){shiftCBHB2$Comp_index_ID=paste(shiftCBHB2$Comp_index_ID,"HB2")}
    shiftCBHB3<-merge(shiftCB,shiftHB3,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCBHB3)))){shiftCBHB3$Comp_index_ID=paste(shiftCBHB3$Comp_index_ID,"HB3")}
    shiftCBH<-rbind(rbind(rbind(shiftCBHB,shiftCBHB1),shiftCBHB2),shiftCBHB3)
    #shiftCBH$type="CB"

    shiftCGHG<-merge(shiftCG,shiftHG,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCGHG)))){shiftCGHG$Comp_index_ID=paste(shiftCGHG$Comp_index_ID,"HG")}
    shiftCGHG2<-merge(shiftCG,shiftHG2,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCGHG2)))){shiftCGHG2$Comp_index_ID=paste(shiftCGHG2$Comp_index_ID,"HG2")}
    shiftCGHG3<-merge(shiftCG,shiftHG3,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCGHG3)))){shiftCGHG3$Comp_index_ID=paste(shiftCGHG3$Comp_index_ID,"HG3")}
    shiftCG1HG11<-merge(shiftCG1,shiftHG11,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCG1HG11)))){shiftCG1HG11$Comp_index_ID=paste(shiftCG1HG11$Comp_index_ID,"HG11")}
    shiftCG1HG12<-merge(shiftCG1,shiftHG12,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCG1HG12)))){shiftCG1HG12$Comp_index_ID=paste(shiftCG1HG12$Comp_index_ID,"HG12")}
    shiftCG1HG13<-merge(shiftCG1,shiftHG13,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCG1HG13)))){shiftCG1HG13$Comp_index_ID=paste(shiftCG1HG13$Comp_index_ID,"HG13")}
    shiftCG2HG21<-merge(shiftCG2,shiftHG21,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCG2HG21)))){shiftCG2HG21$Comp_index_ID=paste(shiftCG2HG21$Comp_index_ID,"HG21")}
    shiftCG2HG22<-merge(shiftCG2,shiftHG22,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCG2HG22)))){shiftCG2HG22$Comp_index_ID=paste(shiftCG2HG22$Comp_index_ID,"HG22")}
    shiftCG2HG23<-merge(shiftCG2,shiftHG23,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCG2HG23)))){shiftCG2HG23$Comp_index_ID=paste(shiftCG2HG23$Comp_index_ID,"HG23")}
    shiftCGH<-rbind(rbind(rbind(rbind(rbind(rbind(rbind(rbind(shiftCGHG,shiftCGHG2),shiftCGHG3),shiftCG1HG11),shiftCG1HG12),shiftCG1HG13),shiftCG2HG21),shiftCG2HG22),shiftCG2HG23)
    #shiftCGH$type="CG"

    shiftCDHD2<-merge(shiftCD,shiftHD2,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCDHD2)))){shiftCDHD2$Comp_index_ID=paste(shiftCDHD2$Comp_index_ID,"HD2")}
    shiftCDHD3<-merge(shiftCD,shiftHD3,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCDHD3)))){shiftCDHD3$Comp_index_ID=paste(shiftCDHD3$Comp_index_ID,"HD3")}
    shiftCD1HD1<-merge(shiftCD1,shiftHD1,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCD1HD1)))){shiftCD1HD1$Comp_index_ID=paste(shiftCD1HD1$Comp_index_ID,"HD1")}
    shiftCD2HD2<-merge(shiftCD2,shiftHD2,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCD2HD2)))){shiftCD2HD2$Comp_index_ID=paste(shiftCD2HD2$Comp_index_ID,"HD2")}
    shiftCD1HD11<-merge(shiftCD1,shiftHD11,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCD1HD11)))){shiftCD1HD11$Comp_index_ID=paste(shiftCD1HD11$Comp_index_ID,"HD11")}
    shiftCD1HD12<-merge(shiftCD1,shiftHD12,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCD1HD12)))){shiftCD1HD12$Comp_index_ID=paste(shiftCD1HD12$Comp_index_ID,"HD12")}
    shiftCD1HD13<-merge(shiftCD1,shiftHD13,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCD1HD13)))){shiftCD1HD13$Comp_index_ID=paste(shiftCD1HD13$Comp_index_ID,"HD13")}
    shiftCD2HD21<-merge(shiftCD2,shiftHD21,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCD2HD21)))){shiftCD2HD21$Comp_index_ID=paste(shiftCD2HD21$Comp_index_ID,"HD21")}
    shiftCD2HD22<-merge(shiftCD2,shiftHD22,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCD2HD22)))){shiftCD2HD22$Comp_index_ID=paste(shiftCD2HD22$Comp_index_ID,"HD22")}
    shiftCD2HD23<-merge(shiftCD2,shiftHD23,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCD2HD23)))){shiftCD2HD23$Comp_index_ID=paste(shiftCD2HD23$Comp_index_ID,"HD23")}
    shiftCDH<-rbind(rbind(rbind(rbind(rbind(rbind(rbind(rbind(rbind(shiftCDHD2,shiftCDHD3),shiftCD1HD1),shiftCD2HD2),shiftCD1HD11),shiftCD1HD12),shiftCD1HD13),shiftCD2HD21),shiftCD2HD22),shiftCD2HD23)
    #shiftCDH$type="CD"

    shiftCEHE2<-merge(shiftCE,shiftHE2,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCEHE2)))){shiftCEHE2$Comp_index_ID=paste(shiftCEHE2$Comp_index_ID,"HE2")}
    shiftCEHE3<-merge(shiftCE,shiftHE3,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCEHE3)))){shiftCEHE3$Comp_index_ID=paste(shiftCEHE3$Comp_index_ID,"HE3")}
    shiftCE1HE1<-merge(shiftCE1,shiftHE1,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCE1HE1)))){ shiftCE1HE1$Comp_index_ID=paste(shiftCE1HE1$Comp_index_ID,"HE1")}
    shiftCE2HE2<-merge(shiftCE2,shiftHE2,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCE2HE2)))){shiftCE2HE2$Comp_index_ID=paste(shiftCE2HE2$Comp_index_ID,"HE2")}
    shiftCE3HE3<-merge(shiftCE3,shiftHE3,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCE3HE3)))){shiftCE3HE3$Comp_index_ID=paste(shiftCE3HE3$Comp_index_ID,"HE3")}
    shiftCEH<-rbind(rbind(rbind(rbind(shiftCEHE2,shiftCEHE3),shiftCE1HE1),shiftCE2HE2),shiftCE3HE3)
    #shiftCEH$type="CE"

    shiftCZHZ<-merge(shiftCZ,shiftHZ,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCZHZ)))){shiftCZHZ$Comp_index_ID=paste(shiftCZHZ$Comp_index_ID,"HZ")}
    shiftCZ2HZ2<-merge(shiftCZ2,shiftHZ2,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCZ2HZ2)))){shiftCZ2HZ2$Comp_index_ID=paste(shiftCZ2HZ2$Comp_index_ID,"HZ2")}
    shiftCZ3HZ3<-merge(shiftCZ3,shiftHZ3,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
    if (!(all(is.na(shiftCZ3HZ3)))){shiftCZ3HZ3$Comp_index_ID=paste(shiftCZ3HZ3$Comp_index_ID,"HZ3")}
    shiftCZH<-rbind(rbind(shiftCZHZ,shiftCZ2HZ2),shiftCZ3HZ3)
    #shiftCZH$type="CZ"

    shift_pairs<-rbind(rbind(rbind(rbind(rbind(shiftCAH,shiftCBH),shiftCGH),shiftCDH),shiftCEH),shiftCZH)
    outdat<-shift_pairs[,c("Entry_ID","Comp_index_ID","Entity_ID","Assigned_chem_shift_list_ID","Comp_ID.x","Comp_ID.y","Atom_ID.x","Atom_ID.y","Val.x","Val.y")]
    names(outdat)[names(outdat)=="Comp_ID.x"]<-"Comp_ID_C"
    names(outdat)[names(outdat)=="Comp_ID.y"]<-"Comp_ID_H"
    names(outdat)[names(outdat)=="Atom_ID.x"]<-"Atom_ID_C"
    names(outdat)[names(outdat)=="Atom_ID.y"]<-"Atom_ID_H"
    names(outdat)[names(outdat)=="Val.x"]<-"C"
    names(outdat)[names(outdat)=="Val.y"]<-"H"
  }
  return(outdat)
  })
}

#'NMR Chemical shifts list for a given atom from BMRB
#'
#'Downloads the full list of chemical shifts from BMRB macromolecular/metabolomics database for a given atom
#'@param atom atom name in NMR-STAR atom nomenclature ; Example: CA,CB2
#'@param db macromolecules, metabolomics (optional, by default will fetch from macromolecules database)
#'@return R data frame that contains full chemical shift list for a given atom
#'@export fetch_atom_chemical_shifts
#'@examples
#'df<-fetch_atom_chemical_shifts('CG2','macromolecules')
#'# Downloads CB2 chemical shifts from macromolecules database at BMRB
#'df<-fetch_atom_chemical_shifts('C1','metabolomics')
#'# Downloads C1 chemical shifts from metabolomics database at BMRB
#'@seealso \code{\link{fetch_entry_chemical_shifts}},\code{\link{filter_residue}} and \code{\link{chemical_shift_corr}}
fetch_atom_chemical_shifts<-function(atom,db='macromolecules'){
  bmrb_api<-"http://webapi.bmrb.wisc.edu/"
  raw_data<-httr::GET(bmrb_api,path=paste0("/v1/rest/chemical_shifts/",atom,"/",db))
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
  dat_frame$Val<-suppressWarnings(as.numeric(dat_frame$Val))
  dat_frame$Val_err<-suppressWarnings(as.numeric(dat_frame$Val_err))
  }
  return(dat_frame)
}


#'NMR Chemical shifts list for a given atom from BMRB
#'
#'Downloads the full list of chemical shifts from BMRB macromolecular database for a given residue
#'@param res residue name in NMR-STAR atom nomenclature ; Example: ALA,GLY
#'@param atm atom name in NMR-STAR nomenclautre ; Example :CA,HB2
#'@return R data frame that contains full chemical shift list for a given atom
#'@export fetch_res_chemical_shifts
#'@examples
#'df<-fetch_res_chemical_shifts('ALA','CA')
#'# Downloads all atom chemical shifts of ALA from macromolecules database at BMRB
#'@seealso \code{\link{fetch_entry_chemical_shifts}},\code{\link{filter_residue}} and \code{\link{chemical_shift_corr}}
fetch_res_chemical_shifts<-function(res='*',atm='*'){
  res=gsub('[*]','%',res)
  atm=gsub('[*]','%',atm)
  if (res=='*'){res="%"}
  if (atm=='*'){atm="%"}
  bmrb_apiurl_json<-"http://webapi.bmrb.wisc.edu/v1/jsonrpc"
  query=rjson::toJSON(list(method='select',
                           jsonrpc='2.0',
                           params=list(database='macromolecules',
                                       query=list(where=list(Comp_ID=res,Atom_ID=atm),
                                                  select=list('Entry_ID','Entity_ID','Entity_assembly_ID','Comp_index_ID','Comp_ID','Atom_ID','Atom_type','Val','Val_err','Assigned_chem_shift_list_ID'),
                                                  hash='false',
                                                  from='Atom_chem_shift')),
                           id=1))
  rawdata<-httr::POST(bmrb_apiurl_json,encode='json',body=query)
  dat<-httr::content(rawdata,'parsed')
  if (length(dat$result)==0){
    warning('Atom or db wrong')
    dat_frame<-NA
  }else{
    dat_frame<-data.table::as.data.table(dat$result)
    names(dat_frame)[names(dat_frame)=="Atom_chem_shift.Entry_ID"]="Entry_ID"
    names(dat_frame)[names(dat_frame)=="Atom_chem_shift.Entity_ID"]="Entity_ID"
    names(dat_frame)[names(dat_frame)=="Atom_chem_shift.Entity_assembly_ID"]="Entity_assembly_ID"
    names(dat_frame)[names(dat_frame)=="Atom_chem_shift.Comp_index_ID"]="Comp_index_ID"
    names(dat_frame)[names(dat_frame)=="Atom_chem_shift.Comp_ID"]="Comp_ID"
    names(dat_frame)[names(dat_frame)=="Atom_chem_shift.Atom_ID"]="Atom_ID"
    names(dat_frame)[names(dat_frame)=="Atom_chem_shift.Atom_type"]="Atom_type"
    names(dat_frame)[names(dat_frame)=="Atom_chem_shift.Val"]="Val"
    names(dat_frame)[names(dat_frame)=="Atom_chem_shift.Val_err"]="Val_err"
    names(dat_frame)[names(dat_frame)=="Atom_chem_shift.Assigned_chem_shift_list_ID"]="Assigned_chem_shift_list_ID"
    dat_frame$Val=suppressWarnings(as.numeric(dat_frame$Val))
    dat_frame$Atom_ID=suppressWarnings(as.character(dat_frame$Atom_ID))
    dat_frame$Comp_ID=suppressWarnings(as.character(dat_frame$Comp_ID))
    }
  return(dat_frame)
}

#'Chemical shift histogram (or) density plot
#'
#'Plots the histogram (or) density of a chemical shift of given atom
#'@param res residue name in NMR-STAR atom nomenclature ; Example: ALA,GLY
#'@param atm atom name in NMR-STAR nomenclautre ; Example :CA,HB2
#'@param type count or density
#'@param interactive TRUE/FALSE default TRUE
#'@return R plot object
#'@export chemical_shift_hist
#'@examples
#'df<-chemical_shift_hist('ALA')
#'#plots the histogram of all atoms of ALA
#'@seealso \code{\link{fetch_entry_chemical_shifts}},\code{\link{filter_residue}} and \code{\link{chemical_shift_corr}}
chemical_shift_hist<-function(res='*',atm='*',type='count',bw=0.01,interactive=TRUE){
  cs_dat2<-fetch_res_chemical_shifts(res,atm)
  if (all(is.na(cs_dat2))){
    return(NA)
  }else{
    for (atom in unique(cs_dat2$Atom_ID)){
      ss<-sd(subset(cs_dat2,Atom_ID==atom)$Val)
      m<-mean(subset(cs_dat2,Atom_ID==atom)$Val)
      min_val<-m-(ss*10)
      max_val<-m+(ss*10)
      if (exists('cs_dat')){
        cs_dat<-rbind(cs_dat,subset(cs_dat2,Atom_ID==atom & Val>=min_val & Val<=max_val))
      }else{
        cs_dat<-subset(cs_dat2,Atom_ID==atom & Val>=min_val & Val<=max_val)
      }
    }
    if (type=='count'){
      plt<-ggplot2::ggplot(cs_dat)+
        ggplot2::geom_histogram(ggplot2::aes(x=Val,color=Atom_ID,fill=Atom_ID),binwidth=bw,position = 'identity',alpha=0.5)+
        ggplot2::xlab("Chemical shift")+
        ggplot2::ylab("Count")+
        ggplot2::labs(color="",fill="")
    } else{
      plt<-ggplot2::ggplot(cs_dat)+
        ggplot2::geom_density(ggplot2::aes(x=Val,color=Atom_ID,fill=Atom_ID),alpha=0.5,trim=T)+
        ggplot2::xlab("Chemical shift")+
        ggplot2::ylab("Density")+
        ggplot2::labs(color="",fill="")
    }
    if (interactive){
      plt2<-plotly::plotly_build(plt)
    }else{
      plt2<-plt
    }
  }
  return(plt2)
}


#'Simulates H1-N15 HSQC spectra for a given entry or list of entries from BMRB
#'
#'Simulates H1-N15 HSQC(Hetronuclear Single Quantum Coherence) spectra directly from BMRB database. Default plot type will be 'scatter'.Peaks from different spectra(entries) can be connected based on residue numbers by specifying plot type as 'line'.
#'By default it will generate interactive graphics using plotly library
#'@param idlist list of bmrb ids in csv
#'@param type scatter/line default=scatter
#'@param interactive TRUE/FALSE default=TRUE
#'@return R plot object
#'@export HSQC_15N
#'@examples
#'plot_hsqc<-HSQC_15N(c(17074,17076,17077))
#'#simulates N15-HSQC spectra for the given 3 entreis
#'plot_hsqc<-HSQC_15N(18857,'line')
#'#simulates the N15-HSQC spectra from many chemical shift lists from a single entry
#'plot_hsqc<-HSQC_15N(c(17074,17076,17077),interactive=FALSE)
#'#example for non interactive plots
#'@seealso \code{\link{HSQC_13C}} and \code{\link{TOCSY}}
HSQC_15N<-function(idlist,type='scatter',interactive=TRUE){
  cs_data<-fetch_entry_chemical_shifts(idlist)
  if (all(is.na(cs_data))){
    return(NA)
  }else{
  hsqc_data<-convert_cs_to_n15hsqc(cs_data)
  with(hsqc_data,{
  if (all(is.na(hsqc_data))){
    warning('No data')
    plt2<-NA}
  else{
  hsqc_data$Info=NA
  hsqc_data$Info=paste(hsqc_data$Comp_index_ID,hsqc_data$Entity_ID,hsqc_data$Comp_ID_H,hsqc_data$Assigned_chem_shift_list_ID,sep=",")
  if (type=='scatter'){
    if (length(idlist)>1){
      if (interactive){
    plt<-ggplot2::ggplot(hsqc_data)+
    ggplot2::geom_point(ggplot2::aes(x=H,y=N,color=Entry_ID,label=Info))+ggplot2::labs(color="")
      #ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()
      }else{
        plt<-ggplot2::ggplot(hsqc_data)+
          ggplot2::geom_point(ggplot2::aes(x=H,y=N,color=Entry_ID,label=Info))+
          ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+ggplot2::labs(color="")
      }
    }else{
      if (interactive){
      plt<-ggplot2::ggplot(hsqc_data)+
        ggplot2::geom_point(ggplot2::aes(x=H,y=N,color=Comp_ID_H,label=Info))+ggplot2::labs(color="")
        #ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()
      }else{
        plt<-ggplot2::ggplot(hsqc_data)+
          ggplot2::geom_point(ggplot2::aes(x=H,y=N,color=Comp_ID_H,label=Info))+
          ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+ggplot2::labs(color="")
      }
    }

  } else {
    if (interactive){
    plt<-ggplot2::ggplot(hsqc_data)+
    ggplot2::geom_line(ggplot2::aes(x=H,y=N,group=Comp_index_ID,label=Info))+
    ggplot2::geom_point(ggplot2::aes(x=H,y=N,color=Entry_ID,label=Info))+ggplot2::labs(color="")
     # ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()
    }else{
      plt<-ggplot2::ggplot(hsqc_data)+
        ggplot2::geom_line(ggplot2::aes(x=H,y=N,group=Comp_index_ID,label=Info))+
        ggplot2::geom_point(ggplot2::aes(x=H,y=N,color=Entry_ID,label=Info))+
        ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+ggplot2::labs(color="")
    }
  }
  if (interactive){
    plt2<-plotly::plotly_build(plt)
    plt2$x$layout$xaxis$autorange = "reversed"
    plt2$x$layout$yaxis$autorange = "reversed"}
  else{
    plt2<-plt
  }
  }
  return(plt2)})}
}

#'Simulates H1-C13 HSQC spectra for a given entry or list of entries from BMRB
#'
#'Simulates H1-C13 HSQC(Hetronuclear Single Quantum Coherence) spectra directly from BMRB database. By default it will generate interactive graphics using plotly library
#'@param idlist list of bmrb ids in csv
#'@param type scatter/line default=scatter
#'@param interactive TRUE/FALSE default=TRUE
#'@return R plot object
#'@export HSQC_13C
#'@examples
#'plot_hsqc<-HSQC_13C(c(17074,17076,17077))
#'#Simulates C13-HSQC spectra form the given list of entries
#'plot_hsqc<-HSQC_13C(c(17074,17076,17077),'line')
#'#Simulates C13-HSQC and connects the peaks with same sequence number
#'plot_hsqc<-HSQC_13C(c(17074,17076,17077),interactive=FALSE)
#'#Example for non interactive plot
#'@seealso \code{\link{HSQC_15N}} and \code{\link{TOCSY}}
#'
#'
HSQC_13C<-function(idlist,type='scatter',interactive=TRUE){
  cs_data<-fetch_entry_chemical_shifts(idlist)
  if (all(is.na(cs_data))){
    return(NA)
  }else{
    hsqc_data<-convert_cs_to_c13hsqc(cs_data)
    with(hsqc_data,{
      if (all(is.na(hsqc_data))){
        warning('No data')
        plt2<-NA}
      else{
        hsqc_data$Info=NA
        hsqc_data$Info=paste(hsqc_data$Comp_index_ID,hsqc_data$Entity_ID,hsqc_data$Comp_ID_H,hsqc_data$Assigned_chem_shift_list_ID,sep=",")
        if (type=='scatter'){
          if (length(idlist)>1){
            if (interactive){
              plt<-ggplot2::ggplot(hsqc_data)+
                ggplot2::geom_point(ggplot2::aes(x=H,y=C,color=Entry_ID,label=Info))+ggplot2::labs(color="")
              #ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()
            }else{
              plt<-ggplot2::ggplot(hsqc_data)+
                ggplot2::geom_point(ggplot2::aes(x=H,y=C,color=Entry_ID,label=Info))+
                ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+ggplot2::labs(color="")
            }
          }else{
            if (interactive){
              plt<-ggplot2::ggplot(hsqc_data)+
                ggplot2::geom_point(ggplot2::aes(x=H,y=C,color=Comp_ID_H,label=Info))+ggplot2::labs(color="")
              #ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()
            }else{
              plt<-ggplot2::ggplot(hsqc_data)+
                ggplot2::geom_point(ggplot2::aes(x=H,y=C,color=Comp_ID_H,label=Info))+
                ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+ggplot2::labs(color="")
            }
          }

        } else {
          if (interactive){
            plt<-ggplot2::ggplot(hsqc_data)+
              ggplot2::geom_line(ggplot2::aes(x=H,y=C,group=Comp_index_ID,label=Info))+
              ggplot2::geom_point(ggplot2::aes(x=H,y=C,color=Entry_ID,label=Info))+ggplot2::labs(color="")
            # ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()
          }else{
            plt<-ggplot2::ggplot(hsqc_data)+
              ggplot2::geom_line(ggplot2::aes(x=H,y=C,group=Comp_index_ID,label=Info))+
              ggplot2::geom_point(ggplot2::aes(x=H,y=C,color=Entry_ID,label=Info))+
              ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+ggplot2::labs(color="")
          }
        }
        if (interactive){
          plt2<-plotly::plotly_build(plt)
          plt2$x$layout$xaxis$autorange = "reversed"
          plt2$x$layout$yaxis$autorange = "reversed"}
        else{
          plt2<-plt
        }
      }
      return(plt2)})}
}


#'Chemical shift correlation between any two atoms from a single residue
#'
#'Plots the correlated chemical shift distribution of any two atoms in a single residue for the 20 standard amino acids from BMRB database.
#'By default it will generate interactive graphics using plotly library
#'@param atom1 atom name in NMR-STAR nomenclature like CA,CB2
#'@param atom2 atom name in NMR_STAR nomenclature like HA,HB2
#'@param res residue name like ALA,GLY (optional by default includes all possible amino acids)
#'@param type 'c' for contour plot and 's' for scatter plot default 'c'.scatter plot will be slow and heavy for large data set
#'@param interactive TRUE/FALSE default=TRUE
#'@return plot object
#'@export chemical_shift_corr
#'@examples
#'plt<-chemical_shift_corr('HE21','HE22')
#'#plots the chemical shift distribution between HE21 and HE22
#'@seealso \code{\link{fetch_atom_chemical_shifts}}
chemical_shift_corr<-function(atom1,atom2,res=NA,type="c",interactive=TRUE){
  at1_cs<-fetch_atom_chemical_shifts(atom1)
  at2_cs<-fetch_atom_chemical_shifts(atom2)
  if (all(is.na(at1_cs)) | all(is.na(at2_cs))){
    warning('No data or not a valid atom name')
    plt2<-NA
    return(plt2)}
  else{
  with(at1_cs,{
    with(at2_cs, {
  if (is.na(res)){
    cs1<-filter_residue(at1_cs)
    cs2<-filter_residue(at2_cs)
  }
  else{
    cs1<-subset(at1_cs,Comp_ID==res)
    cs2<-subset(at2_cs,Comp_ID==res)
  }
  cs<-merge(cs1,cs2,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))[,c("Entry_ID","Comp_index_ID","Entity_ID","Assigned_chem_shift_list_ID","Comp_ID.x","Comp_ID.y","Atom_ID.x","Atom_ID.y","Val.x","Val.y")]
  names(cs)[names(cs)=="Comp_ID.x"]<-"Comp_ID_1"
  names(cs)[names(cs)=="Comp_ID.y"]<-"Comp_ID_2"
  names(cs)[names(cs)=="Atom_ID.x"]<-"Atom_ID_1"
  names(cs)[names(cs)=="Atom_ID.y"]<-"Atom_ID_2"
  xl1=mean(cs$Val.x)-5*sd(cs$Val.x)
  xl2=mean(cs$Val.x)+5*sd(cs$Val.x)
  yl1=mean(cs$Val.y)-5*sd(cs$Val.y)
  yl2=mean(cs$Val.y)+5*sd(cs$Val.y)
  cs<-subset(cs,(Val.x>xl1 & Val.x<xl2 & Val.y>yl1 & Val.y<yl2))
  #names(cs)[names(cs)=="Val.x"]<-atom1
  #names(cs)[names(cs)=="Val.y"]<-atom2
  cs$Info=paste(cs$Comp_index_ID,cs$Entity_ID,cs$Atom_ID_1,cs$Atom_ID_2,cs$Assigned_chem_shift_list_ID,sep=",")
  if (type=="c"){
    if (interactive){
  plt<-ggplot2::ggplot(cs)+
    #ggplot2::geom_density_2d(ggplot2::aes(x=eval(as.name(atom1)),y=eval(as.name(atom2)),color=Comp_ID_1))+
    #ggplot2::stat_density_2d(geom='polygon',ggplot2::aes(x=eval(as.name(atom1)),y=eval(as.name(atom2)),fill= Comp_ID_1),bins=500)+
    ggplot2::geom_density_2d(ggplot2::aes(x=Val.x,y=Val.y,color=Comp_ID_1),bins=100)+
    #ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+
    ggplot2::xlab(atom1)+
    ggplot2::ylab(atom2)+ggplot2::labs(color="")
    }else{
      plt<-ggplot2::ggplot(cs)+
        #ggplot2::geom_density_2d(ggplot2::aes(x=eval(as.name(atom1)),y=eval(as.name(atom2)),color=Comp_ID_1))+
        #ggplot2::stat_density_2d(geom='polygon',ggplot2::aes(x=eval(as.name(atom1)),y=eval(as.name(atom2)),fill= Comp_ID_1),bins=500)+
        ggplot2::geom_density_2d(ggplot2::aes(x=Val.x,y=Val.y,color=Comp_ID_1),bins=100)+
        ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+
        ggplot2::xlab(atom1)+
        ggplot2::ylab(atom2)+ggplot2::labs(color="")
    }
  }
  else{
    if (interactive){
    plt<-ggplot2::ggplot(cs)+
      ggplot2::geom_point(ggplot2::aes(x=Val.x,y=Val.y,color=Comp_ID_1))+
      #ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+
      ggplot2::xlab(atom1)+
      ggplot2::ylab(atom2)+ggplot2::labs(color="")
    }else{
      plt<-ggplot2::ggplot(cs)+
        ggplot2::geom_point(ggplot2::aes(x=Val.x,y=Val.y,color=Comp_ID_1))+
        ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+
        ggplot2::xlab(atom1)+
        ggplot2::ylab(atom2)+ggplot2::labs(color="")
    }
  }
  if (interactive){
    plt2<-plotly::plotly_build(plt)
    plt2$x$layout$xaxis$autorange = "reversed"
    plt2$x$layout$yaxis$autorange = "reversed"}
  else{
    plt2<-plt
  }
  return(plt2)
  })})
  }
}

#'Filter for standard 20 amino acids
#'
#'Filters out non standard amino acids using Comp_ID. The data frame should contain three letter anio acid code in COMP_ID column.
#'@param df data frame with amino acid information in Comp_ID column
#'@return R data frame that contains information from only standard 20 amino acids.
#'@export filter_residue
#'@examples
#'df<-filter_residue(fetch_atom_chemical_shifts("CG2"))
#'#Downloads all CG2 chemical shifts and removes non standard amino acids
#'@seealso \code{\link{fetch_atom_chemical_shifts}}
filter_residue<-function(df){
  if (is.data.frame(df)){
  with(df,{
    if ("Comp_ID" %in% colnames(df)){
  out_dat<-subset(df,Comp_ID=="ALA" |
                    Comp_ID=="ARG" |
                    Comp_ID=="ASP" |
                    Comp_ID=="ASN" |
                    Comp_ID=="CYS" |
                    Comp_ID=="GLU" |
                    Comp_ID=="GLN" |
                    Comp_ID=="GLY" |
                    Comp_ID=="HIS" |
                    Comp_ID=="ILE" |
                    Comp_ID=="LEU" |
                    Comp_ID=="LYS" |
                    Comp_ID=="MET" |
                    Comp_ID=="PHE" |
                    Comp_ID=="PRO" |
                    Comp_ID=="SER" |
                    Comp_ID=="THR" |
                    Comp_ID=="TRP" |
                    Comp_ID=="TYR" |
                    Comp_ID=="VAL")
  return(out_dat)
    }
  else{
    warning("Comp_ID column not found")
    out_dat<-NA
    return(out_dat)
  }
  })
  }
  else{
    warning("Argument is not a valid data frame")
    out_dat<-NA
    return(out_dat)
  }
}



#'Simulates TOCSY spectra for a given entry or a list of entries from BMRB
#'
#'Simulates TOCSY(TOtal Correlation SpectroscopY) spectra directly from BMRB database. By default it will generate interactive graphics using plotly library
#'@param idlist list of bmrb ids c(17074,17076,17077)
#'@param interactive TRUE/FALSE default=TRUE
#'@return plot object
#'@export TOCSY
#'@examples
#'plot_tocsy<-TOCSY(c(17074,17076,17077))
#'#Simulates TOCSY spectra for the given 3 entries
#'plot_tocsy<-TOCSY(c(17074,17076,17077),interactive=FALSE)
#'# Example to disable interactive plot feature
#'@seealso \code{\link{HSQC_15N}} and \code{\link{HSQC_13C}}
TOCSY<-function(idlist,interactive=TRUE){
  cs<-fetch_entry_chemical_shifts(idlist)
  with(cs,{cs_h<-subset(cs,Atom_type=="H")
  tocsy_data<-merge(cs_h,cs_h,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))
  with(tocsy_data,{
  if (all(is.na(tocsy_data))){
    warning('No data')
    plt2<-NA}
  else{
    tocsy_data$Info=NA
    tocsy_data$Info=paste(tocsy_data$Comp_index_ID,tocsy_data$Entity_ID,tocsy_data$Comp_ID.x,tocsy_data$Atom_ID.x,tocsy_data$Atom_ID.y,tocsy_data$Assigned_chem_shift_list_ID,sep=",")
    if (length(idlist)>1){
      if (interactive){
      plt<-ggplot2::ggplot(tocsy_data)+
        ggplot2::geom_point(ggplot2::aes(x=Val.x,y=Val.y,color=Entry_ID,label=Info))+
        #ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+
        ggplot2::xlab("H")+ggplot2::ylab("H")+ggplot2::labs(color="")
      }else{
        plt<-ggplot2::ggplot(tocsy_data)+
          ggplot2::geom_point(ggplot2::aes(x=Val.x,y=Val.y,color=Entry_ID,label=Info))+
          ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+
          ggplot2::xlab("H")+ggplot2::ylab("H")+ggplot2::labs(color="")
      }
    }else{
      if (interactive){
      plt<-ggplot2::ggplot(tocsy_data)+
        ggplot2::geom_point(ggplot2::aes(x=Val.x,y=Val.y,color=Comp_ID.x,label=Info))+
        #ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+
        ggplot2::xlab("H")+ggplot2::ylab("H")+ggplot2::labs(color="")
      }else{
        plt<-ggplot2::ggplot(tocsy_data)+
          ggplot2::geom_point(ggplot2::aes(x=Val.x,y=Val.y,color=Comp_ID.x,label=Info))+
          ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+ggplot2::labs(color="")+
          ggplot2::xlab("H")+ggplot2::ylab("H")
      }
    }
    if (interactive){
      plt2<-plotly::plotly_build(plt)
      plt2$x$layout$xaxis$autorange = "reversed"
      plt2$x$layout$yaxis$autorange = "reversed"}
    else{
      plt2<-plt
    }
  }
  return(plt2)
})})
}

