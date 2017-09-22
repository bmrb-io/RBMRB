#'Imports checmial shift table for a given entry id from BMRB data base
#'
#'Downloads NMR chemical shift data from BMRB database for a given Entry ID
#'@param ID sinlge BMRB ID; For macromolecule entries it is just a number without bmrb prefix (example: 15060);
#'For metabolomics entries it should have 'bmse' prefix (example: 'bmse000035')
#'@return R data frame that contains  Atom_chem_shift data for a given entry ID
#'@export fetch_entry_cs
#'@examples
#'df<-fetch_entry_cs(15060)
#'# Downloads NMR chemical shifts of the given entry from macromolecule database
#'df<-fetch_entry_cs('bmse000034')
#'# Downloads data from BMRB metabolomics database
#'@seealso \code{\link{fetch_entry_chemical_shifts}},\code{\link{fetch_atom_chemical_shifts}} and \code{\link{fetch_res_chemical_shifts}}
fetch_entry_cs<-function(ID){
  if (length(ID)>1){
    stop("This Function takes only one ID; Use fetch_entry_chemical_shifts for a list of IDs")
  }
  bmrb_apiurl<-paste0("http://webapi.bmrb.wisc.edu/v2/entry/",ID,"?loop=Atom_chem_shift")
  rawdata<-httr::GET(bmrb_apiurl,httr::add_headers(Application = "RBMRB V2.1.0"))
  c<-rjson::fromJSON(httr::content(rawdata,'text',encoding = 'UTF-8'))
  if (rawdata$status_code==200){
    for (y in c[[1]][1]$`Atom_chem_shift`){
      if (is.null(y)){
        warnings('No Chemical shift list found')}
      else{
        csdata<-data.table::as.data.table(y$data)
        cstags<-as.data.frame(data.table::as.data.table(y$tags))$V1
        csdata<-as.data.frame(data.table::data.table(t(csdata)))
        colnames(csdata)<-cstags
        if (!("Entry_ID" %in% colnames(csdata))){
          csdata$Entry_ID = makeRandomString()
        }
        if (!("Comp_index_ID" %in% colnames(csdata)) & ("Seq_ID" %in% colnames(csdata)) ){
          csdata$Comp_index_ID = csdata$Seq_ID
        }
        if (!("Entity_ID" %in% colnames(csdata))){
          csdata$Entity_ID=1
        }
        if (!("Assigned_chem_shift_list_ID" %in% colnames(csdata))){
          csdata$Assigned_chem_shift_list_ID=1
        }
        if (exists('cs_data') & exists('csdata')){
          if (length(colnames(cs_data)) != length(colnames(csdata))){
            common_col<-intersect(colnames(cs_data),colnames(csdata))
            cs_data<-subset(cs_data,select=common_col)
            csdata<-subset(csdata,select=common_col)
            warning("Entries have differnet columns;mismatch will be removed")
          }
          cs_data<-rbind(cs_data,csdata)}
        else{
          cs_data<-csdata}
      }
    }
    if (exists('cs_data')){
      cs_data$Val<-suppressWarnings(as.numeric(cs_data$Val))
      cs_data$Val_err<-suppressWarnings(as.numeric(cs_data$Val_err))
    }
    else{
      warning('No data')
      cs_data<-NA}
  }
  else{
    warning(c$error)
    cs_data<-NA
  }
  return (cs_data)
}




#'Imports chemical shift table for a given entry or list of entries from BMRB data base
#'
#'Downloads NMR chemical shift data from BMRB database for a given Entry ID or list of Entry IDs
#'@param IDlist sinlge BMRB ID (or) list of BMRB IDs in csv format; For macromolecule entries it is just a number without bmrb prefix (example: c(15060,15000,18867));
#'For metabolomics entries it should have 'bmse' prefix (example: c('bmse000035','bmse000035','bmse000036'))
#'@return R data frame that contains  Atom_chem_shift data for a given list of entries
#'@export fetch_entry_chemical_shifts
#'@examples
#'df<-fetch_entry_chemical_shifts(15060)
#'# Downloads NMR chemical shifts of a single entry from BMRB
#'df<-fetch_entry_chemical_shifts(c(17074,17076,17077))
#'# Downloads NMR chemical shifts of multiple entries from BMRB
#'df<-fetch_entry_chemical_shifts(c('bmse000034','bmse000035','bmse000036'))
#'# Downloads data from BMRB metabolomics database
#'@seealso \code{\link{fetch_atom_chemical_shifts}},\code{\link{fetch_entry_cs}} and \code{\link{fetch_res_chemical_shifts}}
fetch_entry_chemical_shifts<-function(IDlist){
  for (bid in IDlist){
    csdata<-fetch_entry_cs(bid)
    if (all(is.na(csdata))){
      warning("Entry not found")
    }
    else{
      if (exists('cs_data') & exists('csdata')){
        if (length(colnames(cs_data)) != length(colnames(csdata))){
          common_col<-intersect(colnames(cs_data),colnames(csdata))
          cs_data<-subset(cs_data,select=common_col)
          csdata<-subset(csdata,select=common_col)
          warning("Entries have differnet columns;mismatch will be removed")
        }
        cs_data<-rbind(cs_data,csdata)}
      else{
        cs_data<-csdata}
    }}
  if (exists('cs_data')==F){
    cs_data=NA
  }

  return(cs_data)
}

#'Generates random string of fixed length(for internal use in RBMRB)
#'
#'Local files may not have Entry_ID, in that case random Entry_ID is assigned using this function. It is an internal function used only by RBMRB package
makeRandomString <- function()
{
  n = 1
  lenght = 6
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    lenght, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}


#'Exports NMR-STAR file to BMRB API server
#'
#'Exports NMR-STAR file to BMRB API server for data visualization.  This function will return a tocken, which can be used like a pseudo BMRB ID. The tocken will expire after 7 days
#'@param filename filename with correct path
#'@return Temporary tocken to access the file
#'@export export_star_data
#'@examples
#'# ent_id <- export_star_data('/nmrdata/hpr.str')
#'# Exports hpr.str file to BMRB API server and gets a temporary tocken
#'@seealso \code{\link{fetch_atom_chemical_shifts}}, \code{\link{fetch_entry_chemical_shifts}} \code{\link{fetch_res_chemical_shifts}}
export_star_data<-function(filename){
  bmrb_apiurl<-"http://webapi.bmrb.wisc.edu/v2/entry/"
  query=rjson::toJSON(list(method='store',jsonrpc='2.0',params=list(data=readChar(filename, file.info(filename)$size)),id=1))
  rawdata<-httr::POST(bmrb_apiurl,body=readChar(filename, file.info(filename)$size),httr::add_headers(Application = "RBMRB V2.1.0"))
  c<-rjson::fromJSON(httr::content(rawdata,'text',encoding = 'UTF-8'))
  if (length(c)!=0){
    ent_id<-c$entry_id
    print(paste("Please note down the ID",ent_id,sep=":"))
    print(paste("ID will expire on ", as.POSIXct(c$expiration, origin = "1970-01-01"),sep=" "))
  }
  else{
    warning('Entry not found')
    ent_id<-NA
  }
  return (ent_id)
}


#'Imports all chemical shifts of a given atom from BMRB database
#'
#'Downloads the full chemical shift data from BMRB macromolecules/metabolomics database for a given atom
#'@param atom atom name in NMR-STAR atom nomenclature ; Example: CA,CB2; default * (all atoms)
#'@param db macromolecules, metabolomics (optional, by default will fetch from macromolecules database)
#'@return R data frame that contains full chemical shift list for a given atom
#'@export fetch_atom_chemical_shifts
#'@examples
#'#df<-fetch_atom_chemical_shifts('CG2','macromolecules')
#'# Downloads CB2 chemical shifts from macromolecules database at BMRB
#'#df<-fetch_atom_chemical_shifts('C1','metabolomics')
#'# Downloads C1 chemical shifts from metabolomics database at BMRB
#'@seealso \code{\link{fetch_entry_chemical_shifts}},\code{\link{fetch_res_chemical_shifts}},\code{\link{filter_residue}} and \code{\link{chem_shift_corr}} and \code{\link{atom_chem_shift_corr}}
fetch_atom_chemical_shifts<-function(atom="*",db='macromolecules'){
  bmrb_api<-paste0("http://webapi.bmrb.wisc.edu/v2/search/chemical_shifts?atom_id=",atom,"&database=",db)
  raw_data<-httr::GET(bmrb_api,httr::add_headers(Application = "RBMRB V2.1.0"))
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

#'Imports chemical shift data for a given amino acid/nucleic acid
#'
#'Downloads chemical shift data from BMRB macromolecular database for a given amino acid (or) nucleic acid. Optionally particular atom can be specified in the parameter
#'@param res residue name in NMR-STAR atom nomenclature ; Example: ALA,GLY ; default '*' (all residues)
#'@param atm atom name in NMR-STAR nomenclautre ; Example :CA,HB2; default * (all atoms)
#'@return R data frame that contains full chemical shift list for a given atom
#'@export fetch_res_chemical_shifts
#'@examples
#'#df<-fetch_res_chemical_shifts('GLY')
#'# Downloads chemical shift data of all atoms of GLY
#'#df<-fetch_res_chemical_shifts('ALA','CA')
#'# Downloads C alpha chemical shifts of ALA from macromolecules database at BMRB
#'@seealso \code{\link{fetch_atom_chemical_shifts}},\code{\link{filter_residue}} and \code{\link{chemical_shift_hist}}
fetch_res_chemical_shifts<-function(res='*',atm='*'){
  if (res=="*"){
    dat_frame<-fetch_atom_chemical_shifts(atm)
  }else{
  db='macromolecules'
  bmrb_api<-paste0("http://webapi.bmrb.wisc.edu/v2/search/chemical_shifts?comp_id=",res,"&atom_id=",atm,"&database=",db)
  raw_data<-httr::GET(bmrb_api,httr::add_headers(Application = "RBMRB V2.1.0"))
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
  }
  return(dat_frame)
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
#'Reformats the output dataframe from \link{fetch_entry_chemical_shifts} into a simple dataframe that contains proton and carbon chemical shifts in two columns.
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


#'Plots chemical shift distribution of all atoms of a given amino acid
#'
#'Plots the histogram (or) density of  chemical shift distribution of  all atoms of a given amino acid (or) nucleic acid from BMRB database.
#'@param res residue name in NMR-STAR atom nomenclature ; Example: ALA,GLY
#'@param type count ; other than count will assume density plot
#'@param cutoff values not with in the cutoff time standard deviation from both sides ofthe mean will be excluded from the plot;default value 8
#'@param interactive TRUE/FALSE default TRUE
#'@return R plot object
#'@export chemical_shift_hist_res
#'@examples
#'#plt<-chemical_shift_hist_res('ALA')
#'#plots the histogram of all atoms of ALA
#'#plt<-chemical_shift_hist('GLY',type='density')
#'#plots the density plot
#'@importFrom plotly %>%
#'@seealso \code{\link{fetch_res_chemical_shifts}},\code{\link{filter_residue}} and \code{\link{chem_shift_corr}} and \code{\link{atom_chem_shift_corr}}
chemical_shift_hist_res<-function(res='*',type='count',cutoff=8,interactive=TRUE){
    cs_dat2<-fetch_res_chemical_shifts(res)
  with(cs_dat2,{
    cs_dat2$Mean=NA
    cs_dat2$SD=NA

    if (all(is.na(cs_dat2))){
      return(NA)
    }else{
      for (atom in unique(cs_dat2$Atom_ID)){
        ss<-stats::sd(subset(cs_dat2,Atom_ID==atom)$Val)
        m<-mean(subset(cs_dat2,Atom_ID==atom)$Val)
        cs_dat2$Mean[cs_dat2$Atom_ID == atom] = m
        cs_dat2$SD[cs_dat2$Atom_ID == atom] = ss
        min_val<-m-(ss*cutoff)
        max_val<-m+(ss*cutoff)

        if (exists('cs_dat')){
          cs_dat<-rbind(cs_dat,subset(cs_dat2, Atom_ID == atom & Val>=min_val & Val<=max_val))
        }else{
          cs_dat<-subset(cs_dat2, Atom_ID == atom & Val>=min_val & Val<=max_val)
        }
      }
      cs_dat$Info = paste0("Mean=",round(cs_dat$Mean,2),";","SD=",round(cs_dat$SD,2))
      if (type=='count'){
        if (interactive){
          plt<-plotly::plot_ly(x = cs_dat$Val,type = "histogram", color = cs_dat$Atom_ID,alpha = 0.7, text = cs_dat$Info) %>% plotly::layout(barmode = "overlay")
        }
        else{
         plt<-ggplot2::ggplot(cs_dat)+
           ggplot2::geom_histogram(ggplot2::aes(x=Val,color=Atom_ID,fill=Atom_ID),binwidth=0.1,position = 'identity',alpha=0.7)+
           ggplot2::xlab("Chemical shift")+
           ggplot2::ylab("Count")+
           ggplot2::labs(color="",fill="")
        }
      } else{
        if (interactive) {
          plt<-plotly::plot_ly(x = cs_dat$Val,type = "histogram",color = cs_dat$Atom_ID,histnorm = "probability",alpha = 0.7, text = cs_dat$Info) %>% plotly::layout(barmode = "overlay")
        }
        else{
         plt<-ggplot2::ggplot(cs_dat)+
           ggplot2::geom_density(ggplot2::aes(x=Val,color=Atom_ID,fill=Atom_ID),alpha=0.7,trim=T)+
           ggplot2::xlab("Chemical shift")+
           ggplot2::ylab("Density")+
           ggplot2::labs(color="",fill="")
        }
      }

    }
    return(plt)
  })
}









#'Plots chemical shift distribution
#'
#'Plots the histogram (or) density of  chemical shift distribution of a given atom from amino acid (or) nucleic acid from BMRB database.  Optionally particular atom can be specified in the parameter
#'@param res residue name in NMR-STAR atom nomenclature ; Example: ALA,GLY ; default '*' (includes everything)
#'@param atm atom name in NMR-STAR nomenclautre ; Example :CA,HB2 default '*' (includes all atoms)
#'@param type count ; other than count will assume density plot
#'@param bw binwith for histogram; default value 0.1ppm
#'@param cutoff values not with in the cutoff time standard deviation from both sides ofthe mean will be excluded from the plot;default value 8
#'@return R plot object
#'@export chemical_shift_hist
#'@examples
#'#plt<-chemical_shift_hist('ALA')
#'#plots the histogram of all atoms of ALA
#'#plt<-chemical_shift_hist("*","CB*")
#'#plots  CB chemical shift distribution of standard amino acids
#'#plt<-chemical_shift_hist('GLY',type='density')
#'#plots the density plot
#'@importFrom plotly %>%
#'@seealso \code{\link{fetch_res_chemical_shifts}},\code{\link{filter_residue}} and \code{\link{chem_shift_corr}} and \code{\link{atom_chem_shift_corr}}
chemical_shift_hist<-function(res='*',atm='*',type='count',bw=0.1,cutoff=8){
  if (res=="*"){
    cs_dat2<-filter_residue(fetch_atom_chemical_shifts(atm))
  }else{
  cs_dat2<-filter_residue(fetch_res_chemical_shifts(res,atm))
  }
  with(cs_dat2,{
    cs_dat2$Mean=NA
    cs_dat2$SD=NA
    cs_dat2$bins=NA
    cs_dat2$UID=paste(cs_dat2$Comp_ID,"-",cs_dat2$Atom_ID)
  if (all(is.na(cs_dat2))){
    return(NA)
  }else{
    if (atm == "*" & res == "*"){
      for (atom in unique(cs_dat2$UID)){
        ss<-stats::sd(subset(cs_dat2,UID==atom)$Val)
        m<-mean(subset(cs_dat2,UID==atom)$Val)
        cs_dat2$Mean[cs_dat2$UID == atom] = m
        cs_dat2$SD[cs_dat2$UID == atom] = ss
        min_val<-m-(ss*cutoff)
        max_val<-m+(ss*cutoff)
        cs_dat2$bins[cs_dat2$UID == atom]=(max_val-min_val)/bw
        if (exists('cs_dat')){
          cs_dat<-rbind(cs_dat,subset(cs_dat2, UID == atom & Val>=min_val & Val<=max_val))
        }else{
          cs_dat<-subset(cs_dat2, UID == atom & Val>=min_val & Val<=max_val)
        }
      }
    }
    else if (atm == "*"){
    for (atom in unique(cs_dat2$Atom_ID)){
      ss<-stats::sd(subset(cs_dat2,Atom_ID==atom)$Val)
      m<-mean(subset(cs_dat2,Atom_ID==atom)$Val)
      cs_dat2$Mean[cs_dat2$Atom_ID == atom] = m
      cs_dat2$SD[cs_dat2$Atom_ID == atom] = ss
      min_val<-m-(ss*cutoff)
      max_val<-m+(ss*cutoff)
      cs_dat2$bins[cs_dat2$Atom_ID == atom]=(max_val-min_val)/bw
      if (exists('cs_dat')){
        cs_dat<-rbind(cs_dat,subset(cs_dat2, Atom_ID == atom & Val>=min_val & Val<=max_val))
      }else{
        cs_dat<-subset(cs_dat2, Atom_ID == atom & Val>=min_val & Val<=max_val)
      }
    }
    }else if (res == "*"){

        for (atom in unique(cs_dat2$Comp_ID)){
          ss<-stats::sd(subset(cs_dat2,Comp_ID==atom)$Val)
          m<-mean(subset(cs_dat2,Comp_ID==atom)$Val)
          cs_dat2$Mean[cs_dat2$Comp_ID == atom] = m
          cs_dat2$SD[cs_dat2$Comp_ID == atom] = ss
          min_val<-m-(ss*cutoff)
          max_val<-m+(ss*cutoff)
          cs_dat2$bins[cs_dat2$Comp_ID == atom]=(max_val-min_val)/bw
          if (exists('cs_dat')){
            cs_dat<-rbind(cs_dat,subset(cs_dat2, Comp_ID == atom & Val>=min_val & Val<=max_val))
          }else{
            cs_dat<-subset(cs_dat2, Comp_ID == atom & Val>=min_val & Val<=max_val)
          }
        }


    }else{
      ss<-stats::sd(cs_dat2$Val)
      m<-mean(cs_dat2$Val)
      cs_dat2$Mean= m
      cs_dat2$SD= ss
      min_val<-m-(ss*cutoff)
      max_val<-m+(ss*cutoff)
      cs_dat2$bins=(max_val-min_val)/bw
      cs_dat<-subset(cs_dat2, Val>=min_val & Val<=max_val)

    }
    cs_dat$Info = paste0("Mean=",round(cs_dat$Mean,2),";","SD=",round(cs_dat$SD,2))
    if (type=='count'){
      if (res=="*" & atm == "*"){
        plt<-suppressWarnings( plotly::plot_ly(x = cs_dat$Val,type = "histogram", color = cs_dat$UID,alpha = 0.6, nbinsx = cs_dat$bins,text = cs_dat$Info, visible = F) %>% plotly::layout(barmode = "overlay",xaxis = list(title="Chemical shift",zeroline = F),yaxis = list(title="Count",zeroline = F)))
      }else if (atm=="*") {
      plt<-suppressWarnings(plotly::plot_ly(x = cs_dat$Val,type = "histogram", color = cs_dat$Atom_ID,alpha = 0.6, nbinsx = cs_dat$bins,text = cs_dat$Info) %>% plotly::layout(barmode = "overlay",xaxis = list(title="Chemical shift",zeroline = F),yaxis = list(title="Count",zeroline = F)))
      }else if (res=="*") {
        plt<-suppressWarnings(plotly::plot_ly(x = cs_dat$Val,type = "histogram", color = cs_dat$Comp_ID,alpha = 0.6, nbinsx = cs_dat$bins,text = cs_dat$Info) %>% plotly::layout(barmode = "overlay",xaxis = list(title="Chemical shift",zeroline = F),yaxis = list(title="Count",zeroline = F)))
      }else{
        plt<-suppressWarnings(plotly::plot_ly(x = cs_dat$Val,type = "histogram",alpha = 0.6, nbinsx = cs_dat$bins,text = cs_dat$Info) %>% plotly::layout(barmode = "overlay",xaxis = list(title="Chemical shift",zeroline = F),yaxis = list(title="Count",zeroline = F)))
      }
      # plt<-ggplot2::ggplot(cs_dat)+
      #   ggplot2::geom_histogram(ggplot2::aes(x=Val,color=Atom_ID,fill=Atom_ID),binwidth=bw,position = 'identity',alpha=0.7)+
      #   ggplot2::xlab("Chemical shift")+
      #   ggplot2::ylab("Count")+
      #   ggplot2::labs(color="",fill="")
    } else{
      if (res=="*" & atm == "*"){
        plt<-suppressWarnings(plotly::plot_ly(x = cs_dat$Val,type = "histogram",color = cs_dat$UID,histnorm = "probability",alpha = 0.6,visible = F) %>% plotly::layout(barmode = "overlay",xaxis = list(title="Chemical shift",zeroline = F),yaxis = list(title="Probability density",zeroline = F)))
      }else if (atm == "*") {
      plt<-suppressWarnings(plotly::plot_ly(x = cs_dat$Val,type = "histogram",color = cs_dat$Atom_ID,histnorm = "probability",alpha = 0.6) %>% plotly::layout(barmode = "overlay",xaxis = list(title="Chemical shift",zeroline = F),yaxis = list(title="Probability density",zeroline = F)))
      }else if (res == "*"){
        plt<-suppressWarnings(plotly::plot_ly(x = cs_dat$Val,type = "histogram",color = cs_dat$Comp_ID,histnorm = "probability",alpha = 0.6) %>% plotly::layout(barmode = "overlay",xaxis = list(title="Chemical shift",zeroline = F),yaxis = list(title="Probability density",zeroline = F)))
      }else{
        plt<-suppressWarnings(plotly::plot_ly(x = cs_dat$Val,type = "histogram",histnorm = "probability",alpha = 0.6) %>% plotly::layout(barmode = "overlay",xaxis = list(title="Chemical shift",zeroline = F),yaxis = list(title="Probability density",zeroline = F)))
      }
      # plt<-ggplot2::ggplot(cs_dat)+
      #   ggplot2::geom_density(ggplot2::aes(x=Val,color=Atom_ID,fill=Atom_ID),alpha=0.7,trim=T)+
      #   ggplot2::xlab("Chemical shift")+
      #   ggplot2::ylab("Density")+
      #   ggplot2::labs(color="",fill="")
    }
    plt2<-plt
  }
  return(plt2)
  })
}

#'Plots chemical shift distribution for a list of atoms
#'
#'Plots the histogram (or) density of  chemical shift distribution of a given list of atoms. Atoms from different residues cam be specified as "residue-atom". Exammple "ALA-CA","GLN-HE21","GLN-HE*"
#'@param atm list Example: c("ALA-CA","GLY-CA")
#'@param type count ; other than count will assume density plot
#'@param bw binwith for histogram; default value 0.1ppm
#'@param cutoff values not with in the cutoff time standard deviation from both sides ofthe mean will be excluded from the plot;default value 8
#'@param interactive TRUE/FALSE default TRUE
#'@return R plot object
#'@export chemical_shift_hists
#'@examples
#'#plt<-chemical_shift_hists(c('ALA-C*'))
#'#plots the histogram of all atoms of ALA
#'#plt<-chemical_shift_hists(c("GLY-H*","ALA-HA"),type='density')
#'#plots the density plot
#'@seealso \code{\link{fetch_res_chemical_shifts}},\code{\link{filter_residue}} and \code{\link{chem_shift_corr}} and \code{\link{atom_chem_shift_corr}}
chemical_shift_hists<-function(atm=NA,type='count',bw=0.1,cutoff=8,interactive=TRUE){
  for (atom in atm){
    res<-strsplit(atom,"-")[[1]][1]
    at<-strsplit(atom,"-")[[1]][2]
    if (exists('cs_dat2')){
      cs_dat2<-rbind(cs_dat2,fetch_res_chemical_shifts(res,at))
    }else{
      cs_dat2<-fetch_res_chemical_shifts(res,at)
    }
  }
  with(cs_dat2,{
  #cs_dat2<-fetch_res_chemical_shifts(res,atm)

  if (all(is.na(cs_dat2))){
    return(NA)
  }else{
    for (atom in unique(cs_dat2$Atom_ID)){
      ss<-stats::sd(subset(cs_dat2,Atom_ID==atom)$Val)
      m<-mean(subset(cs_dat2,Atom_ID==atom)$Val)

      min_val<-m-(ss*cutoff)
      max_val<-m+(ss*cutoff)
      if (exists('cs_dat')){
        cs_dat<-rbind(cs_dat,subset(cs_dat2, Val>=min_val & Val<=max_val))
      }else{
        cs_dat<-subset(cs_dat2, Val>=min_val & Val<=max_val)
      }
    }
    if (type=='count'){
      plt<-ggplot2::ggplot(cs_dat)+
        ggplot2::geom_histogram(ggplot2::aes(x=Val,color=Atom_ID,fill=Comp_ID),binwidth=bw,position = 'identity',alpha=0.7)+
        ggplot2::xlab("Chemical shift")+
        ggplot2::ylab("Count")+
        ggplot2::labs(color="",fill="")
    } else{
      plt<-ggplot2::ggplot(cs_dat)+
        ggplot2::geom_density(ggplot2::aes(x=Val,color=Atom_ID,fill=Comp_ID),alpha=0.7,trim=T)+
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
  })
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
#'@importFrom plotly %>%
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
                plt<- plotly::plot_ly(x=hsqc_data$H,y=hsqc_data$N, color = hsqc_data$Entry_ID, text = hsqc_data$Info, type = "scatter", mode = "markers") %>%
                  plotly::layout(xaxis = list(autorange = "reversed",title="H",zeroline = F),yaxis = list(autorange = "reversed",title="N",zeroline = F) )
                #plt<-suppressWarnings(ggplot2::ggplot(hsqc_data)+
                # ggplot2::geom_point(ggplot2::aes(x=H,y=N,color=Entry_ID,label=Info))+ggplot2::labs(color=""))
                # ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()
                }
              else{
                plt<-suppressWarnings(ggplot2::ggplot(hsqc_data)+
                                        ggplot2::geom_point(ggplot2::aes(x=H,y=N,color=Entry_ID,label=Info))+
                                        ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+ggplot2::labs(color=""))
                }
              }
            else{
              if (interactive){
                plt<-plotly::plot_ly(x=hsqc_data$H,y=hsqc_data$N, color = hsqc_data$Comp_ID_H, text = hsqc_data$Info, type = "scatter", mode = "markers") %>%
                  plotly::layout(xaxis = list(autorange = "reversed",title="H",zeroline = F),yaxis = list(autorange = "reversed",title="N",zeroline = F))
              #plt<-suppressWarnings(ggplot2::ggplot(hsqc_data)+
              # ggplot2::geom_point(ggplot2::aes(x=H,y=N,color=Comp_ID_H,label=Info))+ggplot2::labs(color=""))
              # ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()
                }
              else{
                plt<-suppressWarnings(ggplot2::ggplot(hsqc_data)+
                                        ggplot2::geom_point(ggplot2::aes(x=H,y=N,color=Comp_ID_H,label=Info))+
                                        ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+ggplot2::labs(color=""))
                }
            }
            }
          else{
            if (interactive){
              plt<- plotly::plot_ly(x=hsqc_data$H,y=hsqc_data$N,color = hsqc_data$Entry_ID,text = hsqc_data$Info, type = "scatter", mode = "markers") %>%
                plotly::add_lines(split = hsqc_data$Comp_index_ID,color="x" ,showlegend=F) %>%
                plotly::layout(xaxis = list(autorange = "reversed",title="H",zeroline = F),yaxis = list(autorange = "reversed",title="N",zeroline = F) )
                #plt<-suppressWarnings(ggplot2::ggplot(hsqc_data)+
                #ggplot2::geom_line(ggplot2::aes(x=H,y=N,group=Comp_index_ID,label=Info))+
                #ggplot2::geom_point(ggplot2::aes(x=H,y=N,color=Entry_ID,label=Info))+ggplot2::labs(color=""))
                # ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()
              }
            else{
              plt<-suppressWarnings(ggplot2::ggplot(hsqc_data)+
                                      ggplot2::geom_line(ggplot2::aes(x=H,y=N,group=Comp_index_ID,label=Info))+
                                      ggplot2::geom_point(ggplot2::aes(x=H,y=N,color=Entry_ID,label=Info))+
                                      ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+ggplot2::labs(color=""))
            }
            }

        }
        plt2<-plt
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
    return(NA)}
  else{
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
              plt<- plotly::plot_ly(x=hsqc_data$H,y=hsqc_data$C, color = hsqc_data$Entry_ID, text = hsqc_data$Info, type = "scatter", mode = "markers") %>%
                plotly::layout(xaxis = list(autorange = "reversed",title="H",zeroline = F),yaxis = list(autorange = "reversed",title="C",zeroline = F) )
              #plt<-suppressWarnings(ggplot2::ggplot(hsqc_data)+
               #                       ggplot2::geom_point(ggplot2::aes(x=H,y=C,color=Entry_ID,label=Info))+ggplot2::labs(color=""))
                                      #ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()
              }
            else{
              plt<-suppressWarnings(ggplot2::ggplot(hsqc_data)+
                                      ggplot2::geom_point(ggplot2::aes(x=H,y=C,color=Entry_ID,label=Info))+
                                      ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+ggplot2::labs(color=""))
              }
            }
          else{
            if (interactive){
              plt<-plotly::plot_ly(x=hsqc_data$H,y=hsqc_data$C, color = hsqc_data$Comp_ID_H, text = hsqc_data$Info, type = "scatter", mode = "markers") %>%
                plotly::layout(xaxis = list(autorange = "reversed",title="H",zeroline = F),yaxis = list(autorange = "reversed",title="C",zeroline = F))
              #plt<-suppressWarnings(ggplot2::ggplot(hsqc_data)+
                                      #ggplot2::geom_point(ggplot2::aes(x=H,y=C,color=Comp_ID_H,label=Info))+ggplot2::labs(color=""))
                                      #ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()
              }
            else{
              plt<-suppressWarnings(ggplot2::ggplot(hsqc_data)+
                                      ggplot2::geom_point(ggplot2::aes(x=H,y=C,color=Comp_ID_H,label=Info))+
                                      ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+ggplot2::labs(color=""))
              }
            }
          }
        else{
          if (interactive){
            plt<- plotly::plot_ly(x=hsqc_data$H,y=hsqc_data$C,color = hsqc_data$Entry_ID,text = hsqc_data$Info, type = "scatter", mode = "markers") %>%
              plotly::add_lines(split = hsqc_data$Comp_index_ID,color="x" ,showlegend=F) %>%
              plotly::layout(xaxis = list(autorange = "reversed",title="H",zeroline = F),yaxis = list(autorange = "reversed",title="C",zeroline = F) )
            #plt<-suppressWarnings(ggplot2::ggplot(hsqc_data)+
                                    #ggplot2::geom_line(ggplot2::aes(x=H,y=C,group=Comp_index_ID,label=Info))+
                                    #ggplot2::geom_point(ggplot2::aes(x=H,y=C,color=Entry_ID,label=Info))+ggplot2::labs(color=""))\
            # ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()
            }
          else{
            plt<-suppressWarnings(ggplot2::ggplot(hsqc_data)+
                                    ggplot2::geom_line(ggplot2::aes(x=H,y=C,group=Comp_index_ID,label=Info))+
                                    ggplot2::geom_point(ggplot2::aes(x=H,y=C,color=Entry_ID,label=Info))+
                                    ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+ggplot2::labs(color=""))
          }
        }
        plt2<-plt
      }
      return(plt2)})}
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
          plt<- plotly::plot_ly(x=tocsy_data$Val.x,y=tocsy_data$Val.y, color = tocsy_data$Entry_ID, text = tocsy_data$Info, type = "scatter", mode = "markers") %>%
            plotly::layout(xaxis = list(autorange = "reversed",title="H",zeroline = F),yaxis = list(autorange = "reversed",title="H",zeroline = F) )
          #plt<-suppressWarnings(ggplot2::ggplot(tocsy_data)+
           #                       ggplot2::geom_point(ggplot2::aes(x=Val.x,y=Val.y,color=Entry_ID,label=Info))+
                                  #ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+
            #                      ggplot2::xlab("H")+ggplot2::ylab("H")+ggplot2::labs(color=""))
        }else{
          plt<-suppressWarnings(ggplot2::ggplot(tocsy_data)+
                                  ggplot2::geom_point(ggplot2::aes(x=Val.x,y=Val.y,color=Entry_ID,label=Info))+
                                  ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+
                                  ggplot2::xlab("H")+ggplot2::ylab("H")+ggplot2::labs(color=""))
        }
      }else{
        if (interactive){
          plt<-plotly::plot_ly(x=tocsy_data$Val.x,y=tocsy_data$Val.y, color = tocsy_data$Comp_ID_H, text = tocsy_data$Info, type = "scatter", mode = "markers") %>%
            plotly::layout(xaxis = list(autorange = "reversed",title="H",zeroline = F),yaxis = list(autorange = "reversed",title="H",zeroline = F))
         # plt<-suppressWarnings(ggplot2::ggplot(tocsy_data)+
          #                        ggplot2::geom_point(ggplot2::aes(x=Val.x,y=Val.y,color=Comp_ID.x,label=Info))+
                                  #ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+
           #                       ggplot2::xlab("H")+ggplot2::ylab("H")+ggplot2::labs(color=""))
        }else{
          plt<-suppressWarnings(ggplot2::ggplot(tocsy_data)+
                                  ggplot2::geom_point(ggplot2::aes(x=Val.x,y=Val.y,color=Comp_ID.x,label=Info))+
                                  ggplot2::scale_y_reverse()+ggplot2::scale_x_reverse()+ggplot2::labs(color="")+
                                  ggplot2::xlab("H")+ggplot2::ylab("H"))
        }
      }
     plt2<-plt
    }
    return(plt2) })})

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
#'@export chem_shift_corr
#'@examples
#'#plt<-chem_shift_corr('HE21','HE22')
#'#plots the chemical shift distribution between HE21 and HE22
#'@seealso \code{\link{fetch_atom_chemical_shifts}} and \code{\link{atom_chem_shift_corr}}
chem_shift_corr<-function(atom1,atom2,res=NA,type="c",interactive=TRUE){
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
  xl1=mean(cs$Val.x)-5*stats::sd(cs$Val.x)
  xl2=mean(cs$Val.x)+5*stats::sd(cs$Val.x)
  yl1=mean(cs$Val.y)-5*stats::sd(cs$Val.y)
  yl2=mean(cs$Val.y)+5*stats::sd(cs$Val.y)
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




#'Chemical shift correlation between given pair of atoms in a given amino acid (or) nucleic acid
#'
#'Plots the correlated chemical shift distribution of given pair of atoms in a single residue  from BMRB database.
#'By default it will generate interactive graphics using plotly library
#'@param atom1 atom name in NMR-STAR nomenclature like CA,CB2
#'@param atom2 atom name in NMR-STAR nomenclature like HA,HB2
#'@param res residue name in NMR-STAR nomenclature like ALA
#'@return plot object
#'@export atom_chem_shift_corr
#'@examples
#'#plt<-atom_chem_shift_corr('HE21','HE22','GLN')
#'#plots the chemical shift distribution between HE21 and HE22
#'@seealso \code{\link{fetch_res_chemical_shifts}} and \code{\link{chem_shift_corr}}
atom_chem_shift_corr<-function(atom1,atom2,res=NA){
  if (is.na(res)){
    stop("Residue must be specified")
  }
  else{
    css2<-fetch_res_chemical_shifts(res)
  with(css2,{if (grepl("\\*",atom1)){
    cs1<-subset(css2,grepl(paste0("^",sub("\\*","",atom1)),Atom_ID))
  }else{
      cs1<-subset(css2,Atom_ID == atom1)
  }
  if (grepl("\\*",atom2)){
    cs2<-subset(css2,grepl(paste0("^",sub("\\*","",atom2)),Atom_ID))
  }else{
  cs2<-subset(css2,Atom_ID == atom2)
  }
  if (all(is.na(cs1)) | all(is.na(cs2))){
    stop('One of the input atom has no data')
    plt2<-NA
    return(plt2)}
  else{
    with(cs1,{
      with(cs2, {
        cs<-merge(cs1,cs2,by=c('Entry_ID','Entity_ID','Comp_index_ID','Assigned_chem_shift_list_ID'))[,c("Entry_ID","Comp_index_ID","Entity_ID","Assigned_chem_shift_list_ID","Comp_ID.x","Comp_ID.y","Atom_ID.x","Atom_ID.y","Val.x","Val.y")]
        names(cs)[names(cs)=="Comp_ID.x"]<-"Comp_ID_1"
        names(cs)[names(cs)=="Comp_ID.y"]<-"Comp_ID_2"
        names(cs)[names(cs)=="Atom_ID.x"]<-"Atom_ID_1"
        names(cs)[names(cs)=="Atom_ID.y"]<-"Atom_ID_2"
        xl1=mean(cs$Val.x)-5*stats::sd(cs$Val.x)
        xl2=mean(cs$Val.x)+5*stats::sd(cs$Val.x)
        yl1=mean(cs$Val.y)-5*stats::sd(cs$Val.y)
        yl2=mean(cs$Val.y)+5*stats::sd(cs$Val.y)
        cs<-subset(cs,(Val.x>xl1 & Val.x<xl2 & Val.y>yl1 & Val.y<yl2))
        #names(cs)[names(cs)=="Val.x"]<-atom1
        #names(cs)[names(cs)=="Val.y"]<-atom2
        cs$Info=paste(cs$Comp_index_ID,cs$Entity_ID,cs$Atom_ID_1,cs$Atom_ID_2,cs$Assigned_chem_shift_list_ID,sep=",")
            plt<- plotly::subplot(
              plotly::plot_ly(x = cs$Val.x,type = "histogram" , text = cs$Info, alpha = 0.7,showlegend = F),
              suppressWarnings( plotly::plotly_empty(type = "scatter", mode = "markders")),
              plotly::plot_ly(x = cs$Val.x, y = cs$Val.y, type = "histogram2dcontour", showscale = FALSE, contours = list(coloring = "lines")) %>%
                plotly::layout(xaxis = list(title=atom1,zeroline = F),yaxis = list(title=atom2,zeroline = F)),
              plotly::plot_ly(y = cs$Val.y, type = "histogram", text = cs$Info, alpha = 0.7,showlegend = F),
              nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), margin = 0,
              shareX = TRUE, shareY = TRUE, titleX = TRUE, titleY = TRUE
            )
        return(plt)
      })})
  }})
  }
}








#'Filter for standard 20 amino acids
#'
#'Filters out non standard amino acids using Comp_ID. The data frame should contain three letter anio acid code in COMP_ID column.
#'@param df data frame with amino acid information in Comp_ID column
#'@return R data frame that contains information from only standard 20 amino acids.
#'@export filter_residue
#'@examples
#'#df<-filter_residue(fetch_atom_chemical_shifts("CG2"))
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



