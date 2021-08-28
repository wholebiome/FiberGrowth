#!/usr/bin/env Rscript

suppressWarnings(library(docopt))
suppressWarnings(library(knitr))
suppressWarnings(library(data.table))
suppressWarnings(library(ggplot2))
suppressWarnings(library(DT))
suppressWarnings(library(rmarkdown))
suppressWarnings(library(magrittr))
suppressWarnings(library(rhmmer))

doc <- '
FiberGrowth

Usage:
  FiberGrowth.R [-t <threads>] [--gffFeatureName=<gffFeatureName>] [--gffAttributeKey=<gffAttributeKey>] [--lib=<path_to_library>] --out=<output_folder> --gff=<gff_file> --proteins=<faa_file>
  FiberGrowth.R [-t <threads>]  [--lib=<path_to_library>] --out=<output_folder> --genome=<fasta_file>
  FiberGrowth.R (-h | --help)

Options:
  --gff=<gff_file>                          Gene locations in gff format.
  --proteins=<faa_file>                     Amino acid sequences in fasta format.
  --genome=<fasta_file>                     Genome in fasta format.
  --lib=<path_to_library>                   Path to library of PUL models 
                                            (default is PUL_models in installation directory).
  --out=<output_folder>                     Output directory (will be created).
  --gffFeatureName=<gffFeatureName>         Name of feature to use in gff file [default: CDS]
  --gffAttributeKey=<gffAttributeKey>       Name of attribute key to use in gff file [default: Name]
  -t --threads=<threads>                    Number of threads to run [default: 1].
  -h, --help                                Show this.
'

if(length(commandArgs(trailingOnly = TRUE)) == 0L) {
  docopt:::help(doc)
  quit()
}

arguments <- docopt(doc)


in_gff <- arguments$gff
in_proteins <- arguments$proteins
out_dir <- normalizePath(arguments$out)
gff_attribute_key <- arguments$gffAttributeKey
gff_feature_name <- arguments$gffFeatureName


##
# Functions


#' puldog: FiberGrowth main function (Polysaccharides Utilization Loci Detection On Genomes)
#'
#' @param hmmerResultFile
#' @param gffFile
#' @param pathwayFile a tab delimited file defining pathway members, order, orientation, and hmm_ID
#' @param gffAttributeKey The attribute name in the gff file that matches the header in the fasta file used for hmmer search. Defaults to locus_tag
#' @param gffFeatureName Feature name of rows on the gff file containing location of open reading frames. Defaults to CDS
#' @return data.table of clustered genes
#' @export
puldog <- function(hmmscan_tab,gffFile,pathwayTable,raw_prodigal_format=F, gffAttributeKey='locus_tag', gffFeatureName='CDS'){

  tabGff <- gff_tab
  pathwayTable[, RefPos := .I]
  tabPfamResults <- hmmscan_tab %>%
    .[,.(hmm_name,gene_id,sequence_score,description)] %>%
    #only keep entries that are in pathwayTable:
    .[ hmm_name %in% pathwayTable$hmm_name ]
  
  #Return -1 if hmmscan file does not contain any hits of pathwayTable:
  if(nrow(tabPfamResults) == 0){
    return(-1)
  }

  #Join HMMER results and gff table:
  # Keep only the top scores for each query name.
  tabPfamResults <- tabPfamResults[, .SD[which.max(sequence_score)], by="gene_id"]
  # tabPfamResults$gene_id %>% uniqueN == nrow(tabPfamResults)

  # Join HMMER and gff
  tabGffPfam <- tabGff[tabPfamResults, on = c(gene_id = "gene_id")]

  # label hits to pathway of interest:
  # tabGffPfamTarg <- tabGffPfam[pathwayTable, on = c(hmm_name = "hmm_name")]
  tabGffPfamTarg <- pathwayTable[tabGffPfam, on = c(hmm_name = "hmm_name")]

  # Standardize location to be "left" side of gene regardless of strand
  tabGffPfamTarg[, leftSide := min(start, end), by = "hashID"]
  setorder(tabGffPfamTarg, contig_ID, leftSide)

  # Define the number of entries per sequence
  tabGffPfamTarg[, nGenesPerSeq := .N, by = "contig_ID"]
  # tabGffPfamTarg %>% setorder(-nGenesPerSeq) %>% head


  #Return unclustered table if the maximum number of pathway genes per contig is 1:
  if(max(tabGffPfamTarg[,nGenesPerSeq]) == 1){
    return(tabGffPfamTarg)
  }

  # Compute location-based cluster based on distance
  # to other hits within same contig/sequence.
  # Every entry in Colocation is meaningful *within* its DNA sequence and cluster definition.
  tabGffPfamTarg[(nGenesPerSeq >= 2),
                 Colocation := dna_distance_clustering(loc = leftSide),
                 by = c("contig_ID", "PathwayName")]
  # Count number of co-located genes in each cluster
  tabGffPfamTarg[!is.na(Colocation),
                 nGenesPerColocation := .N,
                 by = c("contig_ID", "PathwayName", "Colocation")]

  # Count number of co-located genes in each cluster
  tabGffPfamTarg[!is.na(Colocation),
                 nUniqueGenesPerColocation := uniqueN(hmm_name),
                 by = c("contig_ID", "PathwayName", "Colocation")]

  # Exit if no co-located genes were found:
  if(all(is.na(tabGffPfamTarg$Colocation))){
    tabGffPfamTarg[,Synteny:=NA] # add Synteny column for downstream compatibility.
    return(tabGffPfamTarg)
  }

  #Calculate co-location:
  tabGffPfamTarg <-
    tabGffPfamTarg %>%
    # Only consider strand of co-location clusters
    .[!is.na(Colocation),
      cbind(
        hashID,
        check_orientation_match(
          qhmm_name = hmm_name,
          qStrand = strand,
          rStrand = RefStrand)
      ), by = c("contig_ID", "Colocation")] %>%
    # Drop variables from group-by return
    .[, c("contig_ID", "Colocation") := NULL] %>%
    # Join back with original table
    .[tabGffPfamTarg, on = c("hashID")]


  # Order full data table by position within each colocation cluster
  setorder(tabGffPfamTarg, contig_ID, Colocation, leftSide)
  tabGffPfamTarg[
    !is.na(Colocation),
    OrderMatchLeft := order_match(
      Side = "left",
      RefPos = RefPos,
      pathwayOrientation = pathwayOrientation),
    by = .(contig_ID, Colocation)]
  tabGffPfamTarg[
    !is.na(Colocation),
    OrderMatchRight := order_match(
      Side = "right",
      RefPos = RefPos,
      pathwayOrientation = pathwayOrientation),
    by = .(contig_ID, Colocation)]
  # Define a cluster-wise score for gene-order matching
  tabGffPfamTarg[
    !is.na(Colocation),
    OrderMatchScore := sum(OrderMatchLeft, OrderMatchRight, na.rm = TRUE),
    by = .(contig_ID, Colocation)]

  # Synteny (Order + Orientation):
  tabGffPfamTarg[, Synteny := StrandMatch & OrderMatchRight & OrderMatchLeft]

  tabGffPfamTarg[
    !is.na(Colocation) ,
    Growth := predict_growth_c(
      fiber=gsub('^[gnp]+_|_metabolism$','',PathwayName),
      geneName=GeneName,
      minSize=max(MinSize),
      pathwayTable=pathwayTable
    ),
    by = .(contig_ID, Colocation)]

  

  return(tabGffPfamTarg)
}


#' Parse gff file
#'
#' @param file Input file in gff format
#' @param calc_hash Calculate hash if T.
#' @param key The attribude in the gff file that matches to the header on the fasta file. Defaults to locus_tag
#' @param feature_name Feature name of rows containing location of open reading frames. Defaults to CDS
#' @param raw_prodigal_format Are gene_ids in gff in prodigal raw format format?
#' @return data.table of information from gff.
#'
read_gff <- function(file,calc_hash=F,key='locus_lag',feature_name='CDS',raw_prodigal_format=F){

  tabGff <- read.table(file,header = F,sep = '\t',quote = '',stringsAsFactors = F,na.strings = '',fill = T,comment.char = '#',
                       col.names = c('contig_ID','source','feature','start','end','score','strand','frame','attribute')) %>%
    data.table %>%
    na.omit()


  if(raw_prodigal_format){
    tabGff <- tabGff[, gene_id:=gsub(';.+','',attribute)] %>%
      .[,gene_id:=paste0(contig_ID,"_",gsub('ID=[0-9]+_','',
                                            gsub(';.+','',attribute)))] %>%
      .[,.(contig_ID,gene_id,feature,start,end,strand,attribute)]
  }else{

    tabGff <- tabGff[ feature == feature_name, ] %>%
      .[grepl(paste0(key,'='),attribute)] %>%
      .[,gene_id:=gsub(';.+','',
                       gsub(paste0('.*',key,'='),'',attribute))] %>%
      .[,.(contig_ID,gene_id,start,end,strand,attribute)]
  }

  if(calc_hash){
    tabGff$hashID <- sapply(tabGff$attribute, digest::digest)
  }

  return(tabGff[, !"attribute"])
}



#' Single-linkage clustering on DNA position of hits within each contig
#'
#' @param loc integer vector of gene location within a contig.
#' @param discriminationDistance The distance between genes to define a link.
#' @param method The agglomerative clustering method to use, passed to [hclust].
#'   Default is single-linkage clustering. Not recommended to change.
#'   This kind of cluster-definition is naturally greedy.
#' @return Returns a character vector of cluster-membership with length equal to `loc`/input.
#'
dna_distance_clustering = function(loc, discriminationDistance = 5000, method = "single"){
  x <-
    loc %>%
    dist() %>%
    hclust(method = method) %>%
    cutree(h = discriminationDistance) %>%
    paste0("cl", .)
  # Keep only colocation-clusters with two or more members
  keepCl <- x %>% table %>% is_greater_than(1) %>% which() %>% names
  x[!(x %in% keepCl)] <- NA_character_
  return(x)
}


#' Determine if Orietentation (strand) matches the reference pathway
#'
#' @param qhmm_name vector of pfam IDs that match entries in `pathwayTable`.
#' @param qStrand vector of query strand indicator as "+" or "-".
#' @param rStrand vector of cluster-reference strand indicator as "+" or "-".
#' @return data.table with logical vector of the same length as `qhmm_name` and `qStrand`,
#' indicating whether the orientation (strand) of each gene matches the reference,
#' and a second column indicating whether the whole cluster is in sense
#' or anti-sense relative to the reference.
check_orientation_match = function(qhmm_name, qStrand, rStrand){
  tabTest <-
    data.table(
      hmm_name = qhmm_name,
      Strand = qStrand,
      RefStrand = rStrand)
  tabTest[, sense := Strand == RefStrand]
  tabTest[, antisense := Strand != RefStrand]
  winningStrand <-
    c(
      antisense = tabTest$antisense %>% sum(na.rm = TRUE),
      sense = tabTest$sense %>% sum(na.rm = TRUE)) %>%
    which.max %>%
    names()
  return(
    data.table(
      StrandMatch = tabTest[[winningStrand]],
      pathwayOrientation = winningStrand
    )
  )
}


#' Determine if order of gene cluster matches reference cluster.
#'
#' @param RefPos position of genes in reference cluster.
#' @param pathwayOrientation indication if pathway is on sense or antisense strand.
#' @param Side side of neighboring cluster genes to check.
#' @return data.table
order_match = function(RefPos, pathwayOrientation, Side = "left"){
  lagStem <- RefPos %>% diff(lag = 1)
  orderMatch <-
    switch(
      EXPR = Side,
      right = c(lagStem, 0),
      left = c(0, lagStem)
    ) %>%
    multiply_by(
      e2 = (
        switch(pathwayOrientation[1],
               antisense = -1,
               sense = 1)
      ))
  orderMatch[orderMatch!=1L] <- 0
  return(orderMatch)
}
predict_growth_c <- function(fiber,geneName,minSize=0,pathwayTable=NULL){
  
  growth <- F
  cluster <- geneName
  if (fiber == "inulin"){
    growth <- any(grepl("GH", cluster)) & any(grepl("fructokinase", cluster)) & any(grepl("transporter", cluster)) & length(geneName) >= minSize #3
  }  else
    if (fiber == "arabinoxylan"){ 
      growth <- any(grepl("GH43_1", cluster)) & any(grepl("GH10", cluster)) & any(grepl("carbohydrate_esterase", cluster)) & length(geneName) >= minSize #5
    } else
      if (fiber == "heparin") {
        growth <- any(grepl("sulfatase", cluster)) & any(grepl("GH88", cluster)) & any(grepl("PL15", cluster)) & (any(grepl('transporter',cluster)) | any(grepl('TonB',cluster)) )  & length(geneName) >= minSize #7
      } else
        if (fiber == "laminarin") {
          growth <- any(grepl("^GH16$", cluster)) & any(grepl("^GH3", cluster)) & any(grepl("^TonB", cluster)) & length(geneName) >= minSize #3
        } else
          if (fiber == "levan"){
            growth <- any(grepl("^GH32$", cluster)) & any(grepl("GH32Btheta_like", cluster)) & any(grepl("glucose_transporter", cluster)) & any(grepl("fructokinase", cluster)) & length(geneName) >= minSize #5
          } else
            if (fiber == "starch"){
              growth <- any(grepl("TonB", cluster)) & any(grepl("SusD_RagB", cluster)) & any(grepl("GH13", cluster)) & length(geneName) >= minSize #5
            }else
              if (fiber == "mucin"){
                growth <- any(grepl("TonB", cluster)) & any(grepl("SusD", cluster)) & any(grepl("GH18", cluster)) & any(grepl("glucanase", cluster)) & length(geneName) >= minSize #5
              }else{
                #custom:
                if(sum(pathwayTable$Weight == 0) == 0){
                  growth <- (length(geneName) >= minSize)
                }else{
                  growth <- (length(geneName) >= minSize && all(pathwayTable[ Weight != 0, GeneName] %in% cluster))
                }
              }
  
  return(growth)
}
# modified version of the `read_tblout` function of the `rhmmer` package by Zebulun Arendsee (https://github.com/arendsee/rhmmer)
#'
#' @param file hmmscan tblout file.
#' @return data.table
read_tblout_quietly <- function(file){
  col_types <- readr::cols(domain_name = readr::col_character(), domain_accession = readr::col_character(), 
                           query_name = readr::col_character(), query_accession = readr::col_character(), 
                           sequence_evalue = readr::col_double(), sequence_score = readr::col_double(), 
                           sequence_bias = readr::col_double(), best_domain_evalue = readr::col_double(), 
                           best_domain_score = readr::col_double(), best_domain_bis = readr::col_double(), 
                           domain_number_exp = readr::col_double(), domain_number_reg = readr::col_integer(), 
                           domain_number_clu = readr::col_integer(), domain_number_ov = readr::col_integer(), 
                           domain_number_env = readr::col_integer(), domain_number_dom = readr::col_integer(), 
                           domain_number_rep = readr::col_integer(), domain_number_inc = readr::col_character(), 
                           description = readr::col_character())
  
  N <- length(col_types$cols)
  readr::read_lines(file, progress = F) %>% 
    sub(pattern = sprintf("(%s) *(.*)",paste0(rep("\\S+", N - 1), collapse = " +")), replacement = "\\1\t\\2", perl = TRUE) %>%
    paste0(collapse = "\n") %>% 
    readr::read_tsv(col_names = c("X","description"), comment = "#", na = "-",show_col_types = F, progress = F) %>%
    tidyr::separate(.data$X,head(names(col_types$cols), -1), sep = " +") %>% 
    readr::type_convert(col_types = col_types) %>% 
    data.table()
}
get_this_path <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
exit <- function() {
  invokeRestart("abort") 
} 

# Functions
##


# Check dependencies
dependencies <- Sys.which(c("prodigal","hmmscan","hmmpress")) %>%
  data.table(dependency=names(.),path=.)

if(any(dependencies[,path] == '')){
  cat('ERROR: Missing dependencies: ',paste(dependencies[path == '', dependency],sep='',collapse = ', '),'\n')
  exit()
}




# Create output directory:
dir.create(paste0(out_dir,'/hmmscan'),recursive = T,showWarnings = F)
dir.create(paste0(out_dir,'/pulscan'),showWarnings = F)
dir.create(paste0(out_dir,'/fibergrowth'),showWarnings = F)

raw_prodigal_format <- F
if(!is.null(arguments$genome)){
  dir.create(paste0(out_dir,'/prodigal'),showWarnings = F)

  in_gff <- paste0(out_dir,'/prodigal/',basename(arguments$genome),'_genomic.gff')
  in_proteins <- paste0(out_dir,'/prodigal/',basename(arguments$genome),'_prot.faa')
  system(paste0('prodigal -i "',arguments$genome,'" -o ', in_gff ,' -a ',in_proteins,' -f gff -g 11 -p single -m -q'))
  raw_prodigal_format <- T
  gff_attribute_key <- 'ID'
  gff_feature_name <- 'CDS'
}


if(is.null(arguments$lib)){
  pul_library_path <- paste0(dirname(get_this_path()),'/PUL_models')
}else{
  pul_library_path <- arguments$lib
}


# List all available fibers in library:
lib_fibers <- data.table(PUL_table=list.files(pul_library_path,pattern = 'pathway_table.tsv',recursive = T,full.names = T)) %>%
  .[,PUL_model:=basename(dirname(PUL_table))] %>%
  .[,PUL_hmm:=paste(dirname(PUL_table),'protein_families.hmm',sep='/')] %>%
  # .[grepl('^gp_|^gn_',PUL_model)] %>%
  .[,fiber_name:=gsub('^gp_|^gn_','',PUL_model)] %>%
  .[,gram:=apply(.[,.(fiber_name,PUL_model)],1,function(x){gsub(paste0('_',x[1]),'',x[2])})]

if(nrow(lib_fibers) == 0){
  cat(paste0('\n**No PUL models found in ',lib_fibers,'**\n'))
  knitr::knit_exit()
}


for(pul in lib_fibers[,PUL_model]){

  pul_hmm <- lib_fibers[PUL_model == pul,PUL_hmm]
  hmmscan_out <- paste0(out_dir,'/hmmscan/',lib_fibers[PUL_model == pul,PUL_model],'_hmmscan.tblout')

  ## check hmmpress files, (protein_families.hmm.h3m,...)
  if(!file.exists(paste0(pul_hmm,'.h3m'))){
    hmmpress_command <- paste0('hmmpress ',pul_hmm)
    cat(paste0('File not found: ',pul_hmm,'.h3m','\n >running: ',hmmpress_command,'\n'))
    system(hmmpress_command)
  }
    
  system(paste0('hmmscan --tblout ', hmmscan_out ,' --cut_ga --cpu ',arguments$threads,' "',
                pul_hmm,'" "',in_proteins,'">/dev/null 2>&1'))

  # Read hmmscan result:
  hmmscan_tab <- read_tblout_quietly(hmmscan_out) %>%
    data.table() %>%
    setnames(old = c('domain_name','query_name'),new=c('hmm_name','gene_id'))

  # Read gff file:
  gff_tab <- read_gff(file = in_gff,calc_hash = T, key=gff_attribute_key,feature_name = gff_feature_name, raw_prodigal_format=raw_prodigal_format)

  # Read PUL definition:
    pathwayTable <- fread(lib_fibers[PUL_model == pul,PUL_table]) %>%
    setnames(old = c('HmmID'),new = c('hmm_name'))


    #Stop if no genes are found:
  if(nrow(hmmscan_tab)==0){
    cat('\n**No PUL genes found**\n')
    knitr::knit_exit()
  }


  # Run PULDOG:
  pulscan_out <- paste0(out_dir,'/pulscan/',lib_fibers[PUL_model == pul,PUL_model],'_pulscan.tsv')
  pul_tab <- puldog(hmmscan_tab = hmmscan_tab,
                               gffFile = gff_tab,
                               pathwayTable = pathwayTable)

  pul_tab_sub <- pul_tab[!is.na(Colocation)] %>% 
    .[,fiber:=gsub('^[gnp]+_|_metabolism$','',PathwayName)] %>% 
    .[,pul_id:=paste(fiber,contig_ID,Colocation,sep='_')] %>% 
    .[ ,.(fiber,
          pul_id,
          gene_name=GeneName,
          hmm_name,
          position=Position,
          contig_id=contig_ID,
          gene_id,
          start,
          end,
          strand,
          description,
          n_genes=nGenesPerColocation,
          n_unique_genes=nUniqueGenesPerColocation,
          growth=Growth)]
  
  write.table(pul_tab_sub,pulscan_out,col.names=T,row.names=F,quote=F,sep='\t')


  if(nrow(pul_tab_sub) > 1 && pul_tab_sub[growth == T ,.N] > 0){

    fibergrowth_out <- paste0(out_dir,'/fibergrowth/',lib_fibers[PUL_model == pul,PUL_model],'_fibergrowth.tsv')
    pul_tab_sum <- pul_tab_sub[growth == T,
                                       list(genes=paste0('[',paste(gene_name,sep='-',collapse = ']-['),']'),
                                            start=min(start),
                                            end=max(end)),
                                       by=.( contig_id, pul_id,n_genes,growth)]
    write.table(pul_tab_sum,fibergrowth_out,col.names=T,row.names=F,quote=F,sep='\t')
  }
}

this_path <- dirname(get_this_path())
rmarkdown::render(paste0(this_path,"/FiberGrowthReport.Rmd"),
                  output_file=paste0(out_dir,'/FiberGrowthReport.html'),
                  params = list(input_dir=paste0(out_dir,'/pulscan')),
                  quiet = T)
