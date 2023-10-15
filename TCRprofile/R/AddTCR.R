#' Add TCR to Seurat object
#' This function adds Cell Ranger TCR output information with transcriptome, cell hashing, 
#' and others if applicable
#' To output a comprehensive table to contain all information that originate from a single 
#' barcode
#' 
#' Three CSV files generated from Cell Ranger are required as input files for this function
#' @param clono clonotypes.csv from Cell Ranger 
#' @param consensus consensus_annotation.csv from Cell Ranger
#' @param filtered_barcodes filtered_contig_annotations.csv from Cell Ranger
#' @importFrom dplyr bind_rows summarize %>% count group_by
#' @importFrom reshape2 dcast
#' @importFrom stringr str_c
#' @export
#' @return a meta table contains all TCR info



#' Check how they solve duplicated barcodes
#' Using the samples and ID parameters, the 
#' function will add the strings as prefixes to prevent issues with repeated 
#' barcodes
#' They used samples and ID to label barcodes in order to distinguish repeated barcodes
#' among different samples/IDs.
#' However, in my previous cases, duplicated barcodes were observed within one single 
#' run/library

# suppressMessages(library(dplyr))
# suppressMessages(library(reshape2))
# suppressMessages(library(stringr))

AddTCR <- function(clono, 
                   consensus, 
                   filtered_barcodes) {
  #clono<-read.csv("/Users/aliu/Documents/Stuart_Pushpak/10X_11/CD3_1_out/vdj_out/clonotypes.csv", 
  #                header = TRUE, sep=",")
  #format the csv - selecting columns 
  clono<-clono[, c(1:3, 6,7)]
  #consensus<-read.csv("/Users/aliu/Documents/Stuart_Pushpak/10X_11/CD3_1_out/vdj_out/consensus_annotations.csv", 
  #               header = TRUE, sep = ",")
  #format the csv - deleting columns that not informative
  consensus<-consensus[,c(-3,-9,-10,-13:-22)]
  
  #Categorise clonotype_id and chain
  cons_freq<-dcast(count(consensus[,c(1,3)], clonotype_id,chain), clonotype_id ~ chain)
  cons_freq[is.na(cons_freq)]<-0
  
  #For paired clonotypes: 
  #TRA=1 and TRB=1
  #TRA=2 and TRB=1
  paired_clono<-na.omit(cons_freq[cons_freq$TRA==1&cons_freq$TRB==1 | 
                                    cons_freq$TRA==2 & cons_freq$TRB==1,])
  paired_clono<-cons[cons$clonotype_id%in%paired_clono$clonotype_id,]
  
  tra<-paired_clono[paired_clono$chain=="TRA",]
  trb<-paired_clono[paired_clono$chain=="TRB",]
  
  #consolidating all info relating to alpha chain
  tra<-tra%>%group_by(clonotype_id)%>%
    summarize(consensus_id=str_c(consensus_id,collapse = ";"),
              tra_v_gene=str_c(v_gene,collapse = ";" ),
              tra_d_gene=str_c(d_gene,collapse = ";" ),
              tra_j_gene=str_c(j_gene,collapse = ";" ),
              tra_c_gene=str_c(c_gene,collapse = ";" ),
              # tra_fwr1=str_c(fwr1, collapse = ";"),
              # tra_fwr1_nt=str_c(fwr1_nt, collapse = ";"),
              # tra_cdr1=str_c(cdr1, collapse = ";"),
              # tra_cdr1_nt=str_c(cdr1_nt, collapse = ";"),
              # tra_fwr2=str_c(fwr2, collapse = ";"),
              # tra_fwr2_nt=str_c(fwr2_nt, collapse = ";"),
              # tra_cdr2=str_c(cdr2, collapse = ";"),
              # tra_cdr2_nt=str_c(cdr2_nt, collapse = ";"),
              # tra_fwr3=str_c(fwr3, collapse = ";"),
              # tra_fwr3_nt=str_c(fwr3_nt, collapse = ";"),
              tra_cdr3=str_c(cdr3, collapse = ";"),
              tra_cdr3_nt=str_c(cdr3_nt, collapse = ";"))
  # tra_fwr4=str_c(fwr4, collapse = ";"),
  # tra_fwr4_nt=str_c(fwr4_nt, collapse = ";"))
  
  #consolidating all info relating to beta chain
  trb<-trb%>%group_by(clonotype_id)%>%
    summarize(consensus_id=str_c(consensus_id,collapse = ";"),
              trb_v_gene=str_c(v_gene,collapse = ";" ),
              trb_d_gene=str_c(d_gene,collapse = ";" ),
              trb_j_gene=str_c(j_gene,collapse = ";" ),
              trb_c_gene=str_c(c_gene,collapse = ";" ),
              # trb_fwr1=str_c(fwr1, collapse = ";"),
              # trb_fwr1_nt=str_c(fwr1_nt, collapse = ";"),
              # trb_cdr1=str_c(cdr1, collapse = ";"),
              # trb_cdr1_nt=str_c(cdr1_nt, collapse = ";"),
              # trb_fwr2=str_c(fwr2, collapse = ";"),
              # trb_fwr2_nt=str_c(fwr2_nt, collapse = ";"),
              # trb_cdr2=str_c(cdr2, collapse = ";"),
              # trb_cdr2_nt=str_c(cdr2_nt, collapse = ";"),
              # trb_fwr3=str_c(fwr3, collapse = ";"),
              # trb_fwr3_nt=str_c(fwr3_nt, collapse = ";"),
              trb_cdr3=str_c(cdr3, collapse = ";"),
              trb_cdr3_nt=str_c(cdr3_nt, collapse = ";"))
  # trb_fwr4=str_c(fwr4, collapse = ";"),
  # trb_fwr4_nt=str_c(fwr4_nt, collapse = ";"))
  
  #Sanity check
  #Double check same number of clonotypes 
  length(unique(paired_clono$clonotype_id))==length(unique(tra$clonotype_id))
  length(unique(paired_clono$clonotype_id))==length(unique(trb$clonotype_id))
  
  #Merge alpha and beta chains together 
  paired_clono<-merge(tra,trb, by= "clonotype_id")
  paired_clono$consensus_id<-paste0(paired_clono$consensus_id.x, sep=";"
                                    ,paired_clono$consensus_id.y)
  #paired_clono<-paired_clono[,c(-2,-21)]
  #change column range to exclude cols 2 and 9 
  #to keep only one column as consensus_id (excluding tra_consensus_id & trb_consensus_id)
  paired_clono<-paired_clono[,c(-2,-9)]
  
  #Combine with further info
  #Combine with Freq, inkt_evidence and mait_evidence 
  paired_clono<-merge(clono, paired_clono, by="clonotype_id")
  #paired_clono<-paired_clono[,c(1:5,42,6:41)]
  
  #Only one chain
  #TRB only clonotype:
  #TRA = 0, TRB=1
  trb_only<-na.omit(cons_freq[cons_freq$TRA==0&cons_freq$TRB==1,])
  trb_only_clono<-consensus[consensus$clonotype_id%in%trb_only$clonotype_id,]
  #Combine with Freq, inkt_evidence and mait_evidence
  trb_only_clono<-merge(clono, trb_only_clono, by="clonotype_id")
  #exclude "chain" column, as TRB for all clonotypes
  trb_only_clono<-trb_only_clono[,-7]
  #colnames(trb_only_clono)[7:24]<-paste0("trb_",colnames(trb_only_clono)[7:24])
  #column numbers set to 7-12
  colnames(trb_only_clono)[7:12]<-paste0("trb_",colnames(trb_only_clono)[7:12])
  
  #TRA only clonotype:
  #TRA=1, TRB=0
  tra_only<-na.omit(cons_freq[cons_freq$TRA==1&cons_freq$TRB==0,])
  tra_only_clono<-consensus[consensus$clonotype_id%in%tra_only$clonotype_id,]
  #Combine with Freq,inkt_evidence and mait_evidence
  tra_only_clono<-merge(clono, tra_only_clono, by="clonotype_id")
  tra_only_clono<-tra_only_clono[,-7]
  #head(tra_only_clono)
  #colnames(tra_only_clono)[7:24]<-paste0("tra_",colnames(tra_only_clono)[7:24])
  #column numbers set to 7 - 12
  colnames(tra_only_clono)[7:12]<-paste0("tra_",colnames(tra_only_clono)[7:12])
  
  #More than 2 chains
  #Multi TRA/TRB
  #TRA =1 & TRB>1 or TRA >=2 & TRB >1 or TRA >2 &TRB=1
  multi<-na.omit(cons_freq[(cons_freq$TRA==1&cons_freq$TRB>1) | (cons_freq$TRA>=2&cons_freq$TRB>1) |(cons_freq$TRA>2&cons_freq$TRB==1), ])
  multi_clono<-consensus[consensus$clonotype_id%in%multi$clonotype_id,]
  multi_tra<-multi_clono[multi_clono$chain=="TRA",]
  multi_trb<-multi_clono[multi_clono$chain=="TRB",]
  
  multi_tra<-multi_tra%>%group_by(clonotype_id)%>%
    summarize(tra_consensus_id=str_c(consensus_id,collapse = ";"),
              tra_v_gene=str_c(v_gene,collapse = ";" ),
              tra_d_gene=str_c(d_gene,collapse = ";" ),
              tra_j_gene=str_c(j_gene,collapse = ";" ),
              tra_c_gene=str_c(c_gene,collapse = ";" ),
              # tra_fwr1=str_c(fwr1, collapse = ";"),
              # tra_fwr1_nt=str_c(fwr1_nt, collapse = ";"),
              # tra_cdr1=str_c(cdr1, collapse = ";"),
              # tra_cdr1_nt=str_c(cdr1_nt, collapse = ";"),
              # tra_fwr2=str_c(fwr2, collapse = ";"),
              # tra_fwr2_nt=str_c(fwr2_nt, collapse = ";"),
              # tra_cdr2=str_c(cdr2, collapse = ";"),
              # tra_cdr2_nt=str_c(cdr2_nt, collapse = ";"),
              # tra_fwr3=str_c(fwr3, collapse = ";"),
              # tra_fwr3_nt=str_c(fwr3_nt, collapse = ";"),
              tra_cdr3=str_c(cdr3, collapse = ";"),
              tra_cdr3_nt=str_c(cdr3_nt, collapse = ";"))
  # tra_fwr4=str_c(fwr4, collapse = ";"),
  # tra_fwr4_nt=str_c(fwr4_nt, collapse = ";"))
  
  multi_trb<-multi_trb%>%group_by(clonotype_id)%>%
    summarize(trb_consensus_id=str_c(consensus_id,collapse = ";"),
              trb_v_gene=str_c(v_gene,collapse = ";" ),
              trb_d_gene=str_c(d_gene,collapse = ";" ),
              trb_j_gene=str_c(j_gene,collapse = ";" ),
              trb_c_gene=str_c(c_gene,collapse = ";" ),
              # trb_fwr1=str_c(fwr1, collapse = ";"),
              # trb_fwr1_nt=str_c(fwr1_nt, collapse = ";"),
              # trb_cdr1=str_c(cdr1, collapse = ";"),
              # trb_cdr1_nt=str_c(cdr1_nt, collapse = ";"),
              # trb_fwr2=str_c(fwr2, collapse = ";"),
              # trb_fwr2_nt=str_c(fwr2_nt, collapse = ";"),
              # trb_cdr2=str_c(cdr2, collapse = ";"),
              # trb_cdr2_nt=str_c(cdr2_nt, collapse = ";"),
              # trb_fwr3=str_c(fwr3, collapse = ";"),
              # trb_fwr3_nt=str_c(fwr3_nt, collapse = ";"),
              trb_cdr3=str_c(cdr3, collapse = ";"),
              trb_cdr3_nt=str_c(cdr3_nt, collapse = ";"))
  # trb_fwr4=str_c(fwr4, collapse = ";"),
  # trb_fwr4_nt=str_c(fwr4_nt, collapse = ";"))
  multi_clono<-merge(multi_tra,multi_trb, by="clonotype_id")
  
  #Combine with Freq, inkt_evidence and mait_evidence
  multi_clono<-merge(clono, multi_clono, by="clonotype_id")
  multi_clono$consensus_id<-paste0(multi_clono$tra_consensus_id, sep=";" ,multi_clono$trb_consensus_id)
  # multi_clono<-multi_clono[,c(-6,-25)]
  #change columns 6 & 13 to be exlcuded
  multi_clono<-multi_clono[,c(-6,-13)]
  # multi_clono<-multi_clono[,c(1:5,42,6:41)]
  #head(multi_clono)
  
  # Append four categories into one meta table
  # Add an extra column to specify paired or not
  paired_clono$pair<-"paired"
  multi_clono$pair<-"multi"
  tra_only_clono$pair<-"tra"
  trb_only_clono$pair<-"trb"
  clono_meta<-bind_rows(paired_clono, multi_clono, tra_only_clono, trb_only_clono)
  
  #Link with barcodes
  #barcodes<-read.csv("/Users/aliu/Documents/Stuart_Pushpak/10X_11/CD3_1_out/vdj_out/filtered_contig_annotations.csv", head=TRUE, sep=",")
  
  
  filtered_barcodes<-filtered_barcodes[!duplicated(filtered_barcodes$barcode),]#matches up with cellranger vdj-t estimated number of cells
  
  filtered_barcodes<-filtered_barcodes[,c("barcode", "raw_clonotype_id")]
  names(filtered_barcodes)[names(filtered_barcodes) == "raw_clonotype_id"] <- "clonotype_id"
  
  clono_meta<-merge(filtered_barcodes,clono_meta,by="clonotype_id")
  
  #diff in barcodes: one barcode missing clonotype_id
  diff_barcodes<-filtered_barcodes[!filtered_barcodes$clonotype_id%in%clono_meta$clonotype_id,]
  
  rownames(clono_meta)<-clono_meta$barcode
  clono_meta$barcode<-NULL
  
  return(clono_meta)
  
  #skipping adding to seurat object
  #fb_subset_vdj<-AddMetaData(object = fb_subset, metadata = clono_meta)
  #fb_subset_vdj
}