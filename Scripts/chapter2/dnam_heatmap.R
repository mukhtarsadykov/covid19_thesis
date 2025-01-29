### DNAm###
#10/04/2024#

# from DMRicher coverage 1, 11 samples
setwd("~/Desktop/phd/Thesis/ch2/dnam/")
dmrs <- read.table("dmricher/DMRs/DMR_individual_smoothed_methylation.txt", header = T, sep = "\t")
dmrs <- dmrs[,17:length(dmrs)]
names(dmrs) <- c("ICU_01",	"ICU_02",	"ICU_03",	"ICU_04",	"NON_01",	"NON_02",	"NON_03",	"NON_04",	"NON_05",	"NON_06",	"NON_07")

study <- read.csv("~/Desktop/phd/Thesis/ch2/meta_11.csv", row.names=1) # dim is 80 x 19

study <- study[1:11,]

study$AltName <- c("ICU_01", "ICU_02", "ICU_03", "ICU_04", "NON_01", "NON_02",
                   "NON_03", "NON_04", "NON_05", "NON_06", "NON_07")

save_pheatmap_pdf <- function(x, filename, width=9, height=9) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Reodered study according to heatmap clustering
#study <- read.csv("~/Desktop/rnaseq/last/raw/80s_new_metadata/meta_80s.csv", row.names = 1)

mycol <- colorRampPalette(c("blue","white","red"))(20)

#annotation info
annotation <- data.frame(study$AltName, study$Outcome,
                         study$Age, study$Hospital, study$ICU, study$Gender, study$Comorbidities, row.names = 1)
colnames(annotation) <- gsub("study.", "", colnames(annotation))
rownames(annotation) <- colnames(lcpm) # check out the row names of annotation
##colors
Age_col <- as.character(wes_palette("IsleofDogs2",  
                                    40, 
                                    type = "continuous"))
names(Age_col) <- unique(study$Age)
Hospital_col <- as.character(wes_palette("Chevalier1",  
                                         length(unique(study$Hospital)), 
                                         type = "continuous"))
names(Hospital_col) <- unique(study$Hospital)
libsize_col <- as.character(wes_palette("IsleofDogs2",  
                                        76, 
                                        type = "continuous"))
names(libsize_col) <- unique(study$lib.size)

ann_colors = list(
  ICU = c(Yes="#972D15", No ="#D8B70B"),
  Comorbidities = c(none="#E1BD6D", one="#9FA682", two="#274150", three="#F2300F"),
  #Nationality = c(Saudi="#F1BB7B", Yemeni="#FD6467", Other="#5B1A18"),
  Gender = c(Female = "#D783FF", Male = "#00DED8"),
  Outcome = c(Dead="#972D15", Alive ="#D8B70B"),
  Age = Age_col,
  Hospital = Hospital_col
  #Library = c(low="#972D15", high ="#D8B70B")
  #lib.size = libsize_col
)
annotation <- annotation[, c("ICU", "Outcome", "Age", "Hospital", "Gender", "Comorbidities")]

#write.csv(lcpm[i,], "~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/degs_67_bn.csv")
heatmap1 <- pheatmap(dmrs, cluster_cols = F, cluster_rows = T,
                     scale = "row", #kmeans_k = 20,
                     cutree_rows = 2, cutree_cols = 2,
                     color = colorRampPalette(c("blue", "white", "red"))(20),
                     show_rownames = 0, annotation_col = annotation, 
                     #cellwidth = 8, cellheight = 8, border_color = "black",
                     labels_col=study$AltName, annotation_colors = ann_colors,
                     treeheight_row=0
)

save_pheatmap_pdf(heatmap1, "~/Desktop/phd/Thesis/ch2/dnam/dmricher/Heatmap_annotation_11_ordered.pdf")
