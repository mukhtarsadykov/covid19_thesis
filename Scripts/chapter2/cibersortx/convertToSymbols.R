# Convert ENSGID into Gene SYmbol
a1 <- read.table("~/Desktop/phd/Thesis/ch2/cts_11.txt", sep = "\t", header = T)

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Example list of Ensembl Gene IDs
ensg_ids <- a1$Geneid

# Query biomaRt to get the corresponding gene symbols
result <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                filters = 'ensembl_gene_id', 
                values = ensg_ids, 
                mart = ensembl)
names(result) <- c("Geneid", "Symbol")

m1 <- merge(a1, result, by = "Geneid")
m1 <- m1[, c(13, 2:12)]

m1 <- m1 %>% filter(Symbol != "" & !is.na(Symbol))
write.table(m1, "~/Desktop/phd/Thesis/ch2/cts_11_symbols.txt", sep ="\t", quote = FALSE, row.names = F)

