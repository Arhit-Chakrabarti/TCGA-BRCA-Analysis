# READ THE log2(x+1) GENE EXPRESSION DATA FOR THE TCGA-BRCA 
data <- read.table(file = './TCGA-BRCA.htseq_counts.tsv', sep = '\t', header = TRUE)

# READ THE GENE NAMES FOR THE TCGA-BRCA TO MAP GENE IDS WITH GENE NAMES
gene <- read.table(file = './gencode.v22.annotation.gene.probeMap2.txt', sep = '\t', header = TRUE)

# Load the required packages
if (!require(CRDscore)) install.packages("CRDscore", dependencies = TRUE); suppressPackageStartupMessages(library(CRDscore))

data("circadians_selected")
if (!require(tidyverse)) install.packages("tidyverse", dependencies = TRUE); suppressPackageStartupMessages(library(tidyverse))

CRG <- circadians_selected # SELECT THE CRG FROM THE PAPER
gene_name <- gene %>% dplyr::select(id, gene) # SELECT ONLY THE GENE IDs AND GENE NAMES
# RENAME Ensemble_ID to id for joining data with gene_name data
data.new <- data %>% 
  plyr::rename(c("Ensembl_ID" = "id"))
# MERGE THE DATA
data_merged <- left_join(gene_name, data.new, by = "id")
data_merged <- data_merged %>% dplyr::select(-id) # REMOVE variable id
data_merged.new <- data_merged %>% dplyr::select(-gene) # REMOVE gene


# CHANGE data frame to matrix 
data_merged.new = as.matrix(data_merged.new)
rownames(data_merged.new) <- NULL
rownames(data_merged.new) = data_merged$gene

# Calculate Tumor vs Normal Patients
patients <- colnames(data_merged)
Patients <- data.frame(patients = patients, sample = substr(patients, 14, 15) ) 
Patients <- Patients[-1, ]
Patients <- Patients %>% mutate(Tissue = ifelse(sample < 10, "Tumor", "Normal"))

# SELECT THE CLOCK GENES OF INTEREST
data_merged.new.Clock = data_merged.new[rownames(data_merged.new) %in% c("PER1", "PER2", "PER3", "CRY1", "CRY2", "BMAL1", "RORB", "RORC", "NR1D1", "NRID2", "ARNT1", "ARNT2", "CLOCK"), ]

tumor.patients = Patients %>% filter(Tissue == "Tumor") %>% select(patients) %>% pull()
normal.patients = Patients %>% filter(Tissue != "Tumor") %>% select(patients) %>% pull()
data_merged.new.Clock_tumor = data_merged.new.Clock[ ,colnames(data_merged.new.Clock) %in% tumor.patients]
data_merged.new.Clock_normal = data_merged.new.Clock[ ,colnames(data_merged.new.Clock) %in% normal.patients ]

# CALCULATE THE MEAN EXPRESSION FOR TUMOR AND NORMAL PATIENTS
mean_expression_cancer = apply(data_merged.new.Clock_tumor, 1, mean)
mean_expression_normal = apply(data_merged.new.Clock_normal, 1, mean)

# CALCULATE THE SE OF MEAN EXPRESSION FOR TUMOR AND NORMAL PATIENTS
SE_expression_cancer = apply(data_merged.new.Clock_tumor, 1, function(x){sd(x)/sqrt(length(x))})
SE_expression_normal = apply(data_merged.new.Clock_normal, 1, function(x){sd(x)/sqrt(length(x))})
################################################################################
# T-TEST
################################################################################
# CALCULATE THE P-VALUES FOR THE T TEST COMPARING THE MEANS OF TUMOR AND NORMAL PATIENTS
p.values = 0
for(i in 1:nrow(data_merged.new.Clock_normal)){
  p.values[i] = t.test(data_merged.new.Clock_tumor[i, ], data_merged.new.Clock_normal[i, ], var.equal = FALSE)$p.value
}

P = data.frame(Gene = rownames(data_merged.new.Clock_normal), p = p.values, Sig = ifelse(p.values < 0.01, "red3", "black"))

# PREPARE THE DATA FRAME FOR PRODUCING THE DESIRED PLOT
DF <- rbind(data.frame(cells = "Surrounding Tissue", Count = mean_expression_normal, SE =  SE_expression_normal, Genes = names(mean_expression_normal)),
            data.frame(cells = "Tumor", Count = mean_expression_cancer, SE =  SE_expression_cancer, Genes = names(mean_expression_cancer)))
DF$cells <- as.factor(DF$cells)

# DEFINE THE SIGNIFICANCE LEVELS 
Signif = ifelse(P$p < 0.05 & P$p > 0.01,"*", 
                ifelse(P$p < 0.01 & P$p > 0.001, "**",
                       ifelse(P$p < 0.001, "***", "n.s")))

DF <- DF %>% mutate(Signif = c(Signif, rep(NA, length(P$p))))
################################################################################
# PLOT THE DE GENES BY T-TEST
################################################################################
plot1 <- ggplot(DF, aes(x=Genes, y=Count, fill=cells)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("grey", "red3")) + 
  labs(x = "", y = "Normalized Count") + 
  theme_classic() + 
  scale_x_discrete(labels = sort(P$Gene)) + 
  scale_y_continuous(limits = c(0, 15), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1,
                                   colour = "black"),
        legend.position = "top",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 14, face="bold"),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size=18, face="bold", colour = "black"),) + 
  guides(fill = guide_legend(title = "")) +
  geom_errorbar(aes(x=Genes, ymin=Count-SE, ymax=Count+SE), width=0.2, colour="black", alpha=1, size=0.8,
                position = position_dodge(.9)) +
  geom_text(data=DF, aes(x=Genes, y= Count + SE + 1, label=Signif),  
            position = position_dodge(0.9), 
            size=5, vjust=-1, hjust = -0.01, fontface="bold")  


################################################################################
# WILCOXON TEST
################################################################################
# CALCULATE THE P-VALUES FOR THE WILCOXON TEST COMPARING THE MEANS OF TUMOR AND NORMAL PATIENTS

p.values = 0
for(i in 1:nrow(data_merged.new.Clock_normal)){
  p.values[i] = wilcox.test(data_merged.new.Clock_tumor[i, ], data_merged.new.Clock_normal[i, ])$p.value
}

P = data.frame(Gene = rownames(data_merged.new.Clock_normal), p = p.values, Sig = ifelse(p.values < 0.01, "red3", "black"))

# PREPARE THE DATA FRAME FOR PRODUCING THE DESIRED PLOT
DF <- rbind(data.frame(cells = "Surrounding Tissue", Count = mean_expression_normal, SE =  SE_expression_normal, Genes = names(mean_expression_normal)),
            data.frame(cells = "Tumor", Count = mean_expression_cancer, SE =  SE_expression_cancer, Genes = names(mean_expression_cancer)))
DF$cells <- as.factor(DF$cells)

# DEFINE THE SIGNIFICANCE LEVELS 
Signif = ifelse(P$p < 0.05 & P$p > 0.01,"*", 
                ifelse(P$p < 0.01 & P$p > 0.001, "**",
                       ifelse(P$p < 0.001, "***", "n.s")))

DF <- DF %>% mutate(Signif = c(Signif, rep(NA, length(P$p))))

################################################################################
# PLOT THE DE GENES BY WILCOXON TEST
################################################################################
plot2 <- ggplot(DF, aes(x=Genes, y=Count, fill=cells)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("grey", "red3")) + 
  labs(x = "", y = "Normalized Count") + 
  theme_classic() + 
  scale_x_discrete(labels = sort(P$Gene)) + 
  scale_y_continuous(limits = c(0, 15), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1,
                                   colour = "black"),
        legend.position = "top",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 14, face="bold"),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size=18, face="bold", colour = "black"),) + 
  guides(fill = guide_legend(title = "")) +
  geom_errorbar(aes(x=Genes, ymin=Count-SE, ymax=Count+SE), width=0.2, colour="black", alpha=1, size=0.8,
                position = position_dodge(.9)) +
  geom_text(data=DF, aes(x=Genes, y= Count + SE + 1, label=Signif),  
            position = position_dodge(0.9), 
            size=5, vjust=-1, hjust = -0.01, fontface="bold") 
################################################################################
# SAVE THE PLOTS
################################################################################
ggsave("./t_test.pdf", plot1, width = 8, height = 6, units = "in")
ggsave("./wilcoxon_test.pdf", plot2, width = 8, height = 6, units = "in")
