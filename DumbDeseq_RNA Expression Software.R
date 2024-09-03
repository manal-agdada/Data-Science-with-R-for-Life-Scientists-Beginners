############################################ 
#### DumbDeseq: RNA Expression Software ####
############################################

# set working directory, load and view the dataset
View(exp)

# generate a volcano plot ####
library(ggplot2)
library(tidyverse)

ggplot(data = exp, 
       aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point() + 
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red") + 
  theme_minimal()

# determine the upregulated genes (log2FC > 1 and pvalue< 0.01) and downregulated genes (log2FC < -1 and pvalue < 0.01) ####
exp <- exp %>%
  mutate(significance = case_when(
    log2FoldChange >= 1 & pvalue <= 0.01 ~ 'Upregulated',
    log2FoldChange <= -1 & pvalue <= 0.01 ~ 'Downregulated',
    abs(log2FoldChange) < 1 | pvalue > 0.01 ~ 'Not significant'
  ))

# better volcano plot adding color to genes according to their significance 
ggplot(exp,
      aes(x = log2FoldChange,
          y = -log10(pvalue),
          color = significance) 
) + 
  geom_point() + 
  geom_vline(xintercept = c(-1,1), color = 'red') + 
  geom_hline(yintercept = -log10(0.01), color = 'red') + 
  labs(x = 'Log2FC', y = '-Log10(pvalue)',color = 'Significance') + 
  scale_color_manual(values = c('blue','grey90','red'))+
  theme_minimal() 


# functions of the top 5 upregulated genes and top 5 downregulated genes (use genecards) ####
## COMMENT FOR THE REVIEWER: I couldn't use genecards for this task and it's difficult to get gene information without having the metadata of the dataset

# top 5 upregulated genes 
top5_up <- exp %>% 
  filter(significance == "Upregulated") %>% 
  arrange(desc(log2FoldChange)) %>% 
  head(n = 5)
top5_up_list <- top5_up$Gene

#top 5 downregulated genes
top5_down <- exp %>% 
  filter(significance == "Downregulated") %>% 
  arrange(desc(log2FoldChange)) %>% 
  head(n = 5)
top5_down_list <- top5_down$Gene

library(biomaRt)

# connect to the Ensembl database (I assume the gene expression data is from human tissue)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# retrieve gene information for the top 5 upregulated genes
gene_info_up <- getBM(
  attributes = c("external_gene_name", "description"),
  filters = "external_gene_name",
  values = top5_up_list,
  mart = ensembl
)

print(gene_info_up)

# retrieve gene information for the top 5 downregulated genes
gene_info_down <- getBM(
  attributes = c("external_gene_name", "description"),
  filters = "external_gene_name",
  values = top5_down_list,
  mart = ensembl
)

print(gene_info_down)
