# Load packages
library(dplyr)
library(ggplot2)
library(pROC)
library(stringr)
library(tidyr)
library(ggpattern)
library(ggridges)
library(ggrepel)

##### Part 1: Analyze PTAFR expression and select samples with high and low PAFR

metadata <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/patient_counts_processed/metadata.tsv", header = TRUE, sep = "\t")

# I will downlod one table as an example
table <- read.table("~/Desktop/Patients_data/files_counts/0a0b75ab-0e92-4d25-abd0-5ebc075b432b.rna_seq.augmented_star_gene_counts.tsv", header = T, sep ="\t", skip = 1)
table <- table[,c(1:3, 7)]
colnames(table) <- c("gene_id", "gene_name", "gene_type", "tpm" )
table <- subset(table, gene_type == "protein_coding")

# Set the directory and get unique tissue names
directory <- "~/Desktop/Patients_data/files_counts"
tissues <- unique(metadata$tissue)

# Function to create a table for each tissue
table_for_tissue <- function(x, metadata, directory) {
    metadata_files <- metadata %>%
        filter(tissue == x)
    
    # List to store tables
    tables_for_tissue <- list()
    
    # Loop on metadata for each tissue
    for (file_name in metadata_files$File.Name) {
        file_path <- file.path(directory, file_name)
        
        # Check if the file exists before attempting to read it
        if (file.exists(file_path)) {
            file <- read.table(file_path, header = TRUE, sep = "\t")
            
            # Extract relevant columns (gene_name, gene_type, tpm_unstranded) and add to the list
            filtered_data <- file %>%
                select(gene_name, tpm_unstranded) %>%
                filter(gene_name == "PTAFR")
            
            tables_for_tissue[[file_name]] <- filtered_data
            cat("File successfully downloaded:", file_path, "\n")
        } else {
            cat("File not found:", file_path, "\n")
        }
    }
    
    return(tables_for_tissue)
}

# List to store tables for each tissue
table_list <- list()

# Loop over the unique tissues to create separate tables
for (tissue_name in tissues) {
    table_list[[tissue_name]] <- table_for_tissue(tissue_name, metadata, directory)
}

file_names <- data.frame(t(bind_rows(table_list)))
file_names <- file_names[!grepl("gene_name", rownames(file_names)), ]
file_names <- rownames(file_names)

combined_dataframe <- data.frame(t(bind_cols(table_list)))
combined_dataframe <- combined_dataframe[!grepl("gene_name", rownames(combined_dataframe)), ]

PTAFR_counts <- data.frame(PTAFR = combined_dataframe,
                           File.Name = file_names)

PTAFR_counts$PTAFR <- as.numeric(PTAFR_counts$PTAFR)
PTAFR_counts <- PTAFR_counts %>%
    mutate_all(~ gsub(".tpm_unstranded", "", .))

PTAFR_counts_data <- merge(PTAFR_counts, metadata, by = "File.Name")
PTAFR_counts_data$PTAFR <- as.numeric(PTAFR_counts_data$PTAFR)

# Organize tissue informations
PTAFR_counts_data <- PTAFR_counts_data %>%
	filter(!is.na(tissue) & tissue != "intestine" & tissue != "adrenal gland")

PTAFR_counts_data$tissue <- sub("esophagus", "stomach", PTAFR_counts_data$tissue)
PTAFR_counts_data$tissue <- sub("rectum", "colorectal", PTAFR_counts_data$tissue)
PTAFR_counts_data$tissue <- sub("colon", "colorectal", PTAFR_counts_data$tissue)
PTAFR_counts_data$tissue <- sub("bile duct", "liver", PTAFR_counts_data$tissue)

# Calcular os quartis
q1 <- quantile(PTAFR_counts_data$PTAFR, 0.25, na.rm = TRUE)
q2 <- quantile(PTAFR_counts_data$PTAFR, 0.50, na.rm = TRUE)
q3 <- quantile(PTAFR_counts_data$PTAFR, 0.75, na.rm = TRUE)

# Histogram
# Configure PTAFR expression table
PTAFR_neg <- subset(PTAFR_counts_data, PTAFR <= q1)
PTAFR_neg$type <- c(rep("Negative", nrow(PTAFR_neg)))

PTAFR_pos <- subset(PTAFR_counts_data, PTAFR >= q3)
PTAFR_pos$type <- c(rep("Positive", nrow(PTAFR_pos)))

PTAFR_neutral <- subset(PTAFR_counts_data, PTAFR <= q3 & PTAFR >= q1)
PTAFR_neutral$type <- c(rep("Control", nrow(PTAFR_neutral)))

PTAFR_expression <- rbind(PTAFR_neg, PTAFR_pos, PTAFR_neutral)
PTAFR_expression$tissue <- str_to_title(PTAFR_expression$tissue)
PTAFR_expression$tissue <- as.factor(PTAFR_expression$tissue)


# Histogram all
ggplot(PTAFR_expression, aes(x = PTAFR, fill = tissue)) +
	geom_histogram(binwidth = 5, position = "identity", alpha = 0.5) +
	geom_vline(xintercept = q2, color = "black", size = 0.25, linetype = "dashed") +
	geom_vline(xintercept = q1, color = "black", size = 0.25) +
	geom_vline(xintercept = q3, color = "black", size = 0.25) +
	coord_cartesian(xlim = c(0, 100)) +
	labs(x = "PTAFR", y = NULL, fill = "Tissue") +
	theme(
		axis.title.x = element_text(hjust = 0.5),
		axis.title.y = element_blank(),
		axis.text.y = element_text(size = 12),
		legend.position = "none",
		panel.border = element_rect(color = "black", size = 1, fill = NA),
		panel.background = element_rect(fill = "white"))

# Individual histogram
ggplot(PTAFR_expression, aes(x = PTAFR, y = reorder(tissue, desc(tissue)), fill = tissue)) +
	geom_density_ridges(
		scale = 1.5, rel_min_height = 0.05,
		alpha = 0.5,
		jittered_points = TRUE,
		point_shape = "|",
		point_size = 3, size = 0,
		position = position_points_jitter(height = 0),
		bandwidth = 9) +
	coord_cartesian(xlim = c(0, 100)) +
	geom_vline(xintercept = q2, color = "black", size = 0.25, linetype = "dashed") +
	geom_vline(xintercept = q1, color = "black", size = 0.25) +
	geom_vline(xintercept = q3, color = "black", size = 0.25) +
	theme(
		axis.title.x = element_text(hjust = 0.5),
		axis.title.y = element_blank(),
		axis.text.y = element_text(size = 12),
		legend.position = "none",
		panel.border = element_rect(color = "black", size = 1, fill = NA),
		panel.background = element_rect(fill = "white"))

# PTAFR postive and negative for each type of tumor
count <- PTAFR_expression %>%
	filter(type != "Control")

count <- count %>%
	group_by(tissue, type) %>% 
	summarise(count = n())

count <- count %>%
	group_by(tissue) %>%
	mutate(sum = sum(count))

count_percentage <- count
count_percentage$perc <- count_percentage$count / count_percentage$sum

count_percentage <- count_percentage %>%
	select(-count, -sum)

count_percentage$tissue <- as.factor(count_percentage$tissue)
count_percentage$type <- as.factor(count_percentage$type)
count_percentage$perc <- count_percentage$perc * 100
count_percentage$perc <- round(count_percentage$perc)
count_percentage$type <- gsub("Negative", "Low", count_percentage$type)
count_percentage$type <- gsub("Positive", "High", count_percentage$type)

# Plot
ggplot(count_percentage, aes(x = reorder(tissue, desc(tissue)), y = perc, fill = tissue, alpha = type, label = perc)) +
	geom_col(position = "stack") +
	geom_text_repel(aes(label = perc, y = perc), position = position_stack(vjust = 0.5), hjust = 0.5, show.legend = FALSE, size = 5) +
	labs(x = NULL, y = "% samples", alpha = "Type") +
	guides(fill = FALSE) + 
	scale_alpha_manual(values = c(Low = 0.2, High = 0.8)) +
	theme_minimal() +
	theme(
		axis.text.x = element_blank(),
		axis.text.y = element_blank(),
		axis.title.x = element_text(size = 15),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank()
	) +
	coord_flip()


# save metadata table
metadata_ptafr <- PTAFR_expression
dir_path <- "/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/patient_counts_processed/"

write.table(metadata_ptafr, file = paste0 (dir_path, "metadata_ptafr.tsv"), sep = "\t", row.names = FALSE)

# survival rates analysis
survival <- data_frame(tissue = c("Breast", "Uteri", "Colorectal", "Endometrium", "Stomach", "Kidney", "Liver", "Lung", "Pancreas", "Prostate", "Thyroid", "Bladder"),
				   survival = c(90.8, 67.2, 65.0,81.0,28.7,77.6,21.6,25.4,12.5,97.1,98.5,77.9))

sample_size <- PTAFR_expression %>%
	filter(type != "Control") %>%
	group_by(tissue, type) %>%
	summarize(count = n())

sample_size <- sample_size %>%
	group_by(tissue) %>%
	mutate(sample = sum(count))

sample_size <- sample_size[,c(-2,-3)]
sample_size <- unique(sample_size)
	
ptafr_rates <- PTAFR_expression %>%
	group_by(tissue) %>%
	summarize(mean = mean(PTAFR),
			sd_value = sd(PTAFR))

ptafr_rates <- merge(ptafr_rates, sample_size, by = "tissue")

survival <- merge(survival, ptafr_rates, by = "tissue")
survival <- survival %>%
	filter(tissue != "Liver") 

# Crie um gr¨¢fico de dispers??o com a linha de regress??o e mantenha a linha interna do geom_smooth
ggplot(data = survival, aes(x = survival, y = mean)) +
	geom_point(aes(size = sample)) +  
	geom_smooth(method = "lm") +  # Use a regress??o linear para o geom_smooth
	geom_text(aes(label = tissue, size = 100), hjust = -0.2, vjust = -2) +  # Adicione r¨®tulos com nomes de tecido
	labs(x = "Tumor 5-Year Relative Survival (%)", y = "PAFR expression", title = paste0("Correlation test =", cor(survival$survival, survival$mean))) +
	ggthemes::theme_few()


##### Part 2: Extract selected data from metadata_ptafr from different tumor types
# Load packages
library(dplyr)
library(ggplot2)
library(pROC)
library(stringr)
library(tidyr)
library(ggpattern)
library(ggridges)

# Load data
metadata_ptafr <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/patient_counts_processed/metadata_ptafr.tsv", header = TRUE, sep = "\t")

metadata <- metadata_ptafr[metadata_ptafr$type != "Control", ]

# Set the directory and get unique tissue names
directory <- "~/Desktop/Patients_data/files_counts"
tissues <- unique(metadata$tissue)

# Function to create a table for each tissue
table_for_tissue <- function(x, metadata, directory) {
	metadata_files <- metadata %>%
		filter(tissue == x)
	
	# List to store tables
	tables_for_tissue <- list()
	
	# Loop on metadata for each tissue
	for (file_name in metadata_files$File.Name) {
		file_path <- file.path(directory, file_name)
		
		# Check if the file exists before attempting to read it
		if (file.exists(file_path)) {
			file <- read.table(file_path, header = TRUE, sep = "\t", skip = 1)
			
			# Extract relevant columns (gene_type, tpm_unstranded) and add to the list
			filtered_data <- file %>%
				select(gene_type, tpm_unstranded) %>%
				filter(!is.na(tpm_unstranded)) %>%
				filter(gene_type == "protein_coding") %>%
				select(tpm_unstranded)
				
			tables_for_tissue[[file_name]] <- filtered_data
			cat("File sucessfully found:", file_path, "\n")
		} else {
			cat("File not found:", file_path, "\n")
		}
	}
	
	return(tables_for_tissue)
}

# List to store tables for each tissue and loop
table_list <- list() # Here I successfully downloaded the data
for (tissue_name in tissues) {
	table_list[[tissue_name]] <- table_for_tissue(tissue_name, metadata, directory)
}

# First, I'll try to rename tpm_unstranded according to file name
# Loop over the tissues in the table_list and rename 'tpm_unstranded'
table_list <- lapply(table_list, function(tissue_df) {
	lapply(names(tissue_df), function(df_name) {
		df <- tissue_df[[df_name]]
		colnames(df) <- df_name
		return(df)
	})
})

# Inicialize o combined_dataframe com o mesmo n¨²mero de linhas que os dataframes individuais
combined_dataframe <- data.frame(matrix(NA, nrow = nrow(table_list[[1]][[1]]), ncol = 0))

# Loop sobre os dataframes nas listas de listas
for (tissue_df in table_list) {
	# Loop sobre os dataframes em tissue_df
	for (df in tissue_df) {
		# Adicione o dataframe atual como uma nova coluna no combined_dataframe
		combined_dataframe <- cbind(combined_dataframe, df)
		cat("File successfully downloaded:")
	}
}

# Remova as colunas duplicadas (caso existam)
combined_dataframe <- combined_dataframe[, !duplicated(colnames(combined_dataframe))] ### Here I sucessfully separated the selected cancer files

## Create a separated metadata table for each tissue type
for (i in metadata$tissue) {
	filtered_df <- metadata %>%
		filter(tissue == i) %>%
		select(File.Name, case_id, tissue, type)
	
	assign(paste0("metadata.", i), filtered_df, envir = .GlobalEnv)
}

# Create a separated counts table for each tissue type
tissues <- unique(metadata$tissue)
metadata_names <- paste("metadata.", tissues, sep = "")

for (i in metadata_names) {
	current_df <- get(i)
	col_name <- as.character(current_df$File.Name)
	selected_col <- combined_dataframe[col_name]
	assign(paste0("cts_", i), selected_col, envir = .GlobalEnv)
}


### Save new dataframes
## counts tables
tissues <- unique(metadata$tissue)
cts_names <- paste0("cts_metadata.", tissues)
dir_path <- "/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/patient_counts_processed/"

for (i in cts_names) {
	table <- get(i)
	write.table(table, file = paste0 (dir_path, i, ".tsv"), sep = "\t", row.names = FALSE)
}

## Metadata tables
## counts tables
tissues <- unique(metadata$tissue)
metadata_names <- paste0("metadata.", tissues)
dir_path <- "/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/patient_counts_processed/"

for (i in metadata_names) {
	table <- get(i)
	write.table(table, file = paste0 (dir_path, i, ".tsv"), sep = "\t", row.names = FALSE)
}

##### Part 3: Extract tables and perform DESEQ analysis
# Load packages
library(dplyr)
library(ggplot2)
library(pROC)
library(stringr)
library(tidyr)
library(ggpattern)
library(ggridges)
library(DESeq2)

# Load data
# Linux
metadata_ptafr <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/patient_counts_processed/metadata_ptafr.tsv", header = TRUE, sep = "\t")

# Drive
metadata_ptafr <- read.delim("/Users/barbaradalmaso/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/patient_counts_processed/metadata_ptafr.tsv", header = TRUE, sep = "\t")

# Cts files
tissues <- unique(metadata_ptafr$tissue)
file_names <- paste0("cts.", tissues, ".tsv")

for (i in file_names) {
	table <- read.delim(paste0("/Users/barbaradalmaso/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/patient_counts_processed/", i), header = TRUE, sep = "\t")
	assign(gsub(".tsv", "", i), table, envir = .GlobalEnv)
}

for (i in file_names) {
	i <- gsub(".tsv", "", i)
	table <- get(i)
	colnames(table) <- gsub(".rna_seq.augmented_star_gene_counts.tsv", "", colnames(table))
	colnames(table) <- gsub("X", "", colnames(table))
	colnames(table) <- gsub("\\.", "-", colnames(table))
	assign(i, table, envir = .GlobalEnv)
}

# Metadata files
tissues <- unique(metadata_ptafr$tissue)
file_names <- paste0("metadata.", tissues, ".tsv")
for (i in file_names) {
	table <- read.delim(paste0("/Users/barbaradalmaso/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/patient_counts_processed/", i), header = TRUE, sep = "\t")
	table$File.Name <- gsub(".rna_seq.augmented_star_gene_counts.tsv", "", table$File.Name)
	assign(gsub(".tsv", "", i), table, envir = .GlobalEnv)
}

# I will downlod one table as an example
table <- read.table("/Users/barbaradalmaso/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/05ea028e-6b8a-4352-88aa-057494cf5e80.rna_seq.augmented_star_gene_counts.tsv", header = T, sep ="\t", skip = 1)
table <- table[,c(1:3, 7)]
colnames(table) <- c("gene_id", "gene_name", "gene_type", "tpm" )
table <- subset(table, gene_type == "protein_coding")

# Transform cts tables on matrix
tissues <- unique(metadata_ptafr$tissue)
file_names <- paste0("cts.", tissues)

for (i in file_names) {
	file <- get(i)
	file <- as.matrix(file) 
	rownames(file) <- table$gene_name
	assign(i, file, envir = .GlobalEnv)
}

# Get coldata
rownames(metadata.Colorectal) <- colnames(cts.Colorectal)
rownames(metadata.Endometrium) <- colnames(cts.Endometrium)
rownames(metadata.Kidney) <- colnames(cts.Kidney)
rownames(metadata.Liver) <- colnames(cts.Liver)
rownames(metadata.Lung) <- colnames(cts.Lung)
rownames(metadata.Pancreas) <- colnames(cts.Pancreas)
rownames(metadata.Prostate) <- colnames(cts.Prostate)
rownames(metadata.Stomach) <- colnames(cts.Stomach)
rownames(metadata.Thyroid) <- colnames(cts.Thyroid)
rownames(metadata.Uteri) <- colnames(cts.Uteri)

# Deseq analysis
metadata_names <- paste0("metadata.", tissues)
cts_names <- paste0("cts.", tissues)
metadata_names <- sort(metadata_names)
cts_names <- sort(cts_names)
for (i in seq_along(cts_names)) {
	cts <- get(cts_names[i])
	metadata <- get(metadata_names[i])

		if (!is.null(cts) && !is.null(metadata)) {
		res <- DESeqDataSetFromMatrix(countData = round(cts),
								colData = metadata,
								design = ~type)
		deseq_name <- paste0("deseq.", tissues[i])
		assign(deseq_name, res, envir = .GlobalEnv)
		cat("DESeqDataSet criado para", tissues[i], "\n")
	} else {
		cat("Erro: Nao foi poss¨ªvel carregar 'cts.' ou 'metadata.' para", tissues[i], "\n")
	}
}

file_names <- paste0("deseq.", tissues)
for (i in file_names) {
	file <- get(i)
	res <- DESeq(file)
	res <- results(res, contrast = c("type", "Positive", "Negative"))
	res <- data.frame(res)
	res$diffexpressed <- "no"
	res$diffexpressed[res$log2FoldChange > 1 & res$padj < 0.05] <- "up"
	res$diffexpressed[res$log2FoldChange < -1 & res$padj < 0.05] <- "down"
	assign(i, res, envir = .GlobalEnv)
	}

# Volcano-plot (general overview)
file_names <- paste0("deseq.", tissues)
file_names <- sort(file_names)
tissues <- sort(tissues)

for (variable in file_names) {

	file <- get(variable) 
	file <- file[row.names(file) != "PTAFR", ]
	
	plot <- ggplot(data = file,
				aes(x = log2FoldChange,
				    y = -log10(padj),
				    col = diffexpressed)) +
		geom_point() +
		labs(title = variable) +
		scale_color_manual(values = c("#1E90FF", "#a9a9a9", "#FF0000")) +
		theme(panel.background = element_blank(),
			 plot.background = element_blank(),
			 panel.border = element_rect(color = "black", fill = NA, size = 2),
			 legend.position = "none",
			 axis.text.x = element_text(size = 18), 
			 axis.text.y = element_text(size = 18),
			 axis.title.y = element_blank(),
			 axis.title.x = element_blank())
	
	print(plot)
	}

# Count % DEG
file_names 
df <- data.frame()
for (i in file_names) {
	file <- get(i)
	gene_counts <- data.frame(table(file$diffexpressed))
		if (ncol(df) == 0) {
		df <- gene_counts
		colnames(df) <- c("Genes", i)
	} else {
		df[[i]] <- gene_counts$Freq
	}
}
rownames(df) <- df$Genes
df <- df[,-1]
DEG <- data.frame(t(df))
total <- rowSums(DEG)
DEG$total <- total

DEG <- DEG %>%
	mutate(perc.up = 100*up/total,
		  perc.down = 100*down/total)

DEG <- DEG[,c(-1,-2,-3,-4)]
DEG$tissue <- c("Endometrium", "Breast", "Pancreas", "Liver",
			 "Thyroid", "Bladder", "Kidney", "Prostate", 
			 "Colorectal", "Lung", "Uteri", "Stomach")

DEG <- DEG %>%
	arrange(tissue)

DEG_all <- data.frame(perc = c(DEG$perc.up, DEG$perc.down),
				  tissue = c(DEG$tissue, DEG$tissue),
				  type = c(rep.int("up", 12), rep.int("down", 12)))
	
# Plot
ggplot(DEG_all, aes(x=tissue, y = perc, fill = tissue, alpha = type)) +
	geom_col(position = "stack") +
	labs(x = NULL, y = "% DEG", alpha = "Type") +
	guides(fill = FALSE) + 
	scale_alpha_manual(values = c(down = 0.2, up = 0.8)) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

	
# PCA-plot
metadata_names <- paste0("metadata.", tissues)
cts_names <- paste0("cts.", tissues)
metadata_names <- sort(metadata_names)
cts_names <- sort(cts_names)

for (i in cts_names) {
	for (o in metadata_names) {
		cts <- get(i)
		coldata <- get(o)
		
		if (ncol(cts) == nrow(coldata)) { 
			
			# Crie um DESeqDataSetFromMatrix
			res_data <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~type)
			
			# Aplique a transforma????o VST
			res_data <- vst(res_data, blind = FALSE)
			
			# Crie um grafico PCA
			plot <- plotPCA(res_data, intgroup = "type") +
				geom_point(aes(color = type)) + 
				labs(title = i) +
				scale_color_manual(values = c("#00DAE0", "#FF9289")) +
				theme(panel.background = element_blank(),
					 plot.background = element_blank(),
					 panel.border = element_rect(color = "black", fill = NA, size = 2),
					 legend.position = "none",
					 axis.text.x = element_text(size = 18), 
					 axis.text.y = element_text(size = 12),
					 axis.title.y = element_text(size = 18),
					 axis.title.x = element_text(size = 18))
			
			print(plot)
			
		} else {
			print(paste0(i, " n??o tem o mesmo n¨²mero de colunas que ", o))
		}
	}
}

# Radsets
library(radsets)
library(circlize)
library(gplots)
library(purrr)
library(RVenn)
library(ggplot2)

# Get sig up genes 
file_names <- paste0("deseq.", tissues)
for (i in file_names) {
	data <- get(i)
	data <- as.data.frame(data)
	data <- data %>%
		filter(padj <= 0.05, log2FoldChange >= 1)
	data_genes <- rownames(data)
	assign(paste0("up.", i), data_genes, envir = .GlobalEnv)
}

### Vennn
file_names <- paste0("up.", file_names)
new_list <- list()
for (i in file_names) {
	list <- list(get(i))
	tissue <- gsub("up.deseq.", "", i)
	names(list) <- tissue
	new_list <- append(new_list,list)
}

toy <- Venn(new_list)
setmap(toy, element_clustering = T, set_clustering = T)
sig_toy <- overlap(toy)

# Get sig down genes 
file_names <- paste0("deseq.", tissues)
for (i in file_names) {
	data <- get(i)
	data <- as.data.frame(data)
	data <- data %>%
		filter(padj <= 0.05, log2FoldChange <= -1)
	data_genes <- rownames(data)
	assign(paste0("down.", i), data_genes, envir = .GlobalEnv)
}

file_names <- paste0("down.", file_names)
new_list2 <- list()
for (i in file_names) {
	list <- list(get(i))
	tissue <- gsub("down.deseq.", "", i)
	names(list) <- tissue
	new_list2 <- append(new_list2,list)
}

toy2 <- Venn(new_list2)
setmap(toy2, element_clustering = T, set_clustering = T)

library(clusterProfiler)
library(pathview)
library(enrichplot)
library(org.Hs.eg.db)

# ORA Analysis for all
genes <- enrichGO(gene = sig_toy,
	    OrgDb = org.Hs.eg.db, keyType = "SYMBOL", readable = TRUE,
	    ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

barplot(genes, drop = TRUE, showCategory = 15, title = "GO Biological Pathways",
	   font.size = 5) +
	ggtitle("Bar-plot with raw GO enrichment data")

# Fazer o c¨¢lculo da m¨¦dia de log2FoldChange entre as amostras
file_names <- paste0("deseq.", tissues)
fc_all <- data.frame()

for (i in file_names) {
	file <- get(i)
	if (length(fc_all) == 0) {
		fc_all <- file[, 1:2]
		colnames(fc_all) <- c("mean", gsub("deseq.", "", i))
	} else {
		colname <- gsub("deseq.", "", i) # Corrigido: extrair o nome da coluna corretamente
		fc_all[, colname] <- file$log2FoldChange
		cat(i, "adicionado\n") # Adicionado "\n" para uma nova linha na mensagem de sa¨ªda
	}
}

fc_all <- fc_all[,-1]

fd_siggenes <- fc_all[rownames(fc_all) %in% sig_toy, ]
fd_siggenes$mean <- rowMeans(fd_siggenes)

gene_list <- fd_siggenes$mean
names(gene_list) <- rownames(fd_siggenes)
cnetplot(genes, circular = TRUE,showCategory = 4, colorEdge = T, foldChange=gene_list) +
	scale_colour_gradient2(name = "fold change", low = "white", mid = "green", high = "red")

# ORA Analysis individual
file_names <- paste0("deseq.", tissues)
for (i in file_names) {
	data <- get(i)
	data <- as.data.frame(data)
	data <- data %>%
		filter(padj <= 0.05, log2FoldChange >= 1 | log2FoldChange < -1)
	data_genes <- rownames(data)
	assign(paste0("sig.", i), data_genes, envir = .GlobalEnv)
}

file_names <- paste0("deseq.", tissues)
file_names <- paste0("sig.", file_names)
for (i in file_names) {
	file <- get(i)
genes <- enrichGO(gene = file,
			   OrgDb = org.Hs.eg.db, keyType = "SYMBOL", readable = TRUE,
			   ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

genes <- data.frame(genes)
assign(paste0("ora", i), genes, envir = .GlobalEnv)
}

file_names <- paste0("ora", file_names)
df <- data.frame(Description = character(0), padjust = numeric(0), Tissue = character(0))
for (i in file_names) {
	file <- get(i)
	file <- file[1:3, c(3, 7)]
	file$Tissue <- rep(gsub("orasig.deseq.", "", i), nrow(file))
	if (nrow(df) == 0) {
		df <- file
	} else {
		df <- rbind(df, file)
	}
}

df <- df %>% 
	mutate(Description = gsub("extracellular structure organization", 
						 "extracellular matrix organization", Description),
		  Description = gsub("leukocyte chemotaxis", 
		  			    "leukocyte migration", Description),
		  Description = gsub("positive regulation of cell adhesion", 
		  			    "leukocyte cell-cell adhesion", Description),
		  Description = gsub("regulation of cell-cell adhesion", 
		  			    "leukocyte cell-cell adhesion", Description),
		  Description = gsub("organic acid catabolic process", 
		  			    "small molecule catabolic process", Description),
		  p.adjust = -log2(p.adjust))

cores <- c("#FFF303", "#FEBE11", "#F6861E", "#D32027",
		 "#F9B1A0", "#C44F9E", "#7162AB", "#313A5E",
		 "#52A3C0", "#75A240", "#056839", "#443E17")

ggplot(df, aes(x = Tissue, y = p.adjust, fill = Description)) +
	geom_col(position = "dodge") +
	coord_polar() + 
	theme_minimal() +
	scale_fill_manual(values = cores)

##### Part 4: Analysis metadados
library(dplyr)
library(ggplot2)
library(pROC)
library(stringr)
library(tidyr)
library(ggpattern)
library(ggridges)
library(DESeq2) 

# Complete metadata
clinical <- read.delim("~/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/clinical_data/clinical.tsv", header = TRUE, sep = "\t")
pathology <- read.delim("~/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/clinical_data/pathology_detail.tsv", header = TRUE, sep = "\t")
slide <- read.delim("~/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/sample_data/slide.tsv", header = TRUE, sep = "\t")
sample <- read.delim("~/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/gdc_sample_sheet.2023-08-09.tsv", header = TRUE, sep = "\t")

# Breast e bladder metadata
# Breast 
clinical_breast <- read.delim("~/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/metadata_breast/clinical.tsv", header = TRUE, sep = "\t")
pathology_breast <- read.delim("~/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/metadata_breast/pathology_detail.tsv", header = TRUE, sep = "\t")
slide_breast <- read.delim("~/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/metadata_breast/slide.tsv", header = TRUE, sep = "\t")
sample_breast <- read.delim("~/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/metadata_breast/gdc_sample_sheet.2023-09-06.tsv", header = TRUE, sep = "\t")

# Bladder
clinical_bladder <- read.delim("~/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/metadata_bladder/clinical.tsv", header = TRUE, sep = "\t")
pathology_bladder <- read.delim("~/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/metadata_bladder/pathology_detail.tsv", header = TRUE, sep = "\t")
slide_bladder <- read.delim("~/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/metadata_bladder/slide.tsv", header = TRUE, sep = "\t")
sample_bladder <- read.delim("~/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/metadata_bladder/gdc_sample_sheet.2023-09-06.tsv", header = TRUE, sep = "\t")

# Metadata all
clinical <- rbind(clinical, clinical_bladder, clinical_breast)
pathology <- rbind(pathology, pathology_bladder, pathology_breast)
slide <- rbind(slide, slide_bladder, slide_breast)
sample <- rbind(sample, sample_bladder, sample_breast)

# Metadata A
metadata_c <- read.delim("~/Library/CloudStorage/GoogleDrive-barbdalmaso@gmail.com/Meu Drive/DSMZ/Cancer/patient_counts_processed/metadata_ptafr.tsv", header = TRUE, sep = "\t")
metadata_a <- merge(metadata_c, slide[,c(2,16,17,19)], by = "case_id")
pathology_a <- pathology[,c(1,21)]
clinical_a <- clinical[,c(1,27)]
slide_a <- slide[,c(2,16,17,19)]
metadata_a <- merge(metadata_c, pathology_a, by = "case_id", all = TRUE)
metadata_a <- merge(metadata_a, clinical_a, by = "case_id", all = TRUE)
metadata_a <- merge(metadata_a, slide_a, by = "case_id", all = TRUE)
metadata_a <- metadata_a[,-1]
metadata_a <- metadata_a[!duplicated(metadata_a),]
metadata_a <- metadata_a[!is.na(metadata_a$File.Name), ]

metadata_a$ajcc_pathologic_n <- gsub("'--", NA, metadata_a$ajcc_pathologic_n)
metadata_a$ajcc_pathologic_n <- gsub("N0 \\(i-\\)|N0 \\(i\\+\\)|N0 \\(mol\\+\\)", "N0", metadata_a$ajcc_pathologic_n)
metadata_a$ajcc_pathologic_n <- gsub("N1a|N1b|N1c|N1mi", "N1", metadata_a$ajcc_pathologic_n)
metadata_a$ajcc_pathologic_n <- gsub("N2a|N2b", "N2", metadata_a$ajcc_pathologic_n)
metadata_a$ajcc_pathologic_n <- gsub("N3a|N3b|N3c", "N3", metadata_a$ajcc_pathologic_n)

metadata_a$lymph_nodes_positive <- gsub("'--", NA, metadata_a$lymph_nodes_positive)
metadata_a$lymph_nodes_positive <- as.numeric(metadata_a$lymph_nodes_positive)

metadata_a$percent_lymphocyte_infiltration <- gsub("'--", NA, metadata_a$percent_lymphocyte_infiltration)
metadata_a$percent_monocyte_infiltration <- gsub("'--", NA, metadata_a$percent_monocyte_infiltration)
metadata_a$percent_neutrophil_infiltration <- gsub("'--", NA, metadata_a$percent_neutrophil_infiltration)
metadata_a$percent_lymphocyte_infiltration <- as.numeric(metadata_a$percent_lymphocyte_infiltration)
metadata_a$percent_monocyte_infiltration <- as.numeric(metadata_a$percent_monocyte_infiltration)
metadata_a$percent_neutrophil_infiltration <- as.numeric(metadata_a$percent_neutrophil_infiltration)

# Criar um gr¨¢fico de barras com as m¨¦dias e desvios-padr??o
metadata_a <- metadata_a[metadata_a$type != "Control", ]

cores <- c("Negative" = "lightgray", "Positive" = "darkgray")
ggplot(data = metadata_a, aes(x = type, y = percent_lymphocyte_infiltration, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.2, size = 1) +
	labs(x = "Type", y = "Lymphocyte Infiltration (%)") +
	scale_fill_manual(values = cores) +  # Define as cores manualmente
	theme(panel.background = element_blank(),
		 plot.background = element_blank(),
		 panel.border = element_rect(color = "black", fill = NA, size = 2),
		 legend.position = "none",
		 axis.text.x = element_text(size = 10), 
		 axis.text.y = element_text(size = 10),
		 axis.title.y = element_text(size = 13),
		 axis.title.x = element_blank())

library(stats)
grupo_negative <- metadata_a$percent_lymphocyte_infiltration[metadata_a$type == "Negative"]
grupo_positive <- metadata_a$percent_lymphocyte_infiltration[metadata_a$type == "Positive"]
t.test(grupo_negative, grupo_positive)

ggplot(data = metadata_a, aes(x = type, y = percent_monocyte_infiltration, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.2, size = 1) +
	labs(x = "Type", y = "Monocyte Infiltration (%)") +
	scale_fill_manual(values = cores) +  # Define as cores manualmente
	theme(panel.background = element_blank(),
		 plot.background = element_blank(),
		 panel.border = element_rect(color = "black", fill = NA, size = 2),
		 legend.position = "none",
		 axis.text.x = element_text(size = 10), 
		 axis.text.y = element_text(size = 10),
		 axis.title.y = element_text(size = 13),
		 axis.title.x = element_blank())
grupo_negative <- metadata_a$percent_monocyte_infiltration[metadata_a$type == "Negative"]
grupo_positive <- metadata_a$percent_monocyte_infiltration[metadata_a$type == "Positive"]
t.test(grupo_negative, grupo_positive)

ggplot(data = metadata_a, aes(x = type, y = percent_neutrophil_infiltration, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.2, size = 1) +
	labs(x = "Type", y = "Neutrophil Infiltration (%)") +
	scale_fill_manual(values = cores) +  # Define as cores manualmente
	theme(panel.background = element_blank(),
		 plot.background = element_blank(),
		 panel.border = element_rect(color = "black", fill = NA, size = 2),
		 legend.position = "none",
		 axis.text.x = element_text(size = 10), 
		 axis.text.y = element_text(size = 10),
		 axis.title.y = element_text(size = 13),
		 axis.title.x = element_blank())
grupo_negative <- metadata_a$percent_neutrophil_infiltration[metadata_a$type == "Negative"]
grupo_positive <- metadata_a$percent_neutrophil_infiltration[metadata_a$type == "Positive"]
t.test(grupo_negative, grupo_positive)

# Metadata B
metadata_b <- merge(metadata_c, slide[,c(2,18,20,23,24,25)], by = "case_id")
pathology_b <- pathology[,c(1,48)]
clinical_b <- clinical[,c(1,26,28,29)]
metadata_b <- merge(metadata_b, pathology_b, by = "case_id", all = TRUE)
metadata_b <- merge(metadata_b, clinical_b, by = "case_id", all = TRUE)
metadata_b <- metadata_b[,-1]
metadata_b <- metadata_b[!duplicated(metadata_b$File.Name),]
metadata_b <- metadata_b[!is.na(metadata_b$File.Name), ]

metadata_b$percent_necrosis <- gsub("'--", NA, metadata_b$percent_necrosis)
metadata_b$percent_normal_cells <- gsub("'--", NA, metadata_b$percent_normal_cells)
metadata_b$percent_stromal_cells <- gsub("'--", NA, metadata_b$percent_stromal_cells)
metadata_b$percent_tumor_cells <- gsub("'--", NA, metadata_b$percent_tumor_cells)
metadata_b$percent_tumor_nuclei <- gsub("'--", NA, metadata_b$percent_tumor_nuclei)
metadata_b$tumor_largest_dimension_diameter <- gsub("'--", NA, metadata_b$tumor_largest_dimension_diameter)
metadata_b$ajcc_pathologic_m <- gsub("'--", NA, metadata_b$ajcc_pathologic_m)
metadata_b$ajcc_pathologic_stage <- gsub("'--", NA, metadata_b$ajcc_pathologic_stage)
metadata_b$ajcc_pathologic_t <- gsub("'--", NA, metadata_b$ajcc_pathologic_t)
metadata_b[, 5:10] <- lapply(metadata_b[, 5:10], as.numeric)
metadata_b <- metadata_b[metadata_b$type != "Control", ]


ggplot(data = metadata_b, aes(x = type, y = percent_tumor_cells, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.2, size = 1) +
	labs(x = "Type", y = "Tumor cells (%)") +
	scale_fill_manual(values = cores) +  # Define as cores manualmente
	theme(panel.background = element_blank(),
		 plot.background = element_blank(),
		 panel.border = element_rect(color = "black", fill = NA, size = 2),
		 legend.position = "none",
		 axis.text.x = element_text(size = 10), 
		 axis.text.y = element_text(size = 10),
		 axis.title.y = element_text(size = 13),
		 axis.title.x = element_blank())
grupo_negative <- metadata_b$percent_tumor_cells[metadata_b$type == "Negative"]
grupo_positive <-  metadata_b$percent_tumor_cells[metadata_b$type == "Positive"]
t.test(grupo_negative, grupo_positive)

ggplot(data = metadata_b, aes(x = type, y = percent_normal_cells, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.2, size = 1) +
	labs(x = "Type", y = "Normal cells (%)") +
	scale_fill_manual(values = cores) +  # Define as cores manualmente
	theme(panel.background = element_blank(),
		 plot.background = element_blank(),
		 panel.border = element_rect(color = "black", fill = NA, size = 2),
		 legend.position = "none",
		 axis.text.x = element_text(size = 10), 
		 axis.text.y = element_text(size = 10),
		 axis.title.y = element_text(size = 13),
		 axis.title.x = element_blank())
grupo_negative <- metadata_b$percent_normal_cells[metadata_b$type == "Negative"]
grupo_positive <-  metadata_b$percent_normal_cells[metadata_b$type == "Positive"]
t.test(grupo_negative, grupo_positive)

ggplot(data = metadata_b, aes(x = type, y = percent_stromal_cells, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.2, size = 1) +
	labs(x = "Type", y = "Stromal cells (%)") +
	scale_fill_manual(values = cores) +  # Define as cores manualmente
	theme(panel.background = element_blank(),
		 plot.background = element_blank(),
		 panel.border = element_rect(color = "black", fill = NA, size = 2),
		 legend.position = "none",
		 axis.text.x = element_text(size = 10), 
		 axis.text.y = element_text(size = 10),
		 axis.title.y = element_text(size = 13),
		 axis.title.x = element_blank())
grupo_negative <- metadata_b$percent_stromal_cells[metadata_b$type == "Negative"]
grupo_positive <-  metadata_b$percent_stromal_cells[metadata_b$type == "Positive"]
t.test(grupo_negative, grupo_positive)

ggplot(data = metadata_b, aes(x = type, y = percent_necrosis, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.2, size = 1) +
	labs(x = "Type", y = "Necrosis (%)") +
	scale_fill_manual(values = cores) +  # Define as cores manualmente
	theme(panel.background = element_blank(),
		 plot.background = element_blank(),
		 panel.border = element_rect(color = "black", fill = NA, size = 2),
		 legend.position = "none",
		 axis.text.x = element_text(size = 10), 
		 axis.text.y = element_text(size = 10),
		 axis.title.y = element_text(size = 13),
		 axis.title.x = element_blank())
grupo_negative <- metadata_b$percent_necrosis[metadata_b$type == "Negative"]
grupo_positive <-  metadata_b$percent_necrosis[metadata_b$type == "Positive"]
t.test(grupo_negative, grupo_positive)

ggplot(data = metadata_b, aes(x = type, y = percent_tumor_nuclei, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.2, size = 1) +
	labs(x = "Type", y = "Necrosis (%)") +
	scale_fill_manual(values = cores) +  # Define as cores manualmente
	theme(panel.background = element_blank(),
		 plot.background = element_blank(),
		 panel.border = element_rect(color = "black", fill = NA, size = 2),
		 legend.position = "none",
		 axis.text.x = element_text(size = 10), 
		 axis.text.y = element_text(size = 10),
		 axis.title.y = element_text(size = 13),
		 axis.title.x = element_blank())
grupo_negative <- metadata_b$percent_tumor_nuclei[metadata_b$type == "Negative"]
grupo_positive <-  metadata_b$percent_tumor_nuclei[metadata_b$type == "Positive"]
t.test(grupo_negative, grupo_positive)

ggplot(data = metadata_b, aes(x = type, y = tumor_largest_dimension_diameter, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.2, size = 1) +
	labs(x = "Type", y = "Tumor diameter (cm)") +
	scale_fill_manual(values = cores) +  # Define as cores manualmente
	theme(panel.background = element_blank(),
		 plot.background = element_blank(),
		 panel.border = element_rect(color = "black", fill = NA, size = 2),
		 legend.position = "none",
		 axis.text.x = element_text(size = 10), 
		 axis.text.y = element_text(size = 10),
		 axis.title.y = element_text(size = 13),
		 axis.title.x = element_blank())
grupo_negative <- metadata_b$tumor_largest_dimension_diameter[metadata_b$type == "Negative"]
grupo_positive <-  metadata_b$tumor_largest_dimension_diameter[metadata_b$type == "Positive"]
t.test(grupo_negative, grupo_positive)
