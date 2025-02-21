rm(list = ls())
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("Biostrings", "biomaRt", "GenomicRanges", "BSgenome.Hsapiens.UCSC.hg38"))
# library(Biostrings)
# library(biomaRt)
# library(GenomicRanges)
# library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(biomaRt)
library(rtracklayer)
library(dplyr)
library(readr)
################################################################################
############  整理Blast+ 需要的数据
############ GCA_000001405.28_GRCh38.p13_genomic从网上下载，序列在这个代码中获得
################################################################################
expression_data <- read.csv("data/data used for mapping probe/expression_data.csv", stringsAsFactors = FALSE)
probe_sequences <- read.csv("data/data used for mapping probe/probe_sequences.csv", stringsAsFactors = FALSE)
head(expression_data)
head(probe_sequences)


# 从 expression_data 提取探针 ID，去掉前面的 URL 部分，只保留探针编号 (A-MEXP-255.x)
expression_data$Reporter_id <- sub(".*A-MEXP-255\\.", "", expression_data$`Reporter.identifier`)
# 确保探针编号是数值型，以便与 probe_sequences 进行匹配
expression_data$Reporter_id <- as.numeric(expression_data$Reporter_id)

# 替换无效字符为 "N"（例如 '_' 或其他非 ACGT 字符）
probe_sequences$`Reporter.Sequence` <- gsub("[^ACGT]", "N", probe_sequences$`Reporter.Sequence`)
## 将 NA 替换为一个空序列或者 "NNNNN"
# probe_sequences$`Reporter.Sequence`[is.na(probe_sequences$`Reporter.Sequence`)] <- "NNNNN"


# 根据探针编号与探针序列数据进行合并
merged_data <- merge(expression_data, probe_sequences, by.x = "Reporter_id", by.y = "Reporter.Name", all.x = TRUE)
head(merged_data)

# # 假设 merged_data 是您的数据框
# # 创建一个 FASTA 格式的字符串
# fasta_lines <- c()
# for (i in 1:nrow(merged_data)) {
#   # 提取 Reporter_id 和 Reporter.Sequence
#   sequence_id <- paste0(">", merged_data$Reporter_id[i])  # 创建序列名称
#   sequence <- merged_data$Reporter.Sequence[i]  # 提取序列
#   fasta_lines <- c(fasta_lines, sequence_id, sequence)  # 将名称和序列添加到列表
# }
# writeLines(fasta_lines, "data/data used for mapping probe/sequences.fasta")




################################################################################
############  Blast 运行后得到result.txt
############ 包含了GenBank的信息，qseqid的这列和merged_data$Reporter.identifier
############ 是一个意思，顺序可能不同
################################################################################
blast_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", 
                "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

blast_results <- read_tsv("data/data used for mapping probe/results.txt", 
                          col_names = blast_cols, 
                          comment = "#")
head(blast_results)



blast_results_filtered <- blast_results %>%
  distinct(qseqid, .keep_all = TRUE)
print(blast_results_filtered)

print(expression_data)


blast_results_rearranged <- expression_data %>%
  left_join(blast_results_filtered, by = c("Reporter_id" = "qseqid"))
print(blast_results_rearranged)





# 加载所需的R包
library(dplyr)
library(biomaRt)
# # 下载组装报告文件
# assembly_report_url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt"
# download.file(assembly_report_url, destfile = "data/data_used_for_mapping_probe/assembly_report.txt", method = "auto")

# 读取文件，跳过注释行
assembly_data <- read.table("data/data used for mapping probe/assembly_report.txt", 
                            sep = "\t", comment.char = "#", header = FALSE, 
                            stringsAsFactors = FALSE, fill = TRUE)

# 为数据添加列名
colnames(assembly_data) <- c("Sequence_Name", "Sequence_Role", "Assigned_Molecule", 
                             "Assigned_Molecule_Location_Type", "GenBank_Accn", 
                             "Relationship", "RefSeq_Accn", "Assembly_Unit", "Sequence_Length", "UCSC_style_name")

# 提取GenBank Accession和对应染色体的映射关系
genbank_to_chromosome <- assembly_data[, c("GenBank_Accn", "Sequence_Name")]

# # 示例GenBank访问号列表
# genbank_ids <- c("CM000663.2", "CM000664.2", "CM000665.2")
# 
# # 过滤出感兴趣的GenBank访问号
# genbank_to_chr <- genbank_to_chromosome %>%
#   filter(GenBank_Accn %in% genbank_ids)
# colnames(genbank_to_chr) <- c("genbank_id", "chr")
# 
# # 显示结果
# print(genbank_to_chr)

# 假设您的BLAST结果保存在一个数据框中，名为blast_results
blast_results <- blast_results_rearranged

# 将 GenBank ID 映射为染色体编号
blast_results_mapped <- blast_results %>%
  left_join(genbank_to_chromosome, by = c("sseqid" = "GenBank_Accn"))

##############################################
# 提取感兴趣的列用于基因符号查询
blast_subset <- blast_results_mapped %>%
  select(Reporter.identifier, Sequence_Name, sstart, send) # 

load("./output/differential_results.rda")
set1 <- ACE_result$index
set2 <- PFA_result$index
set3 <- BH_result$index
unique_in_ACE <- setdiff(set1, union(set2, set3))


blast_subset <- blast_subset[unique_in_ACE, ] ##################################
blast_subset <- na.omit(blast_subset)

# 修正起始和结束位置
blast_subset <- blast_subset %>%
  mutate(start = pmin(sstart, send),
         end = pmax(sstart, send))

# 连接到Ensembl
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

# 初始化结果数据框
all_results <- data.frame()

# 循环每一行进行查询
for(i in 1:nrow(blast_subset)) {
  seq_name <- blast_subset$Sequence_Name[i]
  start_pos <- blast_subset$start[i]
  end_pos <- blast_subset$end[i]
  
  # 查询基因
  result <- getBM(
    attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"),
    filters = c("chromosome_name", "start", "end"),
    values = list(seq_name, start_pos, end_pos),
    mart = mart
  )
  
  # 添加查询结果，并保留原始行信息
  if(nrow(result) > 0) {
    result$Reporter.identifier <- blast_subset$Reporter.identifier[i]
    all_results$chromosome_name <- as.character(all_results$chromosome_name)
    result$chromosome_name <- as.character(result$chromosome_name)
    all_results <- bind_rows(all_results, result)
  }
}

# 查看合并后的结果
print(all_results)


# 确保数据类型正确并合并 `start` 和 `end` 列
all_results <- all_results %>%
  left_join(blast_subset %>% select(Reporter.identifier, start, end), 
            by = "Reporter.identifier") %>%
  mutate(
    start_position = as.numeric(start_position),
    end_position = as.numeric(end_position),
    start = as.numeric(start),
    end = as.numeric(end)
  )

# 计算距离
all_results <- all_results %>%
  mutate(distance = abs(start_position - start) + abs(end_position - end))

# 找到距离最小的记录，并处理 hgnc_symbol 为空的情况
final_results <- all_results %>%
  group_by(Reporter.identifier) %>%
  arrange(distance) %>%
  # 填充空的 hgnc_symbol
  mutate(
    hgnc_symbol = ifelse(
      hgnc_symbol == "" | is.na(hgnc_symbol),
      # 查找同一 Reporter.identifier 的其他非空 hgnc_symbol
      first(hgnc_symbol[hgnc_symbol != "" & !is.na(hgnc_symbol)], default = NA),
      hgnc_symbol
    )
  ) %>%
  # 保留距离最小的记录
  slice(1) %>%
  ungroup()

# 查看结果
print(final_results)
write.csv(final_results, "data/gene_name.csv")





################################################################################
###### 找到gene_name中Reporter.identifier在expression_data中的坐标  ############
rm(list = ls())
expression_data <- read.csv("data/data used for mapping probe/expression_data.csv", stringsAsFactors = FALSE)
head(expression_data)
gene_name <- read.csv("data/gene_name.csv")
head(gene_name)


load("./output/differential_results.rda")
set1 <- ACE_result$index
set2 <- PFA_result$index
set3 <- BH_result$index
union_ACE_PFA <- union(set1, set2)
reject_only_by_PFA <- setdiff(union_ACE_PFA, set1) # unique miRNA rejected by PFA (donot consider BH!!!!)
reject_only_by_ACE <- setdiff(union_ACE_PFA, set2) # unique miRNA rejected by ACE (donot consider BH!!!!)
unique_in_ACE <- setdiff(set1, union(set2, set3))


# 找到gene_name中Reporter.identifier在expression_data中的坐标
row_indices <- match(gene_name$Reporter.identifier, expression_data$Reporter.identifier)
is_subset <- all(row_indices %in% unique_in_ACE)
if (is_subset) {
  print("这些坐标是unique_in_ACE的子集")
} else {
  print("这些坐标不是unique_in_ACE的子集")
}

expression_data$gene <-  NA
expression_data$gene[row_indices] <- gene_name$hgnc_symbol
write.csv(expression_data, "data/probe_gene_1076_in_10707.csv")
