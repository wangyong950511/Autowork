# 加载所需包
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(writexl)
library(ragg)
library(openxlsx)


rm(list = ls())



####需要手动修改的数据####
# 定义文件路径
filepath_input <- '/Users/wangyong/Desktop/2024-04-07_110658.xls'
# 指定顺序 
sample_order <- c("0")





####实际分析过程####
# 设置工作目录
# setwd("/home/drwang/R/Autoanalyse/RT-qPCR/PCR_data")
# 读取数据，跳过结果文件中前47行非数据部分
data <-
  read_excel(filepath_input,sheet = "Results",range = "A48:Q144",col_names=TRUE)


## 集合目的参数
# 定义结构文件路径
Sample_Gene_input <- '/Users/wangyong/Documents/实验攻略/RT-qPCR工具包/通用加样表.xlsx'
# 读取数据
Sample_data <-
  read_excel(Sample_Gene_input,sheet = 1,skip=1,col_names=FALSE,range = "A2:L9")
Gene_data <-
  read_excel(Sample_Gene_input,sheet = 1,skip=1,col_names=FALSE,range = "M2:X9")
# 将数据框从宽格式转换为长格式，转换成一列
Gene_long <- Gene_data %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Gene")
Sample_long <- Sample_data %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Sample")
# 替换到目标表格数据
data <- data %>%
  mutate(`Target Name` = Gene_long$Gene,`Sample Name` = Sample_long$Sample)
#参考样本设定
refsample <- na.omit(Sample_long$Sample)[1]
number<-"" 



# 标准化数据
colnames(data)[colnames(data) == "Sample Name"] <- "Sample"
colnames(data)[colnames(data) == "Target Name"] <- "Detector"
colnames(data)[colnames(data) == "CT"] <- "Ct"

#其中，"column"是要替换值的列名，"指定值"是要替换的值，"替换值"是要替换成的新值
data$Detector <- sub("Actin", "actin", data$Detector)
data$Detector <- sub("ACTIN", "actin", data$Detector)
data$Detector <- sub("ACTB", "actin", data$Detector)

# 剔除task为空白的行
data <- data %>% filter(Task != "")%>%
  group_by(Sample, Detector) %>%
  mutate(replicate = row_number())

# 按样本和基因分组，并计算平均值
data$Ct<-as.numeric(data$Ct)

#计算actin平均值
actin_mean <- data %>%
  filter(Detector == "actin") %>%
  group_by(Sample) %>%
  summarise(actin_mean = mean(Ct, na.rm = TRUE)) %>%
  ungroup()
# 找到所有样本和基因的组合
samples <- unique(data$Sample)
detectors <- setdiff(unique(data$Detector), "actin")

## 检测给定初始值是否正确 
# 判定 sample_order 和 samples 的元素是否不同（不考虑顺序）
if (!setequal(sample_order, samples)) {
  sample_order <- samples
}


# 初始化结果矩阵
results_matrix <- matrix(NA, nrow = length(samples) * length(detectors), ncol = 5)
dimnames(results_matrix) <- list(NULL, c("Sample", "Detector", "2^-ddCt_1", "2^-ddCt_2", "2^-ddCt_3"))

# 计数器，用来记录结果矩阵中已填入的行数
counter <- 0

# 记录计算失败的样品和基因组合
failed_combinations <- data.frame(Sample = character(), Detector = character(), stringsAsFactors = FALSE)

#进行以下步骤：
#按照样品和基因分组，计算 actin 的平均 Ct 值。
#对于每个样品和基因组合，找到对应的三个 gene Ct 值，计算三个 delta Ct 值。
#对于每个样品和基因组合，找到对应的三个 gene ΔCt 值和参考样品的 ΔCt值的平均值，计算三个 ΔΔCt 值。
#对于每个样品和基因组合，计算三个 2^-ΔΔCt 值。

# 遍历所有样品和基因的组合
for (sample in samples) {
  for (detector in detectors) {
    #test factor 
    #sample<-c('')
    #detector<-c('')
    # 找到该组合的三个 gene Ct 值
    ct_gene <- data %>%
      filter(Sample == sample & Detector == detector) %>%
      pull(Ct)
    # 如果有全部是 NA，则跳过该组合
    if (all(is.na(ct_gene))) {
      failed_combinations <- rbind(failed_combinations, data.frame(Sample = sample, Detector = detector, stringsAsFactors = FALSE))
      next
    }
    # 计算三个 delta Ct 值
    delta_ct <- ct_gene - actin_mean$actin_mean[actin_mean$Sample == sample]
    # 找到参考样品的 Ct 值
    ct_untreated <- data %>%
      filter(Sample == refsample & Detector == detector) %>%
      pull(Ct)
    # 如果参考样品的 Ct 值不存在，则跳过该组合
    if (all(is.nan(as.numeric(ct_untreated)))) {
      failed_combinations <- rbind(failed_combinations, data.frame(Sample = sample, Detector = detector, stringsAsFactors = FALSE))
      next
    }
    # 找到参考样品的 actin 平均 Ct 值
    actin_untreated_mean <- actin_mean$actin_mean[actin_mean$Sample == refsample]
    # 计算三个 ΔΔCt 值
    delta_delta_ct <- delta_ct - mean(ct_untreated - actin_untreated_mean, na.rm = TRUE)
    # 计算三个 2^-ΔΔCt 值
    two_to_minus_delta_delta_ct <- 2^-delta_delta_ct
    # 将结果填入结果矩阵
    counter <- counter + 1
    results_matrix[counter, "Sample"] <- sample
    results_matrix[counter, "Detector"] <- detector
    results_matrix[counter, "2^-ddCt_1"] <- two_to_minus_delta_delta_ct[1]
    results_matrix[counter, "2^-ddCt_2"] <- two_to_minus_delta_delta_ct[2]
    results_matrix[counter, "2^-ddCt_3"] <- two_to_minus_delta_delta_ct[3]
  }
}

# 去除 NA 行
results <- as.data.frame(results_matrix[rowSums(is.na(results_matrix)) < 3, ])
#合并原始Ct和ddct
# 将"2^-ddCt_1""2^-ddCt_2""2^-ddCt_3"这三列转换为一列
results_long <- results %>%
  pivot_longer(cols = starts_with("2^-ddCt"),
               names_to = "replicate",
               values_to =  "2^-ddCt") %>%
  group_by(Sample, Detector) %>%
  mutate(replicate = row_number())
# 转换数据类型
results_long$Sample <- coalesce(as.character(results_long$Sample), "")
results_long$Detector <- as.character(results_long$Detector)
# 按照sample和detector列进行匹配
merged_data <- left_join(data, results_long, by = c("Sample", "Detector", "replicate")) %>%
  arrange(Sample, Detector, replicate)
# 选择需要的列
final_data <- merged_data[merged_data$Detector != "actin",] %>% dplyr::select(Sample, Detector, Ct, `2^-ddCt`)%>%
  arrange(Sample, Detector)
filename1 <- paste0("result-", number, ".csv")  # 拼接文件名
write.csv(final_data, file = filename1, row.names = FALSE)
# 输出计算失败的样品和基因组合
failed_combinations <- unique(failed_combinations)
filename2 <- paste0("failed_combinations-", number, ".csv")  # 拼接文件名
write.csv(failed_combinations, file = filename2, row.names = FALSE)
# 删除 results_long 中某一列为空值的行
results_long <- results_long %>%
  filter(!is.na(`2^-ddCt`))

# library(ggplot2)
# # 将Sample列的样本名转换为数值并排序
# results_long$Sample <- factor(results_long$Sample, levels = unique(results_long$Sample))
#
# # 使用ggplot2绘制RT-qPCR图形
# ggplot(results_long, aes(x = Sample,
#                          y = `2^-ddCt`,
#                          fill = as.factor(Sample))) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_wrap(~Detector, scales = "free") +  # 每个基因一个子图
#   labs(title = "RT-qPCR Analysis",
#        x = "Sample",
#        y = "2^-ddCt Value",
#        fill = "Sample") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 旋转x轴标签


# 创建一个新的数据框，其中包含每个样本的平均值和标准差
summary_data <- results_long %>%
  group_by(Sample, Detector) %>%
  mutate(mean_value = mean(as.numeric(`2^-ddCt`), na.rm = TRUE),
         sd_value = sd(as.numeric(`2^-ddCt`), na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(Sample, Detector, .keep_all = TRUE) %>%
  filter(!is.na(mean_value))

# 确保 Sample 列按指定顺序排列
summary_data$Sample <- factor(summary_data$Sample, levels = sample_order)

# 使用ggplot2绘制RT-qPCR基因表达量图（柱状图）
p <- ggplot(summary_data, aes(x = Sample,
                              y = mean_value,
                              fill = Detector)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7,
           color = "black") +
  geom_errorbar(aes(ymin = mean_value - sd_value,
                    ymax = mean_value + sd_value),
                width = 0.25,
                position = position_dodge(width = 0.8),
                color = "black") +
  labs(title = "RT-qPCR Gene Expression Analysis",
       x = number,
       y = "2^-ddCt Value",
       fill = "Detector") +
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 14,
                                   color = "black"),
        axis.text.y = element_text(size = 14,
                                   color = "black"),
        axis.title = element_text(size = 16,
                                  color = "black"),
        legend.title = element_text(size = 14,
                                    color = "black"),
        legend.text = element_text(size = 12,
                                   color = "black"),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 1)) +  # 使用linewidth替代size
  scale_fill_brewer(palette="Set2") +  # 使用Brewer颜色调色板
  facet_wrap(~ Detector, scales = "free_y")

# 打印图表
print(p)

####保存文档####
# 定义输出文件名
output_filename <- 'output.xlsx'
plot_name <- 'Qpcr.svg'
output2_filename <- '加样表.xlsx'
# 调整格式   
results_long$replicate <- NULL
results_long <- results_long %>%
  mutate(Sample = factor(Sample, levels = sample_order)) %>%
  arrange(Detector, Sample)
transposed_results <- t(results_long)
# 使用write_xlsx函数保存数据框到指定路径
# 提取文件夹路径
folder_path <- dirname(filepath_input)
# 创建输出文件路径
filepath_output <- file.path(folder_path, output_filename)
write_xlsx(as.data.frame(transposed_results), filepath_output)
# 保存图片
# 导出 SVG 文件
file_path <- file.path(folder_path, plot_name)
ggsave(file_path, plot = p, device = "svg", width = 10, height = 5, units = "in")
# 保存加样表格图
wb <- loadWorkbook(Sample_Gene_input) 
filepath_output <- file.path(folder_path, output2_filename)
saveWorkbook(wb, filepath_output, overwrite = TRUE)
