# 加载所需包
library(readxl)
library(tidyr)
library(ggplot2)



####需要手动修改的数据####
# 定义文件路径
filepath_input <-
  '/Users/wangyong/Documents/实验数据分析/BCA/Method 1_20241124_194723.xlsx'
# 获得蛋白体积
Volum = 200
# 是否新做标准曲线
New_Standardcurve=TRUE
# 既往标准曲线系数
x_factor=28.192847
constant=-4.170041
# 样本名
# group_labels <- c("NC1","NC2","NC3","R1","R2","R3","T1","T2","T3","R+T1","R+T2","R+T3","R+T4")
# group_labels <- c("1","2","3","4","5","6","7","8","9")
# group_labels <- c("NC1","NC2","NC3","C2w_1","C2w_2","C2w_3","C4w_1","C4w_2","C4w_3")
group_labels <- c("NC-1","NC-2","NC-3","WD-1","WD-2","WD-3","IRT-1","IRT-2","IRT-3","M+I-1","M+I-2","M+I-3")


####标准曲线####
# 读取数据，跳过结果文件中前47行非数据部分
data <-
  read_excel(filepath_input,sheet = "Result sheet",col_names=TRUE,range = "B43:M51")
#计算标准曲线#
if (New_Standardcurve){
  x <- as.numeric(data[1, 1:8])
  y=c(0,0.05,0.1,0.2,0.4,0.6,0.8,1)
  # 拟合标准曲线（线性回归）
  fit <- lm(y ~ x)
  # 查看拟合结果
  summary_fit <- summary(fit)
  # 提取 R^2 值
  r_squared <- summary_fit$r.squared
  # 绘制标准曲线
  plot(x, y, main = "Standardcurve", xlab = "Absorbance", ylab = "C")
  abline(fit, col = "red")
  # 在图表上显示 R^2 值
  r_text <- paste("R² = ", round(r_squared, 4))  # 将 R² 值保留 4 位小数
  text(x = min(x), y = max(y), labels = r_text, pos = 4, col = "blue")  # 显示在左上角
  # 输出拟合方程的系数
  coefficients(fit)
  coef_values <- coef(fit)
  constant <- as.numeric(coef_values[1])
  x_factor <- as.numeric(coef_values[2])
  # 删除数据框的第一行
  data <- data[-1, ]
} else {
  # 如果不计算标准曲线，可以执行其他操作
  cat("不计算标准曲线\n")
}



####分析数据####
# 定义每组的列数
group_size <- 3
# 计算需要的组数
num_groups <- ncol(data) / group_size
# 创建一个列表，用于存储分组后的数据框
data_split <- vector("list", length = num_groups)
# 将数据框按照每3列一组分成4组
for (i in 1:num_groups) {
  cols <- ((i - 1) * group_size + 1):(i * group_size)
  data_split[[i]] <- data[, cols]
  colnames(data_split[[i]]) <- c("1", "2", "3")
}
# 组合分组后的数据框
combined_data <- rbind(data_split[[1]], data_split[[2]], data_split[[3]], data_split[[4]])
# 删除含有空值的行
data_clean <- combined_data[complete.cases(combined_data), ]
# 定义一个函数，根据标准曲线方程计算新值
compute_new_value <- function(x) {
  return((x_factor * x + constant)*10)
}
# 使用 apply() 函数应用计算函数到数据框每个元素
new_data <- as.data.frame(apply(data_clean, 2, compute_new_value))
# 计算每行的平均值和标准差
row_means <- apply(new_data, 1, mean)
row_sd <- apply(new_data, 1, sd)
# 计算每行的均值和标准差
row_means <- apply(new_data[, -1], 1, mean)
row_sd <- apply(new_data[, -1], 1, sd)
# 创建包含均值和标准差的新数据框
summary_data <- data.frame(
  new_data,
  Mean = row_means,
  SD = row_sd
)



####评估加样质量####
# 新建行号列
summary_data$Sample <- as.factor(1:nrow(summary_data))
# 输出所有SD大于Mean的组
bad_row <- which(summary_data$SD > summary_data$Mean * 0.2)
output_text1 <- if (length(bad_row) > 0) {
  paste(paste(bad_row, collapse = ","), "Need to be redone")
} else {
  "All samples meet the requirements,  "
}
# 使用ggplot2绘制柱状图和误差线
p <- ggplot(summary_data, aes(x = Sample, y = Mean)) +
  geom_bar(stat = "identity", fill = "blue", alpha = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, color = "black") +
  labs(title = "Mean and SD Bar Plot",
       x = "Sample",
       y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 打印图形
print(p)



####计算过程####
# 找到浓度最小值
min_value <- min(summary_data$Mean)
# 计算蛋白加样体积
summary_data$Protein_volum <- as.integer((Volum * min_value)/summary_data$Mean)
# 计算RIPA体积
summary_data$RIPA_volum <- Volum - summary_data$Protein_volum
# 导出新数组
out_data <- summary_data[, c("Mean", "RIPA_volum", "Protein_volum")]
# 添加行号作为组别
out_data$group <- as.factor(1:nrow(out_data))
# 将数据框转换为长格式
data_long <- gather(out_data, key = "variable", value = "value",Protein_volum,RIPA_volum)
data_long$variable <- factor(data_long$variable, levels = c("RIPA_volum", "Protein_volum"))
# 将 'group' 变量的顺序反转使样本顺序从上到下
data_long$group <- factor(data_long$group, levels = rev(levels(data_long$group)))
group_labels_rev <- rev(group_labels)
# 整理文本
output_text2 <- paste("Protein concentration=",round(min_value, digits = 1),"ug/ul,"," Protein volume=",Volum,"µL")
# 拼接标题文字
title_text <- paste("Add Rules -", output_text1, "\n", output_text2)
# 使用ggplot2绘制堆叠柱状图
p <- ggplot(data_long, aes(x = value, y = group, fill = variable)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5), size = 6, color = "black") +
  labs(title = title_text, x = "Volum", y = "Sample") +
  scale_fill_manual(values = c("Protein_volum" = "red", "RIPA_volum" = "yellow")) +
  scale_y_discrete(labels = group_labels_rev) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# 打印图形
print(p)
cat(output_text1,"\n",output_text2)
