#RXNs in steroid-----------------------------------------------
steroid <- read.table("../matlab/HSC metabolism/HSCmodel/forPaper/Reactions for subsystem steroid_metabolism.tsv",header = T,sep = '\t')
steroidout <- out[steroid$Reaction.ID,]

# 计算每个代谢物在早期和晚期病人中出现的次数
early_counts <- rowSums(data_matrix[, 1:10] == 1)
late_counts <- rowSums(data_matrix[, 11:20] == 1)

# 计算每个代谢物在两组中出现的比例
early_proportions <- early_counts / 10  # 早期病人样本数为10
late_proportions <- late_counts / 10    # 晚期病人样本数为10

# 将结果组合成数据框
metabolite_data <- data.frame(
  rxns = rownames(data_matrix),  # 代谢物编号
  Young_Proportion = early_proportions,
  Aging_Proportion = late_proportions,
  delta = late_proportions-early_proportions
)

# 输出结果
steroidrxns <- metabolite_data[metabolite_data$rxns %in% steroid$Reaction.ID,]
newsteroidrxns <- merge(steroidrxns,steroid,by.x='rxns',by.y = 'Reaction.ID')
library(tidyr)
newsteroidrxns <- newsteroidrxns[order(newsteroidrxns$delta,decreasing = T),]
# 将数据从宽格式转换为长格式
df <- melt(newsteroidrxns[,c(2,3,5)], id.vars = "Equation", variable.name = "Group", value.name = "Proportion")
df$value <- ifelse(df$Group=="Young_Proportion",-df$Proportion,df$Proportion)

library(ggplot2)

#颜色
col <- c("#ffc8bc","#ff715e")
top10 <- newsteroidrxns[1:10,]$Equation
plotdf <- df[df$Equation %in% top10,]
plotdf$Equation <- factor(plotdf$Equation ,levels = top10)
#绘图
library(dplyr)

# 使用group_by()函数按照id进行分组，并使用summarize()函数对每个分组中的value进行求和
summarized_data <- plotdf %>% 
  group_by(Equation) %>% 
  summarize(total_value = sum(value, na.rm = TRUE))
plotdfs <- merge(plotdf,summarized_data,by='Equation')
plotdfs$Equation <- factor(plotdfs$Equation,levels = levels(plotdfs$Equation)[c(2:3,1,4:10)])

ggplot(plotdfs,aes(x=Equation,y=value,fill=Group,group=Group))+
  #柱状图
  geom_col(width = 0.7)+
  geom_point(aes(y=total_value),col="black",size=2.5)+
  geom_line(aes(x=Equation, y=total_value),col="black",size =0.8)+
  #坐标转换
  coord_flip()+
  #自定义颜色
  scale_fill_manual(values = col)+
  geom_text(aes(label = Proportion),size = 4, hjust=1.2) +
  #主题设置
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black",size = 14,face = "italic"),
        axis.title = element_text(color = "black",size = 16),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(color = "black",size = 15))+
  #标题设置
  labs(x=NULL,y="rxns Frequency")+
  #y轴范围设置
  scale_y_continuous(breaks = seq(-1, 1, 0.5), 
                     labels = as.character(abs(seq(-1, 1, 0.5))),
                     limits = c(-1, 1))


#Dot plot---------------------------------------------------------------------
filtered_df <- ego@result[grepl("^regulation of.*metabolic process$", ego@result$Description), ]
plotdf <- filtered_df %>% filter(pvalue<0.05,Count>=5)
forward <- as.numeric(sub("/\\d+$", "", plotdf$GeneRatio))
backward <- as.numeric(sub("^\\d+/", "", plotdf$GeneRatio))
## add GeneRatio
plotdf$GeneRatio = forward/backward
rownames(plotdf) <- plotdf$Description

library(tidyverse)
library(viridis)
plotdf%>% 
  ## p value
  arrange(-log10(pvalue)) %>% 
  ggplot(aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  geom_point(aes(color=-log10(pvalue), size = Count)) +
  scale_color_viridis(discrete = F, option = "F")+
  scale_x_continuous(position = "top") + 
  scale_size_continuous(range=c(7,10))+
  labs(y=NULL) +
  ggtitle("")+
  ## theme
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = 15, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = 15, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = 15),
        axis.title.y = element_text(angle=90))



