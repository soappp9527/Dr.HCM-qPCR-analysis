library(data.table)
library(ggplot2)
theme_set(theme_bw()+ theme(text = element_text(size=20),legend.position="none"))

data <- fread("data/AS1.csv", stringsAsFactors = FALSE)
data <- data[`GAPDH CT Mean` < 30,]
data <- data[`正常/癌化` != "",]
data$delta <- data$`AS1 CT Mean`-data$`GAPDH CT Mean`
data[,.N,by = `正常/癌化`]
data[`正常/癌化`=="癌化",]$`正常/癌化` <- "Cancerous"
data[`正常/癌化`=="正常",]$`正常/癌化` <- "Normal"
data$`正常/癌化` <- factor(data$`正常/癌化`, levels = c("Normal", "Cancerous"))

t.test(delta~`正常/癌化`, data = data, var.equal = TRUE)

plt <- data[,.(mean = mean(delta), se = sd(delta)/sqrt(.N)), by = `正常/癌化`]
relative_expression <- 2**-(plt[2, 2]-plt[1, 2])


ggsave("AS1_barplot.png", 
  ggplot(plt)+geom_errorbar(aes(x = `正常/癌化`, ymax = mean + se, ymin = mean - 0.1), size = 0.7, width = 0.3)+
    geom_bar(aes(x = `正常/癌化`, y = mean),stat = "identity", width = 0.6)+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 9))+
    labs(x = "Prostate tissue", y = "ΔCt (ABHD11-AS1)"),
width = 16, height = 16 ,units = "cm", dpi = 600)  

ggsave("AS1_boxplot.png", 
  ggplot(data)+aes(x = `正常/癌化`, y = delta, fill = `正常/癌化`)+geom_boxplot(width = 0.5)+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 14), breaks = seq(from = 0, to = 15, by = 3))+
    scale_fill_manual(values=c("#87CDCF", "#FF6095"))+
    labs(x = "Prostate tissue", y = "ΔCt (ABHD11-AS1)"),
width = 16, height = 16 ,units = "cm", dpi = 600)  