library(data.table)
library(ggplot2)
snd_data <- fread("data/qPCR_data.csv",stringsAsFactors = FALSE)
fst_data <- fread("data/1st_data.csv",stringsAsFactors = FALSE)
fst_data$type <- factor(fst_data$type, labels = c("No", "Mild", "Moderate", "Heavy"),
                        levels = c("No Smoking", "Mild Smoking", "Moderate Smoking", "Heavy Smoking"))
group <- read.csv("data/group.csv",stringsAsFactors = FALSE)
group$type <- factor(group$type, labels = c("No", "Mild", "Moderate", "Heavy"),
                     levels = c("No Smoking", "Mild Smoking", "Moderate Smoking", "Heavy Smoking"))

#1st result####
ct_1st <- cbind(fst_data[,1:2], fst_data[,4:6]-fst_data$GAPDH)
ct_1st <- melt(ct_1st, id.var = c("Sample.Name", "type"), variable.name = "gene")

#anova
aov_1st_XPA <- aov(value ~ type, data = ct_1st[gene == "XPA",])
summary(aov_1st_XPA)
aov_1st_XPC <- aov(value ~ type, data = ct_1st[gene == "XPC",])
summary(aov_1st_XPC)
TukeyHSD(aov_1st_XPC)# Mild Smoking-Heavy Smoking 
aov_1st_POLD4 <- aov(value ~ type, data = ct_1st[gene == "POLD4",])
summary(aov_1st_POLD4)

#delta ct boxplot
ggplot(ct_1st)+aes(type, value)+geom_boxplot()+facet_grid(.~gene)+
  theme_bw()+theme(text = element_text(size= 25),axis.text.x = element_text(angle = 45, hjust = 0.9))+labs(x = "Smoking Group", y = "¦¤Ct")

#delta ct barplot
ct.bar <- ct_1st[,.(mean = mean(value), se = sd(value)/sqrt(.N)), by = c("gene","type")]

ggsave(filename = "ct1_barplot.jpeg",
       ggplot(ct.bar)+aes(type, mean)+geom_bar(stat = "identity", fill = "#222222", color = 1, width = 0.8)+
         geom_errorbar(aes(ymax = mean+se, ymin = mean-0.01), size = 0.7, width = 0.5)+facet_grid(.~gene)+
         theme_bw()+theme(text = element_text(size= 25),axis.text.x = element_text(angle = 45, hjust = 0.9))+
         scale_y_continuous(limits = c(0,9), expand = c(0,0))+labs(x = "Smoking Group", y = "¦¤Ct"),
       width = 22, height = 16, dpi = 400, units = "cm", device = "jpeg")

#relative expression barplot
re_1st <- ct_1st[,mean(value, na.rm = TRUE), by = c("type", "gene"),]
re_1st <- dcast(re_1st, gene~type)
re_1st <- cbind(re_1st[,1], re_1st[, 3:5] - re_1st$No)
re_1st <- melt(re_1st, id.vars = "gene", variable.name = "type")
re_1st$value <- 2^-(re_1st$value)#2^-¦¤¦¤CT relative expression

ggsave(filename = "re1_barplot.jpeg",
       ggplot(re_1st)+aes(type, value-1)+geom_bar(stat = "identity", fill = "#222222", width = 0.8)+geom_hline(yintercept = 0)+facet_grid(.~gene)+
         scale_y_continuous(labels = function(y) y + 1,limits = c(-1,1.4),expand = c(0,0))+labs(x = "Smoking Group", y = "relative expression")+
         theme_bw()+theme(text = element_text(size= 25),axis.text.x = element_text(angle = 45, hjust = 0.9)),
width = 22, height = 16, dpi = 400, units = "cm", device = "jpeg")

###all smoke group
bind_1st <- ct_1st
bind_1st[type != "No",]$type <- "Smoke"

t.test(value~type, data = bind_1st[gene == "XPA" ], var.equal = TRUE)
t.test(value~type, data = bind_1st[gene == "XPC" ], var.equal = TRUE)
t.test(value~type, data = bind_1st[gene == "POLD4" ], var.equal = TRUE)

bind1_bar <- bind_1st[,.(mean = mean(value), se = sd(value)/sqrt(.N)), by = c("gene","type")]
ggsave(filename = "bind1_barplot.jpeg",
       ggplot(bind1_bar)+aes(type, mean)+geom_bar(stat = "identity", fill = "#222222", color = 1, width = 0.8)+
         geom_errorbar(aes(ymax = mean+se, ymin = mean-0.01), size = 0.7, width = 0.5)+facet_grid(.~gene)+
         theme_bw()+theme(text = element_text(size= 25),axis.text.x = element_text(angle = 45, hjust = 0.9))+
         scale_y_continuous(limits = c(0,9), expand = c(0,0))+labs(x = "Smoking Group", y = "¦¤Ct"),
       width = 22, height = 16, dpi = 400, units = "cm", device = "jpeg")

bind1_re <- dcast(bind1_bar[,1:3], gene~type)
bind1_re$re <- bind1_re$Smoke-bind1_re$No
bind1_re$re <- 2^-(bind1_re$re)

ggsave(filename = "bind1_re.jpeg",
       ggplot(bind1_re)+aes(gene, re-1)+geom_bar(stat = "identity", fill = "#222222", width = 0.8)+geom_hline(yintercept = 0)+
         scale_y_continuous(labels = function(y) y + 1,limits = c(-1,1.4),expand = c(0,0))+labs(x = "gene", y = "relative expression")+
         theme_bw()+theme(text = element_text(size= 25),axis.text.x = element_text(angle = 45, hjust = 0.9)),
       width = 22, height = 16, dpi = 400, units = "cm", device = "jpeg")

#2nd result####
ct_mean <- snd_data[, mean(C§ä), by = .(Sample.Name, Target.Name)]
ct_mean <- dcast(ct_mean, Sample.Name~Target.Name)
ct_mean <- cbind(ct_mean[,1], ct_mean[,3:5]-ct_mean$GAPDH)
ct_mean <- merge(group, ct_mean, by = c("Sample.Name"))
ct_mean <- melt(ct_mean,id.var = c("Sample.Name", "type"), variable.name = "gene")

#anova
aov_POLD1 <- aov(value ~ type, data = ct_mean[gene == "POLD1"])
anova(aov_POLD1)
aov_POLD2 <- aov(value ~ type, data = ct_mean[gene == "POLD2"])
anova(aov_POLD2)
TukeyHSD(aov_POLD2)
aov_POLD3 <- aov(value ~ type, data = ct_mean[gene == "POLD3"])
anova(aov_POLD3)

#delta ct boxplot
ggsave(filename = "ct_boxplot.jpeg",
  ggplot(ct_mean)+aes(type, value)+geom_boxplot()+facet_grid(.~gene)+theme_bw()+
    theme(text = element_text(size= 25),axis.text.x = element_text(angle = 45, hjust = 0.9))+
    scale_y_continuous(limits = c(0,12.4), expand = c(0,0))+labs(x = "Smoking Group", y = "¦¤Ct"),
  width = 22, height = 16, dpi = 400, units = "cm", device = "jpeg")

#delta ct barplot
ct.bar <- ct_mean[,.(mean = mean(value), se = sd(value)/sqrt(.N)), by = c("gene","type")]

ggsave(filename = "ct_barplot.jpeg",
       ggplot(ct.bar)+aes(type, mean)+geom_bar(stat = "identity", fill = "#222222", color = 1, width = 0.8)+
         geom_errorbar(aes(ymax = mean+se, ymin = mean-0.01), size = 0.7, width = 0.5)+facet_grid(.~gene)+
         theme_bw()+theme(text = element_text(size= 25),axis.text.x = element_text(angle = 45, hjust = 0.9))+
         scale_y_continuous(limits = c(0,12.4), expand = c(0,0))+labs(x = "Smoking Group", y = "¦¤Ct"),
       width = 22, height = 16, dpi = 400, units = "cm", device = "jpeg")

#relative expression barplot
re_2nd <- ct_mean[,mean(value, na.rm = TRUE), by = c("type", "gene"),]
re_2nd <- dcast(re_2nd, gene~type)
re_2nd <- cbind(re_2nd[,1], re_2nd[, 3:5] - re_2nd$No)
re_2nd <- melt(re_2nd, id.vars = "gene", variable.name = "type")
re_2nd$value <- 2^-(re_2nd$value)#2^-¦¤¦¤CT relative expression

ggsave(filename = "re_barplot.jpeg",
       ggplot(re_2nd)+aes(type, value-1)+geom_bar(stat = "identity", fill = "#222222", width = 0.8)+geom_hline(yintercept = 0)+facet_grid(.~gene)+
         scale_y_continuous(labels = function(y) y + 1, expand = c(0,0), limits = c(-1, 5.2))+labs(x = "Smoking Group", y = "relative expression")+
         theme_bw()+theme(text = element_text(size= 25),axis.text.x = element_text(angle = 45, hjust = 0.9)),
       width = 22, height = 16, dpi = 400, units = "cm", device = "jpeg")

###all smoke group
bind_2nd <- ct_mean
bind_2nd[type != "No",]$type <- "Smoke"

t.test(value~type, data = bind_2nd[gene == "POLD1"], var.equal = TRUE)
t.test(value~type, data = bind_2nd[gene == "POLD2"], var.equal = TRUE)
t.test(value~type, data = bind_2nd[gene == "POLD3"], var.equal = TRUE)

bind2_bar <- bind_2nd[,.(mean = mean(value), se = sd(value)/sqrt(.N)), by = c("gene","type")]
ggsave(filename = "bind2_barplot.jpeg",
       ggplot(bind2_bar)+aes(type, mean)+geom_bar(stat = "identity", fill = "#222222", color = 1, width = 0.8)+
         geom_errorbar(aes(ymax = mean+se, ymin = mean-0.01), size = 0.7, width = 0.5)+facet_grid(.~gene)+
         theme_bw()+theme(text = element_text(size= 25),axis.text.x = element_text(angle = 45, hjust = 0.9))+
         scale_y_continuous(limits = c(0,12.4), expand = c(0,0))+labs(x = "Smoking Group", y = "¦¤Ct"),
       width = 22, height = 16, dpi = 400, units = "cm", device = "jpeg")

bind2_re <- dcast(bind2_bar[,1:3], gene~type)
bind2_re$re <- bind2_re$Smoke-bind2_re$No
bind2_re$re <- 2^-(bind2_re$re)

ggsave(filename = "bind2_re.jpeg",
       ggplot(bind2_re)+aes(gene, re-1)+geom_bar(stat = "identity", fill = "#222222", width = 0.8)+geom_hline(yintercept = 0)+
         scale_y_continuous(labels = function(y) y + 1,limits = c(-1,5.2),expand = c(0,0))+labs(x = "gene", y = "relative expression")+
         theme_bw()+theme(text = element_text(size= 25),axis.text.x = element_text(angle = 45, hjust = 0.9)),
       width = 22, height = 16, dpi = 400, units = "cm", device = "jpeg")