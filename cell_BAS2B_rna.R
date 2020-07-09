library(data.table)
library(ggplot2)
theme_set(theme_bw()+ theme(text = element_text(size=20),legend.position="none"))

data <- fread("data/BEAS-2B-RNA.csv", stringsAsFactors = FALSE)
data$conc <- factor(data$conc, levels = c(0.0, 1, 2.5, 5, 10, 25))
xpa <- aov(XPA~conc, data = data)
summary(xpa)
xpc <- aov(XPC~conc, data = data)
summary(xpc)
pold4 <- aov(POLD4~conc, data = data)
summary(pold4)


#plot
sta <- data[,.(XPA_sum = mean(XPA), XPA_se = sd(XPA)/sqrt(.N),
               XPC_sum = mean(XPC), XPC_se = sd(XPC)/sqrt(.N),
               POLD4_sum = mean(POLD4), POLD4_se = sd(POLD4)/sqrt(.N)),
            by = conc]
ggsave("XPA.png", 
  ggplot(sta)+geom_errorbar(aes(x = conc, ymax = XPA_sum+XPA_se, ymin = XPA_sum-0.1), size = 0.7, width = 0.3)+
    geom_bar(aes(x = conc, y = XPA_sum),stat = "identity", width = 0.6)+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 17))+
    labs(x = "Benzoapyrene", y = "Δct")
)

ggsave("XPC.png", 
       ggplot(sta)+geom_errorbar(aes(x = conc, ymax = XPC_sum+XPC_se, ymin = XPC_sum-0.1), size = 0.7, width = 0.3)+
         geom_bar(aes(x = conc, y = XPC_sum),stat = "identity", width = 0.6)+
         scale_y_continuous(expand = c(0, 0), limits = c(0, 17))+
         labs(x = "Benzoapyrene", y = "Δct")
)

ggsave("POLD4.png", 
       ggplot(sta)+geom_errorbar(aes(x = conc, ymax = POLD4_sum+XPC_se, ymin = POLD4_sum-0.1), size = 0.7, width = 0.3)+
         geom_bar(aes(x = conc, y = POLD4_sum),stat = "identity", width = 0.6)+
         scale_y_continuous(expand = c(0, 0), limits = c(0, 17))+
         labs(x = "Benzoapyrene", y = "Δct")
)