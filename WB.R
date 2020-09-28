library(ggplot2)
library(data.table)
theme_set(theme_bw()+ theme(axis.title = element_text(size = 12),
                            axis.text = element_text(size = 10),
                            plot.title = element_text(size = 14, hjust = 0.5),
                            legend.position="none"))



imr90 <- read.csv("data/WB_IMR90.csv", stringsAsFactors = FALSE)
imr90[, 4:6] <- imr90[, 4:6] / imr90[, 3]
aov_pold4 <- aov(POLD4~conc, data = imr90)
summary(aov_pold4)
aov_xpa <- aov(XPA~conc, data = imr90)
summary(aov_xpa)
aov_xpc <- aov(XPC~conc, data = imr90)
summary(aov_xpc)


kruskal.test(POLD4~conc, data = imr90)
kruskal.test(XPA~conc, data = imr90)
kruskal.test(XPC~conc, data = imr90)

imr90 <- setDT(imr90)
imr90_sta <- imr90[, .(XPA_sum = mean(XPA), XPA_se = sd(XPA)/sqrt(.N),
                       XPC_sum = mean(XPC), XPC_se = sd(XPC)/sqrt(.N),
                       POLD4_sum = mean(POLD4), POLD4_se = sd(POLD4)/sqrt(.N)),
                   by = conc]
imr90_sta$conc <- factor(imr90_sta$conc, levels = c(0, 1, 2.5, 5, 10, 25))

ggsave("WB_IMR90_XPA.png", 
       ggplot(imr90_sta) + geom_errorbar(aes(x = conc, ymax = XPA_sum+XPA_se, ymin = XPA_sum-0.1), size = 0.7, width = 0.3) +
         geom_bar(aes(x = conc, y = XPA_sum),stat = "identity", width = 0.6) +
         scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
         labs(title = "XPA expression in IMR-90", x = "Benzo[a]pyrene(μM)", y = "XPA"),
       height = 10, width = 10, units = "cm", dpi = 600
)

ggsave("WB_IMR90_XPC.png", 
       ggplot(imr90_sta) + geom_errorbar(aes(x = conc, ymax = XPC_sum+XPC_se, ymin = XPC_sum-0.1), size = 0.7, width = 0.3) +
         geom_bar(aes(x = conc, y = XPC_sum),stat = "identity", width = 0.6) +
         scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
         labs(title = "XPC expression in IMR-90", x = "Benzo[a]pyrene(μM)", y = "XPC"),
       height = 10, width = 10, units = "cm", dpi = 600
)

ggsave("WB_IMR90_POLD4.png", 
       ggplot(imr90_sta) + geom_errorbar(aes(x = conc, ymax = POLD4_sum+XPC_se, ymin = POLD4_sum-0.1), size = 0.7, width = 0.3) +
         geom_bar(aes(x = conc, y = POLD4_sum),stat = "identity", width = 0.6) +
         scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
         labs(title = "POLD4 expression in IMR-90", x = "Benzo[a]pyrene(μM)", y = "POLD4"),
       height = 10, width = 10, units = "cm", dpi = 600
)





beas <- read.csv("data/WB_BEAS2B.csv")
beas[, 4:6] <- beas[, 4:6] / beas[, 3]


kruskal.test(XPA~conc, data = beas)
kruskal.test(XPC~conc, data = beas)
kruskal.test(POLD4~conc, data = beas)



beas <- setDT(beas)
beas_sta <- beas[, .(XPA_sum = mean(XPA), XPA_se = sd(XPA)/sqrt(.N),
                     XPC_sum = mean(XPC), XPC_se = sd(XPC)/sqrt(.N),
                     POLD4_sum = mean(POLD4), POLD4_se = sd(POLD4)/sqrt(.N)),
                   by = conc]
beas_sta$conc <- factor(beas_sta$conc, levels = c(0, 1, 2.5, 5, 10, 25))


beas_t <- beas[conc == 0.0 | conc == 10.0, ]
t.test(XPA~conc, data = beas_t, var.equal = TRUE)
t.test(XPC~conc, data = beas_t, var.equal = TRUE)
t.test(POLD4~conc, data = beas_t, var.equal = TRUE)


ggsave("WB_BEAS2B_XPA.png", 
       ggplot(beas_sta) + geom_errorbar(aes(x = conc, ymax = XPA_sum+XPA_se, ymin = XPA_sum-0.1), size = 0.7, width = 0.3) +
               geom_bar(aes(x = conc, y = XPA_sum),stat = "identity", width = 0.6) +
               scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5)) +
               labs(title = "XPA expression in BEAS-2B", x = "Benzo[a]pyrene(μM)", y = "XPA"),
       height = 10, width = 10, units = "cm", dpi = 600
)

ggsave("WB_BEAS2B_XPC.png", 
       ggplot(beas_sta) + geom_errorbar(aes(x = conc, ymax = XPC_sum+XPC_se, ymin = XPC_sum-0.1), size = 0.7, width = 0.3) +
               geom_bar(aes(x = conc, y = XPC_sum),stat = "identity", width = 0.6) +
               scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5)) +
               labs(title = "XPC expression in BEAS-2B", x = "Benzo[a]pyrene(μM)", y = "XPC"),
       height = 10, width = 10, units = "cm", dpi = 600
)

ggsave("WB_BEAS2B_POLD4.png", 
       ggplot(beas_sta) + geom_errorbar(aes(x = conc, ymax = POLD4_sum+XPC_se, ymin = POLD4_sum-0.1), size = 0.7, width = 0.3) +
               geom_bar(aes(x = conc, y = POLD4_sum),stat = "identity", width = 0.6) +
               scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5)) +
               labs(title = "POLD4 expression in BEAS-2B", x = "Benzo[a]pyrene(μM)", y = "POLD4"),
       height = 10, width = 10, units = "cm", dpi = 600
)






mg132 <- read.csv("data/WB_BEAS2B_MG132.csv", stringsAsFactors = FALSE)
mg132[, 3:4] <- mg132[, 3:4] / mg132[, 2]

wilcox.test(XPA~conc, data = mg132)
wilcox.test(XPC~conc, data = mg132)

t.test(XPA~conc, data = mg132, var.equal = TRUE)
t.test(XPC~conc, data = mg132, var.equal = TRUE)

mg132 <- setDT(mg132)
mg132_sta <- mg132[, .(XPA_sum = mean(XPA), XPA_se = sd(XPA)/sqrt(.N),
                       XPC_sum = mean(XPC), XPC_se = sd(XPC)/sqrt(.N)),
                   by = conc]

ggsave("WB_BEAS2B_MG132_XPA.png", 
       ggplot(mg132_sta) + geom_errorbar(aes(x = conc, ymax = XPA_sum+XPA_se, ymin = XPA_sum-0.1), size = 0.7, width = 0.3) +
               geom_bar(aes(x = conc, y = XPA_sum),stat = "identity", width = 0.6) +
               scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5)) +
               labs(title = "XPA expression after MG132 treatment", x = "treat", y = "XPA"),
       height = 10, width = 10, units = "cm", dpi = 600
)

ggsave("WB_BEAS2B_MG132_XPC.png", 
       ggplot(mg132_sta) + geom_errorbar(aes(x = conc, ymax = XPC_sum+XPC_se, ymin = XPC_sum-0.1), size = 0.7, width = 0.3) +
               geom_bar(aes(x = conc, y = XPC_sum),stat = "identity", width = 0.6) +
               scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5)) +
               labs(title = "XPC expression after MG132 treatment", x = "treat", y = "XPC"),
       height = 10, width = 10, units = "cm", dpi = 600
)