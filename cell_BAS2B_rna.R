library(data.table)
library(ggplot2)
theme_set(theme_bw()+ theme(axis.title = element_text(size = 12),
                            axis.text = element_text(size = 10),
                            plot.title = element_text(size = 14, hjust = 0.5),
                            legend.position="none"))

#BEAS-2B####

#import ABI-7500 data table
bas1 <- fread("data/20200616-BEAS-2B_QPCR_GYP_data.csv", stringsAsFactors = FALSE, skip = 7)
bas2 <- fread("data/20200629-BEAS-2B-GYP_data.csv", stringsAsFactors = FALSE, skip = 7)
bas3 <- fread("data/20200630-BEAS-2B-GYP_data.csv", stringsAsFactors = FALSE, skip = 7)

cols_need <- c("Sample Name", "Target Name", "Cт")
bas1 <- bas1[1:(nrow(bas1)-5), ..cols_need][Cт != "Undetermined",]
bas1$batch <- 1
bas2 <- bas2[1:(nrow(bas2)-5), ..cols_need][Cт != "Undetermined",][, c("batch", "Sample Name") := tstrsplit(`Sample Name`, "-")]
bas3 <- bas3[1:(nrow(bas3)-5), ..cols_need][Cт != "Undetermined",][, c("batch", "Sample Name") := tstrsplit(`Sample Name`, "-")]

bas_b <- rbind(bas1, bas2, bas3)
bas_b <- bas_b[!is.na(`Sample Name`),][, c("conc", "target", "ct") := .SD, .SDcols = cols_need]
bas_b$conc <- factor(bas_b$conc, levels = c(0.0, 1, 2.5, 5, 10, 25))
bas_b$ct <- as.numeric(bas_b$ct)

out_ct <- c(94, 106, 154, 166, 20, 161, 72, 38, 96, 40, 41, 132, 144, 46, 53, 77)
bas_b <- bas_b[-out_ct,]
bas_c <- bas_b[, .(ct_mean = mean(ct)), by = c("batch", "conc", "target")]
bas_c <-  dcast(bas_c, batch + conc ~ target, value.var = "ct_mean")
bas_c[, 4:6] <- bas_c[, 4:6] - bas_c$GAPDH

#write.csv(bas_c, "BEAS-2B_RNA.csv", row.names = FALSE)
kruskal.test(XPA ~ conc, data = bas_c)
kruskal.test(XPC ~ conc, data = bas_c)
kruskal.test(POLD4 ~ conc, data = bas_c)



ggplot(bas_c) + aes(x = conc, y = XPA) + geom_boxplot()
ggplot(bas_c) + aes(x = conc, y = XPC) + geom_boxplot()
ggplot(bas_c) + aes(x = conc, y = POLD4) + geom_boxplot()



bas_sta <- bas_c[, .(XPA_sum = mean(XPA), XPA_se = sd(XPA)/sqrt(.N),
                     XPC_sum = mean(XPC), XPC_se = sd(XPC)/sqrt(.N),
                     POLD4_sum = mean(POLD4), POLD4_se = sd(POLD4)/sqrt(.N)),
                 by = conc]

ggsave("BAS_XPA.png", 
       ggplot(bas_sta) + geom_errorbar(aes(x = conc, ymax = XPA_sum+XPA_se, ymin = XPA_sum-0.1), size = 0.7, width = 0.3) +
         geom_bar(aes(x = conc, y = XPA_sum),stat = "identity", width = 0.6) +
         scale_y_continuous(expand = c(0, 0), limits = c(0, 17)) +
         labs(title = "XPA mRNA expression in BEAS-2B", x = "Benzo[a]pyrene(μM)", y = "ΔCt(XPA)"),
       height = 10, width = 10, units = "cm", dpi = 600
)

ggsave("BAS_XPC.png", 
       ggplot(bas_sta) + geom_errorbar(aes(x = conc, ymax = XPC_sum+XPC_se, ymin = XPC_sum-0.1), size = 0.7, width = 0.3) +
         geom_bar(aes(x = conc, y = XPC_sum),stat = "identity", width = 0.6) +
         scale_y_continuous(expand = c(0, 0), limits = c(0, 17)) +
         labs(title = "XPC mRNA expression in BEAS-2B", x = "Benzo[a]pyrene(μM)", y = "ΔCt(XPC)"),
       height = 10, width = 10, units = "cm", dpi = 600
)

ggsave("BAS_POLD4.png", 
       ggplot(bas_sta) + geom_errorbar(aes(x = conc, ymax = POLD4_sum+XPC_se, ymin = POLD4_sum-0.1), size = 0.7, width = 0.3) +
         geom_bar(aes(x = conc, y = POLD4_sum),stat = "identity", width = 0.6) +
         scale_y_continuous(expand = c(0, 0), limits = c(0, 17)) +
         labs(title = "POLD4 mRNA expression in BEAS-2B", x = "Benzo[a]pyrene(μM)", y = "ΔCt(POLD4)"),
       height = 10, width = 10, units = "cm", dpi = 600
)



#IMR-90####
imr <- fread("data/IMR90-RNA-IMR90-20200723.csv", stringsAsFactors = FALSE)
imr$conc <- factor(imr$浓度, levels = c(0.0, 1, 2.5, 5, 10, 25))
imr <- imr[批次!=1]

write.csv(imr, "IMR90_RNA.csv", row.names = FALSE)

#stat
kruskal.test(`XPA △Ct` ~ conc, data = imr)
kruskal.test(`XPC △Ct` ~ conc, data = imr)
kruskal.test(`POLD4 △Ct` ~ conc, data = imr)

posthoc.kruskal.nemenyi.test(`XPA △Ct` ~ conc, data = imr, method="Tukey")
imr_xpa
TukeyHSD(imr_xpa)
imr_xpc <- aov(`XPC △Ct` ~ conc, data = imr)
summary(imr_xpc)
TukeyHSD(imr_xpc)
imr_pold4 <- aov(`POLD4 △Ct` ~ conc, data = imr)
summary(imr_pold4)
TukeyHSD(imr_pold4)

#plot
imr_sta <- imr[, .(XPA_sum = mean(`XPA △Ct`), XPA_se = sd(`XPA △Ct`)/sqrt(.N),
                   XPC_sum = mean(`XPC △Ct`), XPC_se = sd(`XPC △Ct`)/sqrt(.N),
                   POLD4_sum = mean(`POLD4 △Ct`), POLD4_se = sd(`POLD4 △Ct`)/sqrt(.N)),
               by = conc]

ggsave("IMR_XPA.png", 
       ggplot(imr_sta) + geom_errorbar(aes(x = conc, ymax = XPA_sum+XPA_se, ymin = XPA_sum-0.1), size = 0.7, width = 0.3) +
         geom_bar(aes(x = conc, y = XPA_sum),stat = "identity", width = 0.6) +
         scale_y_continuous(expand = c(0, 0), limits = c(0, 17)) +
         labs(title = "XPA mRNA expression in IMR-90", x = "Benzo[a]pyrene(μM)", y = "ΔCt(XPA)"),
       height = 10, width = 10, units = "cm", dpi = 600
)

ggsave("IMR_XPC.png", 
       ggplot(imr_sta) + geom_errorbar(aes(x = conc, ymax = XPC_sum+XPC_se, ymin = XPC_sum-0.1), size = 0.7, width = 0.3) +
         geom_bar(aes(x = conc, y = XPC_sum),stat = "identity", width = 0.6) +
         scale_y_continuous(expand = c(0, 0), limits = c(0, 17)) +
         labs(title = "XPC mRNA expression in IMR-90", x = "Benzo[a]pyrene(μM)", y = "ΔCt(XPC)"),
       height = 10, width = 10, units = "cm", dpi = 600
)

ggsave("IMR_POLD4.png", 
       ggplot(imr_sta) + geom_errorbar(aes(x = conc, ymax = POLD4_sum+XPC_se, ymin = POLD4_sum-0.1), size = 0.7, width = 0.3) +
         geom_bar(aes(x = conc, y = POLD4_sum),stat = "identity", width = 0.6) +
         scale_y_continuous(expand = c(0, 0), limits = c(0, 17)) +
         labs(title = "POLD4 mRNA expression in IMR-90", x = "Benzo[a]pyrene(μM)", y = "ΔCt(POLD4)"),
       height = 10, width = 10, units = "cm", dpi = 600
)
