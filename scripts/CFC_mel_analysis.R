library(lme4)
library(lmerTest)
library(ggplot2)
library(ggbiplot)
library(gridExtra)
library(pls)
library(MASS)
library(dplyr)
library(writexl)
library(car)
library(stringr)
library(devtools)
library(psych)
library(RESI)
library(performance)
library(aod)


## set pathway for saving figures ##
path = "~/Desktop"


## Figure 2 ##

# code for Figure 2: a) PCA1 & PCA2, b) PCA3 & PCA4

#import data
chc <- read.csv("CHC.csv")
chc$rows <- interaction(chc$genotype, chc$replicate, chc$treatment)
rownames(chc) <- chc[,20]

#create PCA
chc.pca <- prcomp(chc[,c(4:19)], center = TRUE, scale. = TRUE)
summary(chc.pca)

#plot
chc.treatment <- c(rep("Con-perfumed", 20), rep("Het-perfumed", 20))
chc$legend <- substring(chc$rows, 1, 1)


png(filename=paste(path,"/Figure2.png",sep=""), width=13.5, height=6.5,
    units="in", res=300, pointsize=7)

p.PC12 <- ggbiplot(chc.pca, labels=chc$legend, var.axes=FALSE, 
                   groups = chc.treatment, labels.size = 5) +
  scale_colour_manual(name="Perfuming\ntreatment", values= c('gray70', 'slateblue3')) +
  theme_bw() +
  ylim(-2.6,2.3) +
  xlim(-2.6,2.3) +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), panel.grid = element_blank(),
        legend.text = element_text(size = 12), legend.title.align = 0.5, axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11), legend.title = element_text(size = 12), legend.position = "left",
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

p.PC34 <- ggbiplot(chc.pca, choices = 3:4, labels=chc$legend, var.axes=FALSE, 
                   groups = chc.treatment, labels.size = 5) +
  scale_colour_manual(name="Perfuming\ntreatment", values= c('gray70', 'slateblue3')) +
  theme_bw() +
  ylim(-2.7,2.2) +
  xlim(-2.7,2.2) +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), panel.grid = element_blank(),
        legend.text = element_text(size = 12), legend.title.align = 0.5, axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11), legend.title = element_text(size = 12),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))


grid.arrange(p.PC12, p.PC34, ncol=2)

dev.off()


# statistical analysis for Figure 2

#import data
chc <- read.csv("CHC.csv")
chc$rows <- interaction(chc$genotype, chc$replicate, chc$treatment)
rownames(chc) <- chc[,20]

#create PCA
chc.pca <- prcomp(chc[,c(4:19)], center = TRUE, scale. = TRUE)
summary(chc.pca)

# MANOVA
chc.pcst <- as.data.frame(chc.pca$x[,1:4])
chc.pcst <- as.data.frame(cbind(as.numeric(as.character(chc.pcst$PC1)),as.numeric(as.character(chc.pcst$PC2)),as.numeric(as.character(chc.pcst$PC3)),as.numeric(as.character(chc.pcst$PC4))))
chc.pcst <- cbind(chc.pcst, chc$genotype, chc$treatment)
colnames(chc.pcst) <- c("PC1", "PC2", "PC3", "PC4", "genotype", "treatment")

genotype=as.matrix(chc.pcst[,5])
treatment=as.matrix(chc.pcst[,6])
PC=as.matrix(chc.pcst[,1:4])

chc.manova=manova(PC~genotype + treatment)
summary(chc.manova, test="Wilks")
summary.aov(chc.manova)



## Figure 3 ##

# code for Figure 3: morphology, a) thorax length, b) reproductive trait length

#import data
thorax <- read.csv("thorax_length.csv")
length <- read.csv("SR_testis_length.csv")

# name factors & reorder levels
thorax$line <- factor(thorax$line,
                      levels = c("mel_GFP", "mel_C", "mel_G", "mel_M", "mel_T"),
                      labels = c("GFP", "California", "Greece", "Malaysia", "Taiwan"))
thorax$sex <- factor(thorax$sex,
                     levels = c("F", "M"),
                     labels = c("Female", "Male"))

length$line <- factor(length$line,
                      levels = c("mel_GFP", "mel_C", "mel_G", "mel_M", "mel_T"),
                      labels = c("GFP", "California", "Greece", "Malaysia", "Taiwan"))

length$trait <- factor(length$trait,
                       levels = c("SR", "testis", "sperm"),
                       labels = c("SR", "Testis", "Sperm"))

#plot

png(filename=paste(path,"/Figure3.png",sep=""), width=16.5, height=6.5,
    units="in", res=300, pointsize=7)

p.thor <- ggplot(data = thorax, mapping = aes(x = sex, y = thorax_length, fill = line)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1), show.legend = FALSE) +
  theme_bw() +
  geom_vline(xintercept= c(1.5), linetype = 2, color = "gray") +
  xlab("Sex") +
  ylab("Thorax length (mm)") +
  labs(fill = "Line") +
  scale_fill_viridis_d() +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), panel.grid = element_blank(),
        legend.text = element_text(size = 16), legend.title.align = 0.5, axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), legend.title = element_text(size = 16), legend.position = "left",
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

p.len <- ggplot(data = length, mapping = aes(x = trait, y = size, fill = line)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1), show.legend = FALSE) +
  theme_bw() +
  xlab("Trait") +
  ylab("Length (mm)") +
  labs(fill = "Line") +
  scale_fill_viridis_d() +
  geom_vline(xintercept= c(1.5,2.5), linetype = 2, color = "gray") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), panel.grid = element_blank(),
        legend.text = element_text(size = 16), legend.title.align = 0.5, axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), legend.title = element_text(size = 16),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

grid.arrange(p.thor, p.len, ncol=2)

dev.off()


# statistical analysis for Figure 3: morphology

#import data
morph <- read.csv("morphology.csv")

#re-level factors
morph$line <- factor(morph$line, levels=c('mel_C', 'mel_G', 'mel_M', 'mel_T', 'mel_GFP'))
print(levels(morph$line))

#set deviation-coded contrast
options(contrasts=c(unordered="contr.sum", ordered="contr.poly"))

#lm
m.thorm <- lm(thorax_length_M_mm ~ line, data = morph)
summary(m.thorm)

m.tes <- lm(testis_length_mm ~ thorax_length_M_mm + line, data = morph)
summary(m.tes)

m.sperm <- lm(sperm_length_mm ~ nonMRT_mass_mg + line, data = morph)
summary(m.sperm)

m.massm <- lm(nonMRT_mass_mg ~ line, data = morph)
summary(m.massm)

m.MRT <- lm(MRT_mass_mg ~ nonMRT_mass_mg + line, data = morph)
summary(m.MRT)

m.thorf <- lm(thorax_length_F_mm ~ line, data = morph)
summary(m.thorf)

m.SR <- lm(SR_length_mm ~ thorax_length_F_mm + line, data = morph)
summary(m.SR)

m.massf <- lm(nonFRT_mass_mg ~ line, data = morph)
summary(m.massf)

m.FRT <- lm(FRT_mass_mg ~ nonFRT_mass_mg + line, data = morph)
summary(m.FRT)

m.ST <- lm(ST_area_mm2 ~ thorax_length_F_mm + line, data = morph)
summary(m.ST)


## Figure 4 ##

# code for Figure 4: courtship effort

#import data
data.courtship <- read.csv("courtship.csv")
courtship <- subset(data.courtship, died == "N")

#label factors
courtship$treatment <- factor(courtship$treatment, levels=c("control", "treatment"), labels = c("Con-perfumed", "Het-perfumed"))

courtship$second_male <- factor(courtship$second_male,
                                levels = c("B4_California", "A5_Greece", "B7_Malaysia", "A7_Taiwan"),
                                labels = c("California", "Greece", "Malaysia", "Taiwan"))

courtship$female_line <- factor(courtship$female_line,
                                levels = c("B4_California", "A5_Greece", "B7_Malaysia", "A7_Taiwan"),
                                labels = c("California", "Greece", "Malaysia", "Taiwan"))

text_courtship <- data.frame(
  label = c("Female line", "Male line"),
  female_line = factor(c("Greece", "California")),
  second_male = factor(c("Taiwan", "Greece")),
  x = c(3.3, 2.5),
  y = c(0.05, 3.1),
  angle = c(270, 0)
)

courtship_means <- courtship %>%
  group_by(second_male) %>%
  summarise(
    courtship_mean = mean(courtship_over_time)
  )

png(filename=paste(path,"/Figure4.png",sep=""), width=10, height=8,
    units="in", res=300, pointsize=7)

ggplot(data = courtship, mapping = aes(x = treatment)) +
  facet_grid(rows = vars(female_line), cols = vars(second_male)) +
  geom_boxplot(data = courtship, mapping = aes(x = treatment, y = courtship_over_time, fill = treatment), outlier.shape = NA) +
  geom_hline(data = courtship_means, mapping = aes(yintercept = courtship_mean), color = "grey40", linetype = "dashed", linewidth = 0.8) +
  geom_jitter(data = courtship, mapping = aes(x = treatment, y = courtship_over_time, stroke = 0.7), width = 0.1, size = 1.5) +
  theme_bw() +
  ylab("Second male courtship effort") +
  scale_fill_manual(values = rep(c('gray80', 'slateblue3'), 4)) +
  geom_text(data = text_courtship, mapping = aes(x = x, y = y, label = label, angle = angle), size = 5.5) +
  coord_cartesian(xlim = c(1,2), ylim = c(0, 2.25), clip = 'off') +
  labs(fill = "Perfuming\ntreatment") +
  theme(axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.title.y = element_text(size = 20, margin = margin(0,10,0,0)), axis.text.x = element_blank(),
        legend.text = element_text(size = 12), legend.title.align = 0.5,
        axis.text.y = element_text(size = 11), legend.title = element_text(size = 12),
        axis.ticks.x = element_blank(), plot.margin = margin(1,1,0.5,0.5, "cm"),
        strip.text = element_text(size = 12), legend.justification = "top")

dev.off()


# statistical analysis for Figure 4: courtship effort

#import data
data.courtship <- read.csv("courtship.csv")
courtship <- subset(data.courtship, died == "N")

#re-level factors
courtship$female_line <- factor(courtship$female_line, levels=c('A5_Greece', 'B7_Malaysia', 'A7_Taiwan', 'B4_California'))
courtship$second_male <- factor(courtship$second_male, levels=c('A5_Greece', 'B7_Malaysia', 'A7_Taiwan', 'B4_California'))
courtship$treatment <- factor(courtship$treatment, levels=c('treatment', 'control'))

courtship$female_line <- factor(courtship$female_line, levels=c('B4_California', 'A5_Greece', 'B7_Malaysia', 'A7_Taiwan'))
courtship$second_male <- factor(courtship$second_male, levels=c('B4_California', 'A5_Greece', 'B7_Malaysia', 'A7_Taiwan'))

#set deviation-coded contrast
options(contrasts=c(unordered="contr.sum", ordered="contr.poly"))

#lm with female, male, treatment, & female x male effects
m.courtship <- lm(courtship_over_time ~ female_line*second_male + treatment, data = courtship)
summary(m.courtship)

#estimate effect sizes 
set.seed(0826)
courtship_es <- resi(m.courtship)
summary(courtship_es)
anova(courtship_es)

#anova
aov.courtship <- Anova(m.courtship, type = "3")
abs(log10(aov.courtship$`Pr(>F)`))


## Figure 5 ##

# code for Figure 5: premating success

#import data
data.vortex <- read.csv("paternity.csv")
data.trim <- subset(data.vortex, total >= 10)
data.alive <- subset(data.trim, died == "N")
remating <- subset(data.alive, before_remating > 0 | P1_total > 0)

#re-level and label factors
remating$second_male <- factor(remating$second_male,
                               levels = c("B4_California", "A5_Greece", "B7_Malaysia", "A7_Taiwan"),
                               labels = c("California", "Greece", "Malaysia", "Taiwan"))
remating$female_line <- factor(remating$female_line,
                               levels = c("B4_California", "A5_Greece", "B7_Malaysia", "A7_Taiwan"),
                               labels = c("California", "Greece", "Malaysia", "Taiwan"))

remating$treatment <- factor(remating$treatment, levels=c("control", "treatment"), labels = c("Con-perfumed", "Het-perfumed"))

remating$remated <- factor(remating$remated, levels=c("Y", "N"))

text_remating <- data.frame(
  label = c("Female line", "Male line"),
  female_line = factor(c("Greece", "California")),
  second_male = factor(c("Taiwan", "Greece")),
  x = c(3.3, 2.5),
  y = c(0.05, 19),
  angle = c(270, 0)
)

png(filename=paste(path,"/Figure5.png",sep=""), width=10, height=8,
    units="in", res=300, pointsize=7)

ggplot(data = remating, mapping = aes(x = treatment)) +
  facet_grid(rows = vars(female_line), cols = vars(second_male)) +
  geom_bar(data = remating, mapping = aes(x = treatment, fill = treatment, alpha = remated)) +
  scale_fill_manual(values = rep(c('gray70', 'slateblue3'), 4)) +
  scale_alpha_discrete(range = c(1,0.4), guide = "none") +
  theme_bw() +
  ylab("Mating successes/mating failures") +
  geom_text(data = text_remating, mapping = aes(x = x, y = y, label = label, angle = angle), size = 5.5) +
  coord_cartesian(xlim = c(1,2), ylim = c(0, 14), clip = 'off') +
  labs(fill = "Perfuming\ntreatment") +
  theme(axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.title.y = element_text(size = 20, margin = margin(0,10,0,0)), axis.text.x = element_blank(),
        legend.text = element_text(size = 12), legend.title.align = 0.5,
        axis.text.y = element_text(size = 11), legend.title = element_text(size = 12),
        axis.ticks.x = element_blank(), plot.margin = margin(1,1,0.5,0.5, "cm"),
        strip.text = element_text(size = 12), legend.justification = "top")

dev.off()


# statistical analysis for Figure 5: premating success

#import data
remating <- read.csv("remating.csv")

#re-level factors
remating$FemaleLine <- factor(remating$FemaleLine, levels=c('A5_Greece', 'B7_Malaysia', 'A7_Taiwan', 'B4_California'))
remating$SecondMaleLine <- factor(remating$SecondMaleLine, levels=c('A5_Greece', 'B7_Malaysia', 'A7_Taiwan', 'B4_California'))
remating$CHC <- factor(remating$CHC, levels=c('Treatment', 'Control'))

#set deviation-coded contrast
options(contrasts=c(unordered="contr.sum", ordered="contr.poly"))

#glm with female, male, treatment, & female x male effects
remated <- remating[,6:5]
remated <- data.matrix(remated)

m.remating <- glm(remated ~ FemaleLine + SecondMaleLine + CHC, data = remating, family = "binomial")
summary(m.remating)

#estimate effect sizes
set.seed(0826)
remating_es <- resi(m.remating)
summary(remating_es)
anova(remating_es)

#anova
aov.remating <- Anova(m.remating, type = "3")
abs(log10(aov.remating$`Pr(>Chisq)`))


## Figure 6 ##

# code for Figure 6: postmating success

#import data
data.vortex <- read.csv("paternity.csv")
data.trim <- subset(data.vortex, total >= 10)
data.alive <- subset(data.trim, died == "N")
data <- subset(data.alive, remated == "Y")

#re-level and label factors
data$second_male <- factor(data$second_male,
                           levels = c("B4_California", "A5_Greece", "B7_Malaysia", "A7_Taiwan"),
                           labels = c("California", "Greece", "Malaysia", "Taiwan"))
data$female_line <- factor(data$female_line,
                           levels = c("B4_California", "A5_Greece", "B7_Malaysia", "A7_Taiwan"),
                           labels = c("California", "Greece", "Malaysia", "Taiwan"))

data$treatment <- factor(data$treatment, levels=c("control", "treatment"), labels = c("Con-perfumed", "Het-perfumed"))

#plot
text_paternity <- data.frame(
  label = c("Female line", "Male line"),
  female_line = factor(c("Greece", "California")),
  second_male = factor(c("Taiwan", "Greece")),
  x = c(3.3, 2.5),
  y = c(0.05, 1.35),
  angle = c(270, 0)
)

P2_means <- data %>%
  group_by(second_male) %>%
  summarise(
    P2_mean = mean(P2)
  )


png(filename=paste(path,"/Figure6.png",sep=""), width=10, height=8,
    units="in", res=300, pointsize=7)

ggplot(data = data, mapping = aes(x = treatment)) +
  facet_grid(rows = vars(female_line), cols = vars(second_male)) +
  geom_boxplot(data = data, mapping = aes(x = treatment, y = P2, fill = treatment), outlier.shape = NA) +
  geom_hline(data = P2_means, mapping = aes(yintercept = P2_mean), color = "grey40", linetype = "dashed", linewidth = 0.8) +
  geom_jitter(data = data, mapping = aes(x = treatment, y = P2, stroke = 0.7), width = 0.1, size = 1.5) +
  theme_bw() +
  ylab("Second male siring success") +
  scale_fill_manual(values = rep(c('gray80', 'slateblue3'), 4)) +
  geom_text(data = text_paternity, mapping = aes(x = x, y = y, label = label, angle = angle), size = 5.5) +
  coord_cartesian(xlim = c(1,2), ylim = c(0, 1), clip = 'off') +
  labs(fill = "Perfuming\ntreatment") +
  theme(axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.title.y = element_text(size = 20, margin = margin(0,10,0,0)), axis.text.x = element_blank(),
        legend.text = element_text(size = 12), legend.title.align = 0.5,
        axis.text.y = element_text(size = 11), legend.title = element_text(size = 12),
        axis.ticks.x = element_blank(), plot.margin = margin(1,1,0.5,0.5, "cm"),
        strip.text = element_text(size = 12), legend.justification = "top")

dev.off()

# statistical analysis for Figure 6: postmating success

#import data
P2.vortex <- read.csv("paternity.csv")
P2.trim <- subset(P2.vortex, total >= 10)
P2.alive <- subset(P2.trim, died == "N")
P2 <- subset(P2.alive, remated == "Y")

#create data matrix for binomial response
paternity <- P2[,23:22]
paternity <- data.matrix(paternity)

#re-level factors
P2$female_line <- factor(P2$female_line, levels=c('A5_Greece', 'B7_Malaysia', 'A7_Taiwan', 'B4_California'))
P2$second_male <- factor(P2$second_male, levels=c('A5_Greece', 'B7_Malaysia', 'A7_Taiwan', 'B4_California'))
P2$treatment <- factor(P2$treatment, levels=c('treatment', 'control'))

#set deviation-coded contrast
options(contrasts=c(unordered="contr.sum", ordered="contr.poly"))

#glm witih female, male, treatment, female x male, & female x treatment effects
m.paternity <- glm(paternity ~ female_line*second_male + treatment + female_line:treatment, data = P2, family = "binomial")
summary(m.paternity)

#estimate effect sizes
set.seed(0826)
paternity_es <- resi(m.paternity)
summary(paternity_es)
anova(paternity_es)

#anova
aov.paternity <- Anova(m.paternity, type = "3")
abs(log10(aov.paternity$`Pr(>Chisq)`))


## Figure 7 ##

# code for Figure 7: comparing coefficients across datasets

#import data
mean <- read.csv("means.csv")

#split data
mean_con <- subset(mean, treatment == "c")
mean_treat <- subset(mean, treatment == "t")

#plot
png(filename=paste(path,"/Figure7.png",sep=""), width=12, height=18,
    units="in", res=300, pointsize=7)

#male courtship vs. siring success, control
p.comp1 <- ggplot(data = mean_con, mapping = aes(x = courtship_mean, y = P2_mean)) +
  geom_text(data = mean_con, mapping = aes(x = courtship_mean, y = P2_mean, label = cross_type, col = treatment), size = 6, show.legend = FALSE) +
  theme_bw() +
  ylab("Mean second male siring success") +
  xlab("Mean courtship effort per minute") +
  ylim(0, 1) +
  xlim(0, 1.65) +
  scale_color_manual(values = c('gray70')) +
  labs(tag = "(a)") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.grid = element_blank(), axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.tag = element_text(size = 18))

#male courtship vs. siring success, treatment
p.comp2 <- ggplot(data = mean_treat, mapping = aes(x = courtship_mean, y = P2_mean)) +
  geom_text(data = mean_treat, mapping = aes(x = courtship_mean, y = P2_mean, label = cross_type, col = treatment), size = 6, show.legend = FALSE) +
  theme_bw() +
  ylab("Mean second male siring success") +
  xlab("Mean courtship effort per minute") +
  ylim(0, 1) +
  xlim(0, 1.65) +
  scale_color_manual(values = c('slateblue3')) +
  labs(tag = "(b)") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.grid = element_blank(), axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.tag = element_text(size = 18))

#female remating vs. siring success, control
p.comp3 <- ggplot(data = mean_con, mapping = aes(x = remating_prop, y = P2_mean)) +
  geom_text(data = mean_con, mapping = aes(x = remating_prop, y = P2_mean, label = cross_type, col = treatment), size = 6, show.legend = FALSE) +
  theme_bw() +
  ylab("Mean second male siring success") +
  xlab("Percent female remating") +
  ylim(0, 1) +
  xlim(0.15, 1.05) +
  scale_color_manual(values = c('gray70')) +
  labs(tag = "(c)") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.grid = element_blank(), axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.tag = element_text(size = 18))

#female remating vs. siring success, treatment
p.comp4 <- ggplot(data = mean_treat, mapping = aes(x = remating_prop, y = P2_mean)) +
  geom_text(data = mean_treat, mapping = aes(x = remating_prop, y = P2_mean, label = cross_type, col = treatment), size = 6, show.legend = FALSE) +
  theme_bw() +
  ylab("Mean second male siring success") +
  xlab("Percent female remating") +
  ylim(0, 1) +
  xlim(0.15, 1.05) +
  scale_color_manual(values = c('slateblue3')) +
  labs(tag = "(d)") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.grid = element_blank(), axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.tag = element_text(size = 18))

#sperm length vs. paternity, control
p.comp5 <- ggplot(data = mean_con, mapping = aes(x = sperm_mean, y = P2_mean)) +
  geom_text(data = mean_con, mapping = aes(x = sperm_mean, y = P2_mean, label = cross_type, col = treatment), size = 6, show.legend = FALSE) +
  theme_bw() +
  ylab("Mean second male siring success") +
  xlab("Mean sperm length (mm)") +
  ylim(0, 1) +
  xlim(1.55, 1.85) +
  scale_color_manual(values = c('gray70')) +
  labs(tag = "(e)") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.grid = element_blank(), axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.tag = element_text(size = 18))

#sperm length vs. paternity, treatment
p.comp6 <- ggplot(data = mean_treat, mapping = aes(x = sperm_mean, y = P2_mean)) +
  geom_text(data = mean_treat, mapping = aes(x = sperm_mean, y = P2_mean, label = cross_type, col = treatment), size = 6, show.legend = FALSE) +
  theme_bw() +
  ylab("Mean second male siring success") +
  xlab("Mean sperm length (mm)") +
  ylim(0, 1) +
  xlim(1.55, 1.85) +
  scale_color_manual(values = c('slateblue3')) +
  labs(tag = "(f)") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.grid = element_blank(), axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.tag = element_text(size = 18))


grid.arrange(p.comp1, p.comp2, p.comp3, p.comp4, p.comp5, p.comp6, nrow=3, ncol=2)

dev.off()

#statistical analysis for Figure 7: comparisons across datasets

#import data
mean <- read.csv("means.csv")

#split data
mean_con <- subset(mean, treatment == "c")
mean_treat <- subset(mean, treatment == "t")

#Spearman's rank correlation rho
#male courtship vs. female remating
cor.test(mean$courtship_mean, mean$remating_prop, method = "spearman")

cor.test(mean_con$courtship_mean, mean_con$remating_prop, method = "spearman")
cor.test(mean_treat$courtship_mean, mean_treat$remating_prop, method = "spearman")

#male courtship vs. siring success
cor.test(mean$courtship_mean, mean$P2_mean, method = "spearman")

cor.test(mean_con$courtship_mean, mean_con$P2_mean, method = "spearman")
cor.test(mean_treat$courtship_mean, mean_treat$P2_mean, method = "spearman")

x <- r.test(16, 0.535, 0.353)
print(x, digits = 4)

#female remating vs. siring success
cor.test(mean$remating_prop, mean$P2_mean, method = "spearman")

cor.test(mean_con$remating_prop, mean_con$P2_mean, method = "spearman")
cor.test(mean_treat$remating_prop, mean_treat$P2_mean, method = "spearman")

x <- r.test(16, -0.095, 0.502)
print(x, digits = 4)

#sperm length vs. siring success
cor.test(mean$sperm_mean, mean$P2_mean, method = "spearman")

cor.test(mean_con$sperm_mean, mean_con$P2_mean, method = "spearman")
cor.test(mean_treat$sperm_mean, mean_treat$P2_mean, method = "spearman")

#sperm length minus tester sperm vs. siring success
cor.test(mean$sperm_minus_tester, mean$P2_mean, method = "spearman")

cor.test(mean_con$sperm_minus_tester, mean_con$P2_mean, method = "spearman")
cor.test(mean_treat$sperm_minus_tester, mean_treat$P2_mean, method = "spearman")

#testis length vs. siring success
cor.test(mean$testis_mean, mean$P2_mean, method = "spearman")

cor.test(mean_con$testis_mean, mean_con$P2_mean, method = "spearman")
cor.test(mean_treat$testis_mean, mean_treat$P2_mean, method = "spearman")

#SR length vs. siring success
cor.test(mean$SR_mean, mean$P2_mean, method = "spearman")

cor.test(mean_con$SR_mean, mean_con$P2_mean, method = "spearman")
cor.test(mean_treat$SR_mean, mean_treat$P2_mean, method = "spearman")

#SR length minus sperm length vs. siring success
cor.test(mean$SR_minus_sperm, mean$P2_mean, method = "spearman")

cor.test(mean_con$SR_minus_sperm, mean_con$P2_mean, method = "spearman")
cor.test(mean_treat$SR_minus_sperm, mean_treat$P2_mean, method = "spearman")

#MRT mass vs. siring success
cor.test(mean$MRT_mean, mean$P2_mean, method = "spearman")

cor.test(mean_con$MRT_mean, mean_con$P2_mean, method = "spearman")
cor.test(mean_treat$MRT_mean, mean_treat$P2_mean, method = "spearman")

x <- r.test(5, -0.449, -0.679)
print(x, digits = 4)

#FRT mass vs. siring success
cor.test(mean$FRT_mean, mean$P2_mean, method = "spearman")

cor.test(mean_con$FRT_mean, mean_con$P2_mean, method = "spearman")
cor.test(mean_treat$FRT_mean, mean_treat$P2_mean, method = "spearman")



### Supplemental figures

## Figure S1 ##

# code for Figure S1: ST area, non-RT body mass, RT mass

#import data
morph <- read.csv("morphology.csv")

# name factors & reorder levels
morph$line <- factor(morph$line,
                     levels = c("mel_GFP", "mel_C", "mel_G", "mel_M", "mel_T"),
                     labels = c("GFP", "California", "Greece", "Malaysia", "Taiwan"))

#import data
area <- read.csv("area.csv")
mass <- read.csv("mass.csv")

# name factors & reorder levels
mass$line <- factor(mass$line,
                    levels = c("mel_GFP", "mel_C", "mel_G", "mel_M", "mel_T"),
                    labels = c("GFP", "California", "Greece", "Malaysia", "Taiwan"))
mass$trait <- factor(mass$trait,
                     levels = c("FRT", "MRT"),
                     labels = c("Female", "Male"))

area$line <- factor(area$line,
                    levels = c("mel_GFP", "mel_C", "mel_G", "mel_M", "mel_T"),
                    labels = c("GFP", "California", "Greece", "Malaysia", "Taiwan"))

#plot

png(filename=paste(path,"/FigureS1.png",sep=""), width=16.5, height=13,
    units="in", res=300, pointsize=7)

p.body <- ggplot(data = mass, mapping = aes(x = trait, y = total_mass, fill = line)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1), show.legend = FALSE) +
  theme_bw() +
  geom_vline(xintercept= c(1.5), linetype = 2, color = "gray") +
  xlab("Sex") +
  ylab("Total dry body mass (mg)") +
  labs(fill = "Line") +
  scale_fill_viridis_d() +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), panel.grid = element_blank(),
        legend.text = element_text(size = 16), legend.title.align = 0.5, axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), legend.title = element_text(size = 16), legend.position = "left",
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

p.RT <- ggplot(data = mass, mapping = aes(x = trait, y = RT_mass, fill = line)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1), show.legend = FALSE) +
  theme_bw() +
  xlab("Sex") +
  ylab("Reproductive tract mass (mg)") +
  labs(fill = "Line") +
  scale_fill_viridis_d() +
  geom_vline(xintercept= c(1.5), linetype = 2, color = "gray") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), panel.grid = element_blank(),
        legend.text = element_text(size = 16), legend.title.align = 0.5, axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), legend.title = element_text(size = 16),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

p.ST <- ggplot(data = area, mapping = aes(x = trait, y = size, fill = line)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1), show.legend = FALSE) +
  theme_bw() +
  ylab(expression(paste("ST area (mm"^"2", ")"))) +
  labs(fill = "Line") +
  scale_fill_viridis_d() +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20), panel.grid = element_blank(),
        legend.text = element_text(size = 16), legend.title.align = 0.5, axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14), legend.title = element_text(size = 16), legend.position = "left",
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

grid.arrange(p.body, p.RT, p.ST, nrow=2, ncol=2)

dev.off()


## Figure S2 ##

# code for Figure S2: courtship effort, non-truncated y-axis

#import data
data.courtship <- read.csv("courtship.csv")
courtship <- subset(data.courtship, died == "N")

#label factors
courtship$treatment <- factor(courtship$treatment, levels=c("control", "treatment"), labels = c("Con-perfumed", "Het-perfumed"))

courtship$second_male <- factor(courtship$second_male,
                                levels = c("B4_California", "A5_Greece", "B7_Malaysia", "A7_Taiwan"),
                                labels = c("California", "Greece", "Malaysia", "Taiwan"))

courtship$female_line <- factor(courtship$female_line,
                                levels = c("B4_California", "A5_Greece", "B7_Malaysia", "A7_Taiwan"),
                                labels = c("California", "Greece", "Malaysia", "Taiwan"))

text_courtship <- data.frame(
  label = c("Female line", "Male line"),
  female_line = factor(c("Greece", "California")),
  second_male = factor(c("Taiwan", "Greece")),
  x = c(3.3, 2.5),
  y = c(0.05, 5.7),
  angle = c(270, 0)
)

#plot
png(filename=paste(path,"/FigureS2.png",sep=""), width=10, height=8,
    units="in", res=300, pointsize=7)

ggplot(data = courtship, mapping = aes(x = treatment)) +
  facet_grid(rows = vars(female_line), cols = vars(second_male)) +
  geom_boxplot(data = courtship, mapping = aes(x = treatment, y = courtship_over_time, fill = treatment), outlier.shape = NA) +
  geom_hline(data = courtship_means, mapping = aes(yintercept = courtship_mean), color = "grey40", linetype = "dashed", linewidth = 0.8) +
  geom_jitter(data = courtship, mapping = aes(x = treatment, y = courtship_over_time, stroke = 0.7), width = 0.1, size = 1.5) +
  theme_bw() +
  ylab("Second male courtship effort") +
  scale_fill_manual(values = rep(c('gray80', 'slateblue3'), 4)) +
  geom_text(data = text_courtship, mapping = aes(x = x, y = y, label = label, angle = angle), size = 5.5) +
  coord_cartesian(xlim = c(1,2), ylim = c(0, 4.2), clip = 'off') +
  labs(fill = "Perfuming\ntreatment") +
  theme(axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.title.y = element_text(size = 20, margin = margin(0,10,0,0)), axis.text.x = element_blank(),
        legend.text = element_text(size = 12), legend.title.align = 0.5,
        axis.text.y = element_text(size = 11), legend.title = element_text(size = 12),
        axis.ticks.x = element_blank(), plot.margin = margin(1,1,0.5,0.5, "cm"),
        strip.text = element_text(size = 12), legend.justification = "top")

dev.off()


## Figure S3 ##

# code for Figure S3: comparing coefficients, male courtship vs. female remating & MRT mass vs. siring success

#import data
mean <- read.csv("means.csv")

#split data
mean_con <- subset(mean, treatment == "c")
mean_treat <- subset(mean, treatment == "t")

#plot
png(filename=paste(path,"/FigureS3.png",sep=""), width=12, height=12,
    units="in", res=300, pointsize=7)

#male courtship vs. remating, control
p.comp7 <- ggplot(data = mean_con, mapping = aes(x = courtship_mean, y = remating_prop)) +
  geom_text(data = mean_con, mapping = aes(x = courtship_mean, y = remating_prop, label = cross_type, col = treatment), size = 6, show.legend = FALSE) +
  theme_bw() +
  ylab("Percent female remating") +
  xlab("Mean courtship effort per minute") +
  ylim(0, 1) +
  xlim(0.1, 1.6) +
  scale_color_manual(values = c('gray70')) +
  labs(tag = "(a)") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.grid = element_blank(), axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.tag = element_text(size = 18))

#male courtship vs. remating, treatment
p.comp8 <- ggplot(data = mean_treat, mapping = aes(x = courtship_mean, y = remating_prop)) +
  geom_text(data = mean_treat, mapping = aes(x = courtship_mean, y = P2_mean, label = cross_type, col = treatment), size = 6, show.legend = FALSE) +
  theme_bw() +
  ylab("Percent female remating") +
  xlab("Mean courtship effort per minute") +
  ylim(0, 1) +
  xlim(0.1, 1.6) +
  scale_color_manual(values = c('slateblue3')) +
  labs(tag = "(b)") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.grid = element_blank(), axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.tag = element_text(size = 18))

#MRT mass vs. siring success, control
p.comp9 <- ggplot(data = mean_con, mapping = aes(x = MRT_mean, y = P2_mean)) +
  geom_text(data = mean_con, mapping = aes(x = MRT_mean, y = P2_mean, label = cross_type, col = treatment), size = 6, show.legend = FALSE) +
  theme_bw() +
  ylab("Mean second male siring success") +
  xlab("Mean MRT mass (mg)") +
  ylim(0, 1) +
  xlim(0.02, 0.031) +
  scale_color_manual(values = c('gray70')) +
  labs(tag = "(c)") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.grid = element_blank(), axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.tag = element_text(size = 18))

#MRT mass vs. siring success, treatment
p.comp10 <- ggplot(data = mean_treat, mapping = aes(x = MRT_mean, y = P2_mean)) +
  geom_text(data = mean_treat, mapping = aes(x = MRT_mean, y = P2_mean, label = cross_type, col = treatment), size = 6, show.legend = FALSE) +
  theme_bw() +
  ylab("Mean second male siring success") +
  xlab("Mean MRT mass (mg)") +
  ylim(0, 1) +
  xlim(0.02, 0.031) +
  scale_color_manual(values = c('slateblue3')) +
  labs(tag = "(d)") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.grid = element_blank(), axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.tag = element_text(size = 18))

grid.arrange(p.comp7, p.comp8, p.comp9, p.comp10, nrow=2, ncol=2)

dev.off()


