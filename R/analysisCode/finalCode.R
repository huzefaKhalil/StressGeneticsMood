# load libs & data
library(data.table)
library(ggplot2)
library(ggpubr)
library(plotrix)
library(scales)
library(kableExtra)
library(epitools)
library(rstatix)
library(ggthemes)
library(lme4)
library(lmerTest)
library(stargazer)
library(glmmTMB)
library(FactoMineR)
library(factoextra)
library(missMDA)
library(emmeans)
library(ggeffects)
library(caret)
library(ggsci)
library(svglite)
library(systemfonts)
library(zoo)
library(webr)
library(pROC)
library(modelsummary)
library(scales)

source("code/utilityFunctions.R")

combinedData <- filterCompleteness(loadCombinedData())
#combinedData <- filterCompleteness(loadCombinedData(caucasiansOnly = FALSE))

# change it to CT PRS. This is the final PRS we will use.
combinedData$PRS <- combinedData$PRS.C.1
combinedData$rawPRS <- combinedData$rawPRS.C.1

moodData <- loadMoodData(combinedData)

# remove screening from mood data
moodData[[2]] <- droplevels(moodData[[2]][Event != "May 2019"])
moodData[[3]] <- droplevels(moodData[[3]][Event != "May 2020"])

combinedData$PRSHiLo <- ifelse(combinedData$PRS > 0, "High PRS", "Low PRS")

# get the combined mood for 4 follow-up time-points
cMood <- moodData$combined
cMood <- merge(cMood, combinedData[, c("ID", "Sex")], all.x = TRUE)

cMood$tLin <- poly(as.integer(cMood$Event), 3)[, 1]
cMood$tQuad <- poly(as.integer(cMood$Event), 3)[, 2]
cMood$tCube <- poly(as.integer(cMood$Event), 3)[, 3]

# final snps after clumping
snps <- c(4867, 4917, 4222, 3940, 3838, 3681, 3423, 3173, 2892, 3182, 2973,
          3141, 2353, 2113, 1928, 2054, 1870, 1947, 1474, 1673, 961, 964)

# get activity data
mActData <- loadMonthlyActivityData(combinedData)
qActData <- loadQuarterlyActivityData(combinedData)
dActData <- loadDailyActivityData(combinedData)

combinedData$Phase <- factor(combinedData$Phase,
                             levels = c("Phase1", "Phase2", "Phase3"),
                             labels = c("Phase 1", "Phase 2", "Phase 3"))

bestModSurveys <- c("SBState", "SBTrait", "CTQ", "GAD7",
                    "PANSIPos", "NEOPIR", "PSS", "RFQ",
                    "PANSINeg", "PHQ9")

pcOutBest <- PCA(imputePCA(combinedData[, ..bestModSurveys])$completeObs, graph = FALSE)

cData <- cbind(combinedData, data.frame(bestZscore.1 = scale(pcOutBest$ind$coord[, 1]),
                                        bestZscore.2 = scale(pcOutBest$ind$coord[, 2]),
                                        dim.1 = pcOutBest$ind$coord[, 1],
                                        dim.2 = pcOutBest$ind$coord[, 2]))

##########################################################################
# Figure 1
##########################################################################


myComparisons <- list(c("Phase 1", "Phase 2"), c("Phase 2", "Phase 3"), c("Phase 1", "Phase 3"))

# TODO Change label.y and expansion
makeFig1BoxPlotsByPhase <- function(data, y, ylab, fname,
                                    ylim,
                                    label.y = list(c(4.2, 4.7, 5.3),  c(16, 18, 20)),
                                    expansion = list(c(0, 0.065), c(0, 0.065))) {
  
  tf <- as.formula(paste(y, "~ Phase"))
  cm <- compare_means(tf, p.adjust.method = "fdr", data = data)
  cm <- cm[c(1, 3, 2), ]
  
  gBox <- ggboxplot(data = data, x = "Phase", y = y,
                    color = "Phase", palette = "frontiers",
                    xlab = FALSE, ylab = ylab, ylim = c(-1, ylim[2]),
                    legend = "none", add = "jitter", add.params = list(size = 0.01)) +
    ylim(0, ylim[2]) +
    stat_pvalue_manual(cm,
                       label = c("p.adj"),
                       y.position = label.y[[2]],
                       label.sep = " ",
                       hide.ns = FALSE,
                       tip.length = 0.001, size = 1.93) +
    scale_y_continuous(expand = expansion(mult = expansion[[2]]),
                       breaks = seq(from = 0, to = ylim[2], by = 10)) +
    theme(text = element_text(size = 6, family = "Arial", face = "plain"),
          axis.text.x = element_text(angle = 30, vjust = 0.8),
          axis.line = element_line(linewidth = 0.3, color = "black"),
          axis.ticks = element_line(linewidth = 0.3, color = "black"),
          plot.margin = unit(c(4, 6, 0, 6), 'pt'))
  
  ggsave(paste0("finalPlots/", fname, ".Boxplot.tiff"), gBox, width = 1.4, height = 1.4, units = "in", dpi=320)
  
}

makeFig1BoxPlotsByPhase(combinedData, "PHQ9", "Baseline PHQ-9", "Fig1/adjP/BL.PHQ9",
                        label.y = list(c(5.1, 5.7, 6.4), c(18, 20, 22)),
                        ylim = c(8, 25))
makeFig1BoxPlotsByPhase(combinedData[Sex == "Female"], "PHQ9", "Baseline PHQ-9", "Fig1/adjP/BL.PHQ9.Female",
                        label.y = list(c(5.1, 5.7, 6.4), c(18, 20, 22)),
                        ylim = c(8, 25))
makeFig1BoxPlotsByPhase(combinedData[Sex == "Male"], "PHQ9", "Baseline PHQ-9", "Fig1/adjP/BL.PHQ9.Male",
                        label.y = list(c(5.1, 5.7, 6.4), c(18, 20, 22)),
                        ylim = c(8, 25))

makeFig1BoxPlotsByPhase(combinedData, "GAD7", "Baseline GAD-7", "Fig1/adjP/BL.GAD7",
                        label.y = list(c(5.2, 5.7, 6.3), c(17, 19, 21)),
                        ylim = c(7.5, 22))
makeFig1BoxPlotsByPhase(combinedData[Sex == "Female"], "GAD7", "Baseline GAD-7", "Fig1/adjP/BL.GAD7Female",
                        label.y = list(c(5.2, 5.7, 6.3), c(17, 19, 21)),
                        ylim = c(7.5, 22))
makeFig1BoxPlotsByPhase(combinedData[Sex == "Male"], "GAD7", "Baseline GAD-7", "Fig1/adjP/BL.GAD7.Male",
                        label.y = list(c(5.2, 5.7, 6.3), c(17, 19, 21)),
                        ylim = c(7.5, 22))

makeFig1BoxPlotsByPhase(combinedData, "FU4TPPHQ9", "Follow-up PHQ-9", "Fig1/adjP/FU.PHQ9",
                        label.y = list(c(8.5, 9.5, 10.5), c(21.5, 23.75, 26)),
                        ylim = c(12, 27))
makeFig1BoxPlotsByPhase(combinedData[Sex == "Female"], "FU4TPPHQ9", "Follow-up PHQ-9", "Fig1/adjP/FU.PHQ9.Female",
                        label.y = list(c(8.5, 9.5, 10.5), c(21.5, 23.75, 26)),
                        ylim = c(12.2, 27))
makeFig1BoxPlotsByPhase(combinedData[Sex == "Male"], "FU4TPPHQ9", "Follow-up PHQ-9", "Fig1/adjP/FU.PHQ9.Male",
                        label.y = list(c(8.5, 9.5, 10.5), c(21.5, 23.75, 26)),
                        ylim = c(12.2, 27))

makeFig1BoxPlotsByPhase(combinedData, "FU4TPGAD7", "Follow-up GAD-7", "Fig1/adjP/FU.GAD7",
                        label.y = list(c(7.7, 8.7, 9.7), c(20, 22, 24)),
                        ylim = c(11, 25))
makeFig1BoxPlotsByPhase(combinedData[Sex == "Female"], "FU4TPGAD7", "Follow-up GAD-7", "Fig1/adjP/FU.GAD7Female",
                        label.y = list(c(7.7, 8.7, 9.7), c(20, 22, 24)),
                        ylim = c(11, 25))
makeFig1BoxPlotsByPhase(combinedData[Sex == "Male"], "FU4TPGAD7", "Follow-up GAD-7", "Fig1/adjP/FU.GAD7.Male",
                        label.y = list(c(7.7, 8.7, 9.7), c(20, 22, 24)),
                        ylim = c(11, 25))

# now the persistence figure
# let's look at persistence by Phase
moodCombined <- moodData$combined[Event != "Sep BL"]

depTimes <- moodCombined[, list(depTimes = sum(PHQ9 >= 10, na.rm = TRUE),
                                persistence = sum(PHQ9 >= 10, na.rm = TRUE) > 1,
                                anyDep = any(PHQ9 >= 10, na.rm = TRUE)),
                         by = "ID"]

# combine it
tData <- merge(cData, depTimes, all.y = FALSE)

table(tData$persistence, tData$anyDep, tData$Sex, tData$Phase)
table(tData$persistence, tData$Sex, tData$Phase)
table(tData$persistence, tData$Sex)
table(tData$anyDep, tData$Sex, tData$Phase)

persistence <- data.table(
  Phase = c(rep("Phase 1", 2), rep("Phase 2", 2), rep("Phase 3", 2)),
  Sex = rep(c("Female", "Male"), 3),
  PercentPersistent = c(2/47, 2/58, 13/80, 4/36, 20/71, 2/39) * 100,
  PercentOfDep = c(2/10, 2/7, 13/28, 4/7, 20/31, 2/7) * 100,
  nPersistent = c(2, 2, 13, 4, 20, 2)
)

persistence2 <- data.table(
  Phase = c("Phase 1", "Phase 2", "Phase 3"),
  PercentPersistent = c(4/105, 17/116, 22/110) * 100
)

g1 <- ggbarplot(data = persistence2, x = "Phase", y = "PercentPersistent",
                color = "Phase", palette = "frontiers", fill = "Phase",
                xlab = FALSE, ylab = "Percent with Chronic Dep",
                ylim = c(0, 30),
                legend = "none", add.params = list(size = 0.01)) +
  scale_y_continuous(expand = expansion(mult = 0),
                     breaks = seq(from = 0, to = 30, by = 10)) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'))

g2 <- ggbarplot(data = persistence[Sex == "Female"], x = "Phase", y = "PercentPersistent",
                color = "Phase", palette = "frontiers", fill = "Phase",
                xlab = FALSE, ylab = "Percent with Chronic Dep",
                ylim = c(0, 30),
                legend = "none", add.params = list(size = 0.01)) +
  scale_y_continuous(expand = expansion(mult = 0),
                     breaks = seq(from = 0, to = 30, by = 10)) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'))

g3 <- ggbarplot(data = persistence[Sex == "Male"], x = "Phase", y = "PercentPersistent",
                color = "Phase", palette = "frontiers", fill = "Phase",
                xlab = FALSE, ylab = "Percent with Chronic Dep",
                ylim = c(0, 30),
                legend = "none", add.params = list(size = 0.01)) +
  scale_y_continuous(expand = expansion(mult = 0),
                     breaks = seq(from = 0, to = 30, by = 10)) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'))

ggsave("finalPlots/Fig1/PersistenceAll.tiff", g1, width = 1.4, height = 1.4, units = "in", dpi=320)
ggsave("finalPlots/Fig1/PersistenceFemale.tiff", g2,width = 1.4, height = 1.4, units = "in", dpi=320)
ggsave("finalPlots/Fig1/PersistenceMale.tiff", g3, width = 1.4, height = 1.4, units = "in", dpi=320)

epitab(tData[Phase == "Phase1"]$Sex, tData[Phase == "Phase1"]$persistence, rev = "both")
epitab(tData[Phase == "Phase2"]$Sex, tData[Phase == "Phase2"]$persistence, rev = "both")
epitab(tData[Phase == "Phase3"]$Sex, tData[Phase == "Phase3"]$persistence, rev = "both")

epitab(tData$Phase, tData$persistence, rev = "columns")
epitab(tData$Phase, tData$persistence, rev = "both")

epitab(tData[Sex == "Female"]$Phase, tData[Sex == "Female"]$persistence, rev = "columns")
epitab(tData[Sex == "Female"]$Phase, tData[Sex == "Female"]$persistence, rev = "both")

epitab(tData[Sex == "Male"]$Phase, tData[Sex == "Male"]$persistence, rev = "columns")
epitab(tData[Sex == "Male"]$Phase, tData[Sex == "Male"]$persistence, rev = "both")


##########################################################################
# Figure 2
##########################################################################

#Ns
# PHQ & GAD = 116
# Activity = 80
# Sleep = 75

# load pre, post covid data
p2Gad <- fread("data/phase2/combined/gad.csv")
p2Phq <- fread("data/phase2/combined/phq.csv")

# remove screening
p2Gad <- p2Gad[Event != "May 2019"]
p2Phq <- p2Phq[Event != "May 2019"]

p2Gad <- p2Gad[ID %in% combinedData$ID]
p2Phq <- p2Phq[ID %in% combinedData$ID]

p2Gad$covid <- ifelse(p2Gad$Event %in% c("Sep 2019", "Dec 2019"), "Pre-COVID", "During COVID")
p2Phq$covid <- ifelse(p2Phq$Event %in% c("Sep 2019", "Dec 2019"), "Pre-COVID", "During COVID")

p2SubsetGad <- p2Gad[Event %in% c("Sep 2019", "Dec 2019", "Mar 2020", "Jun 2020", "Sep 2020")]
p2SubsetPhq <- p2Phq[Event %in% c("Sep 2019", "Dec 2019", "Mar 2020", "Jun 2020", "Sep 2020")]

p2GadCovid <- p2SubsetGad[, .(maxGad = max(gadSum, na.rm = TRUE)),
                          by = c("ID", "covid")]
p2PhqCovid <- p2SubsetPhq[, .(maxPhq = max(phqSum, na.rm = TRUE)),
                          by = c("ID", "covid")]

# make the box plots
gGad <- ggboxplot(data = p2GadCovid, x = "covid", y = "maxGad",
                  color = "covid", palette = c("darkorange", "purple"),
                  xlab = FALSE, ylab = "GAD-7", ylim = c(-1, 25),
                  legend = "none", add = "jitter", add.params = list(size = 0.01)) +
  stat_compare_means(comparisons = list(c("Pre-COVID", "During COVID")),
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = 21,
                     paired = TRUE,
                     tip.length = 0.01, size = 1.93) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.065)),
                     breaks = seq(from = 0, to = 25, by = 10)) +
  scale_x_discrete(labels = c("Pre-COVID", "During\nCOVID")) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'))

cairo_pdf(paste0("finalPlots/Fig2/GADCovid.Boxplot.pdf"), width = 1.4, height = 1.4, family = "Arial")
print(gGad)
dev.off()

svglite(paste0("finalPlots/Fig2/GADCovid.Boxplot.svg"), width = 1.4, height = 1.4, system_fonts = list(sans = "Arial"))
print(gGad)
dev.off()

ggsave(paste0("finalPlots/Fig2/GADCovid.Boxplot.tiff"), gGad, width = 1.4, height = 1.4, units = "in", dpi=320)

gPhq <- ggboxplot(data = p2PhqCovid, x = "covid", y = "maxPhq",
                  color = "covid", palette = c("darkorange", "purple"),
                  xlab = FALSE, ylab = "PHQ-9", ylim = c(-1, 25),
                  legend = "none", add = "jitter", add.params = list(size = 0.01)) +
  stat_compare_means(comparisons = list(c("Pre-COVID", "During COVID")),
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = 21,
                     paired = TRUE,
                     tip.length = 0.01, size = 1.93) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.065)),
                     breaks = seq(from = 0, to = 25, by = 10)) +
  scale_x_discrete(labels = c("Pre-COVID", "During\nCOVID")) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'))

cairo_pdf(paste0("finalPlots/Fig2/PHQCovid.Boxplot.pdf"), width = 1.4, height = 1.4, family = "Arial")
print(gPhq)
dev.off()

svglite(paste0("finalPlots/Fig2/PHQCovid.Boxplot.svg"), width = 1.4, height = 1.4, system_fonts = list(sans = "Arial"))
print(gPhq)
dev.off()

ggsave(paste0("finalPlots/Fig2/PHQCovid.Boxplot.tiff"), gPhq, width = 1.4, height = 1.4, units = "in", dpi=320)

# now with activity
covidAct <- dActData[Cohort == "Cohort5"]

covidAct$covid <- ifelse(covidAct$Date < "2020-03-01", "Pre-COVID", "Post-COVID")

covidAct <- covidAct[,  list(avgSteps = mean(StepTotal, na.rm = TRUE),
                             avgSleep = mean(TotalMinutesAsleep, na.rm = TRUE),
                             avgTimeInBed = mean(TotalTimeInBed, na.rm = TRUE),
                             stepsCount = length(StepTotal) - sum(is.na(StepTotal)),
                             sleepCount = length(TotalMinutesAsleep) - sum(is.na(TotalMinutesAsleep))),
                     by = c("ID", "covid", "Sex")]

covidAct <- covidAct[ID %in% covidAct[stepsCount > 20]$ID]
covidAct <- covidAct[ID %in% covidAct[avgSteps > 100]$ID]

gSteps <- ggboxplot(data = covidAct, x = "covid", y = "avgSteps",
                    color = "covid", palette = c("darkorange", "purple"),
                    xlab = FALSE, ylab = "Avg Steps / Day", ylim = c(-300, 17500),
                    legend = "none", add = "jitter", add.params = list(size = 0.01)) +
  stat_compare_means(comparisons = list(c("Pre-COVID", "Post-COVID")),
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = 16500,
                     paired = TRUE,
                     tip.length = 0.01, size = 1.93) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.07)),
                     breaks = seq(from = 0, to = 17000, by = 5000)) +
  scale_x_discrete(labels = c("Pre-COVID", "During\nCOVID")) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'))

cairo_pdf(paste0("finalPlots/Fig2/StepsCovid.Boxplot.pdf"), width = 1.4, height = 1.4, family = "Arial")
print(gSteps)
dev.off()

svglite(paste0("finalPlots/Fig2/StepsCovid.Boxplot.svg"), width = 1.4, height = 1.4, system_fonts = list(sans = "Arial"))
print(gSteps)
dev.off()

ggsave(paste0("finalPlots/Fig2/StepsCovid.Boxplot.tiff"), gSteps, width = 1.4, height = 1.4, units = "in", dpi=320)

covidAct <- covidAct[ID %in% covidAct[avgSleep > 120]$ID]

gSleep <- ggboxplot(data = covidAct, x = "covid", y = "avgSleep",
                    color = "covid", palette = c("darkorange", "purple"),
                    xlab = FALSE, ylab = "Avg Min Asleep / Day", ylim = c(100, 650),
                    legend = "none", add = "jitter", add.params = list(size = 0.01)) +
  stat_compare_means(comparisons = list(c("Pre-COVID", "Post-COVID")),
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = 600,
                     paired = TRUE,
                     tip.length = 0.01, size = 1.93) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.07)),
                     breaks = seq(from = 100, to = 650, by = 200)) +
  scale_x_discrete(labels = c("Pre-COVID", "During\nCOVID")) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'))

cairo_pdf(paste0("finalPlots/Fig2/SleepCovid.Boxplot.pdf"), width = 1.4, height = 1.4, family = "Arial")
print(gSleep)
dev.off()

svglite(paste0("finalPlots/Fig2/SleepCovid.Boxplot.svg"), width = 1.4, height = 1.4, system_fonts = list(sans = "Arial"))
print(gSleep)
dev.off()

ggsave(paste0("finalPlots/Fig2/SleepCovid.Boxplot.tiff"), gSleep, width = 1.4, height = 1.4, units = "in", dpi=320)

# get sleep and activity data over time.
dailyAct <- loadDailyActivityData(combinedData)

# keep only Phase 2
dailyAct <- dailyAct[Phase == "Phase 2"]

dailyAct$month <- factor(as.yearmon(dailyAct$Date), ordered = TRUE)
dailyAct$month <- factor(dailyAct$month, ordered = TRUE,
                         labels = paste0(substr(levels(dailyAct$month), 1, 3), "-", substr(levels(dailyAct$month), 7, 9)))

# now let's do it by month
monthlyData <- dailyAct[, list(avgSteps = mean(StepTotal, na.rm = TRUE),
                               missingSteps = sum(is.na(StepTotal)) / length(StepTotal),
                               avgSleep = mean(TotalMinutesAsleep, na.rm = TRUE),
                               missingSleep = sum(is.na(TotalMinutesAsleep)) / length(TotalMinutesAsleep),
                               avgTimeInBed = mean(TotalTimeInBed, na.rm = TRUE)),
                        by = c("ID", "Cohort", "month", "Sex", "Depression", "Anxiety")]

monthlyData <- monthlyData[month != "Sep-19"] # keep or remove Spetempber
monthlyInd <- monthlyData[, list(missingSteps = mean(missingSteps),
                                  missingSleep = mean(missingSleep)),
                           by = c("ID", "Sex", "Depression", "Anxiety")]

idMissingSteps <- monthlyInd[which(monthlyInd$missingSteps > 0.75)]$ID
idMissingSleep <- monthlyInd[which(monthlyInd$missingSleep > 0.75)]$ID

monthlyData[ID %in% idMissingSteps]$avgSteps <- NA
#monthlyData2[ID %in% idMissingSleep]$avgSleep <- NA

mActData <- droplevels(monthlyData)

td <- mActData[, lapply(.SD, mean, na.rm = TRUE),
               by = c("month", "Depression"), .SDcols = c("avgSteps", "avgSleep", "avgTimeInBed")]

td2 <- mActData[, lapply(.SD, std.error, na.rm = TRUE),
                by = c("month", "Depression"), .SDcols = c("avgSteps", "avgSleep", "avgTimeInBed")]

colnames(td2)[c(3, 4, 5)] <- c("stepsSE", "sleepSE", "timeBedSE")

td <- merge(td, td2, by = c("month", "Depression"))

gAct <- ggplot(td, aes(month, avgSteps, group = Depression)) +
  geom_errorbar(aes(ymin = avgSteps-stepsSE, ymax = avgSteps+stepsSE), width = 0.1, linewidth = 0.25) +
  geom_line(aes(color = Depression, linetype = Depression), linewidth = 0.4) +
  geom_point(aes(color = Depression, shape = Depression), show.legend = FALSE, size = 0.75) +
  ylab("Monthly Avg Steps / day") + xlab("") + labs(linetype="", color = "") + 
  scale_color_manual(values = c("red", "blue"), labels = c("Highest PHQ-9 < 10", "Highest PHQ-9 ≥ 10")) +
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("Highest PHQ-9 < 10", "Highest PHQ-9 ≥ 10")) +
  scale_shape_manual(values = c(16, 1), labels = c("Highest PHQ-9 < 10", "Highest PHQ-9 ≥ 10")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(2500, 14500),
                     breaks = c(5000, 7500, 10000, 12500)) +
  # annotate("text", x = 3, y = 10000, label = "*") +
  # annotate("text", x = 4, y = 11300, label = "*") +
  # annotate("text", x = 11, y = 9300, label = "*") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.25),
        text = element_text(size = 6, family = "Arial", face = "plain"),
        axis.text.x = element_text(angle = 30, vjust = 0.8,),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.text = element_text(size = 6, family = "Arial", face = "plain"),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'))

# ggsave(paste0("finalPlots/Fig2/Steps.Time.tiff"), gAct, width = 4.4, height = 2.8, units = "in", dpi=320)
# 
# cairo_pdf(paste0("finalPlots/Fig2/Steps.Time.pdf"), width = 4.4, height = 2.8, family = "Arial")
# print(gAct)
# dev.off()
# 
# svglite(paste0("finalPlots/Fig2/Steps.Time.svg"), width = 4.4, height = 2.8, system_fonts = list(sans = "Arial"))
# print(gAct)
# dev.off()

ggsave(paste0("finalPlots/Fig2/Steps.Time.2.tiff"), gAct, width = 6.5, height = 5.6, units = "cm", dpi=320)

cairo_pdf(paste0("finalPlots/Fig2/Steps.Time.2.pdf"), width = 2.56, height = 2.205, family = "Arial")
print(gAct)
dev.off()

svglite(paste0("finalPlots/Fig2/Steps.Time.2.svg"), width = 2.56, height = 2.205, system_fonts = list(sans = "Arial"))
print(gAct)
dev.off()



# now sleep
gSleep <- ggplot(td, aes(month, avgSleep, group = Depression)) +
  geom_errorbar(aes(ymin = avgSleep-sleepSE, ymax = avgSleep+sleepSE), width = 0.1, linewidth = 0.25) +
  geom_line(aes(color = Depression, linetype = Depression), linewidth = 0.4) +
  geom_point(aes(color = Depression, shape = Depression), show.legend = FALSE, size = 0.75) +
  ylab("Monthly Avg Min Asleep / day") + xlab("") + labs(linetype="", color = "") + 
  scale_color_manual(values = c("red", "blue"), labels = c("Highest PHQ-9 < 10", "Highest PHQ-9 ≥ 10")) +
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("Highest PHQ-9 < 10", "Highest PHQ-9 ≥ 10")) +
  scale_shape_manual(values = c(16, 1), labels = c("Highest PHQ-9 < 10", "Highest PHQ-9 ≥ 10")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(350, 550),
                     breaks = c(400, 450, 500)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.25),
        text = element_text(size = 6, family = "Arial", face = "plain"),
        axis.text.x = element_text(angle = 30, vjust = 0.8,),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.text = element_text(size = 6, family = "Arial", face = "plain"),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'))

# ggsave(paste0("finalPlots/Fig2/Sleep.Time.tiff"), gSleep, width = 4.4, height = 2.8, units = "in", dpi=320)
# 
# cairo_pdf(paste0("finalPlots/Fig2/Sleep.Time.pdf"), width = 4.4, height = 2.8, family = "Arial")
# print(gSleep)
# dev.off()
# 
# svglite(paste0("finalPlots/Fig2/Sleep.Time.svg"), width = 4.4, height = 2.8, system_fonts = list(sans = "Arial"))
# print(gSleep)
# dev.off()

ggsave(paste0("finalPlots/Fig2/Sleep.Time.2.tiff"), gSleep, width = 6.5, height = 5.6, units = "cm", dpi=320)

cairo_pdf(paste0("finalPlots/Fig2/Sleep.Time2..pdf"), width = 2.56, height = 2.205, family = "Arial")
print(gSleep)
dev.off()

svglite(paste0("finalPlots/Fig2/Sleep.Time.2.svg"), width = 2.56, height = 2.205, system_fonts = list(sans = "Arial"))
print(gSleep)
dev.off()

# Make the supplemental tables for the means & medians
p2PhqCovid$measure <- "PHQ9"
p2GadCovid$measure <- "GAD7"

# the activity covid act
act <- covidAct[, c("ID", "covid", "avgSteps")]

covidAct <- covidAct[ID %in% covidAct[avgSleep > 120]$ID]
sleep <- covidAct[, c("ID", "covid", "avgSleep")]

act$covid <- ifelse(act$covid == "Post-COVID", "During COVID", "Pre-COVID")
sleep$covid <- ifelse(sleep$covid == "Post-COVID", "During COVID", "Pre-COVID")

prePostCovid <- Reduce(function(...) merge(..., by = c("ID", "covid"), all = TRUE), list(p2PhqCovid, p2GadCovid, act, sleep))

vars <- c("maxPhq", "maxGad", "avgSteps", "avgSleep")
prePost <- c("Pre-COVID", "During COVID")

outList <- lapply(prePost, function(cv) {
  outList2 <- lapply(vars, function(v) {
    
    wd <- prePostCovid[covid == cv][[v]]
    
    data.frame(
      prePost = cv,
      measure = v,
      mean = mean(wd, na.rm = TRUE),
      median = median(wd, na.rm = TRUE),
      stdDev = sd(wd, na.rm = TRUE)
    )
  })
  rbindlist(outList2)
})

fwrite(rbindlist(outList), "temp/covidMeans.csv")

# get the significance at each timepoint
modAct <- lmer(avgSteps ~ Depression * month + (1|ID), data = mActData)

modActPre <- lmer(avgSteps ~ Depression * month + (1|ID), data = mActData[month < "Mar-20"])
modActPost <- lmer(avgSteps ~ Depression * month + (1|ID), data = mActData[month > "Jun-20"])

eModAct <- emmeans(modAct, ~ Depression | month)
contrast(eModAct, method = "pairwise", adjust = "tukey")

modSleep <- lmer(avgSleep ~ Depression * month + (1|ID), data = mActData)
eModSleep <- emmeans(modSleep, ~ Depression | month)
contrast(eModSleep, method = "pairwise", adjust = "tukey")



##########################################################################
# Figure 3
##########################################################################

phase_labeller <- as_labeller(c(`Phase 1` = "Phase 1         r = 0.25, p = 0.0098",
                                `Phase 2` = "Phase 2           r = 0.2, p = 0.035",
                                `Phase 3` = "Phase 3          r = 0.067, p = 0.49"))

g1 <- ggplot(combinedData, aes(PRS, PHQ9)) +
  geom_point(size = 0.3, shape = 16) +
  geom_smooth(method = "lm", linewidth = 0.3) +
  xlab("MDD-PRS") + ylab("Baseline PHQ-9") +
  scale_y_continuous(limits = c(-0.8, 25),
                     expand = expansion(mult = c(0, 0)),
                     breaks = c(0, 10, 20)) +
  scale_x_continuous(limits = c(-3, 3),
                     expand = expansion(mult = 0),
                     breaks = c(-2:2)) +
  facet_wrap("Phase", nrow = 1, labeller = phase_labeller) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial", face = "plain", size = 6),
    axis.ticks = element_line(linewidth = 0.2),
    strip.background = element_blank(),
    strip.text = element_text(margin = margin(0, 0, 1, 0, "pt"), hjust = 0,
                              family = "Arial", size = 6)
  )

ggsave("finalPlots/Fig3/BL.PHQ.PRS.tiff", g1, width = 4.3, height = 1.4, units = "in", dpi = 320)

cairo_pdf("finalPlots/Fig3/BL.PHQ.PRS.pdf", width = 4.3, height = 1.4, family = "Arial")
print(g1)
dev.off()

svglite("finalPlots/Fig3/BL.PHQ.PRS.svg", width = 4.3, height = 1.4, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

phase_labeller <- as_labeller(c(`Phase 1` = "Phase 1         r = 0.32, p = 0.0009",
                                `Phase 2` = "Phase 2          r = 0.094, p = 0.32",
                                `Phase 3` = "Phase 3          r = 0.059, p = 0.49"))

g1 <- ggplot(combinedData, aes(PRS, FU4TPPHQ9)) +
  geom_point(size = 0.3, shape = 16) +
  geom_smooth(method = "lm", linewidth = 0.3) +
  xlab("MDD-PRS") + ylab("Follow-up PHQ-9") +
  scale_y_continuous(limits = c(-0.8, 25),
                     expand = expansion(mult = c(0, 0)),
                     breaks = c(0, 10, 20)) +
  scale_x_continuous(limits = c(-3, 3),
                     expand = expansion(mult = 0),
                     breaks = c(-2:2)) +
  facet_wrap("Phase", nrow = 1, labeller = phase_labeller) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial", face = "plain", size = 6),
    axis.ticks = element_line(linewidth = 0.2),
    strip.background = element_blank(),
    strip.text = element_text(margin = margin(0, 0, 1, 0, "pt"), hjust = 0,
                              family = "Arial", size = 6)
  )

ggsave("finalPlots/Fig3/FU.PHQ.PRS.tiff", g1, width = 4.3, height = 1.4, units = "in", dpi = 320)

cairo_pdf("finalPlots/Fig3/FU.PHQ.PRS.pdf", width = 4.3, height = 1.4, family = "Arial")
print(g1)
dev.off()

svglite("finalPlots/Fig3/FU.PHQ.PRS.svg", width = 4.3, height = 1.4, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

##########################################################################
# Figure 4
##########################################################################

phase_labeller <- as_labeller(c(`Phase1` = "Phase 1",
                                `Phase2` = "Phase 2",
                                `Phase3` = "Phase 3"))

td <- cMood[, lapply(.SD, mean, na.rm = TRUE),
            by = c("Phase", "Event", "PRSHiLo", "Sex"), .SDcols = c("GAD7", "PHQ9")]

td2 <- cMood[, lapply(.SD, std.error, na.rm = TRUE),
             by = c("Phase", "Event", "PRSHiLo", "Sex"), .SDcols = c("GAD7", "PHQ9")]

colnames(td2)[c(5, 6)] <- c("gadSE", "phqSE")

td <- merge(td, td2, by = c("Phase", "Event", "PRSHiLo", "Sex"))

tStars <- fread("output/PHQ.GAD.Sex.Reg.Contrast.corrected.csv")
tStars$PRSHiLo <- "High PRS"
tStars$stars <- ifelse(tStars$p.value < 0.001, "***",
                       ifelse(tStars$p.value < 0.01, "**",
                              ifelse(tStars$p.value < 0.05, "*", NA)))

# Female PHQ-9
g1 <- ggplot(td[Sex == "Female"], aes(Event, PHQ9, group = PRSHiLo)) +
  geom_errorbar(aes(ymin = PHQ9-phqSE, ymax = PHQ9+phqSE), width = 0.1, linewidth = 0.25) +
  geom_line(aes(color = PRSHiLo, linetype = PRSHiLo), linewidth = 0.35) +
  geom_point(aes(color = PRSHiLo, shape = PRSHiLo), show.legend = FALSE, size = 0.75) +
  facet_wrap(c("Phase"), nrow = 1, labeller = phase_labeller) +
  ylab("PHQ-9") + xlab("") +
  scale_color_manual(values = c("purple", "darkgreen"), labels = c("MDD-PRS > 0", "MDD-PRS ≤ 0")) +
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("MDD-PRS > 0", "MDD-PRS ≤ 0")) +
  scale_shape_manual(values = c(16, 1)) +
  scale_x_discrete(labels = c("Base\nline", "Dec", "Mar", "Jun", "Sep")) +
  scale_y_continuous(expand = expansion(mult = 0), limits = c(0, 10),
                                        breaks = c(2.5, 5, 7.5)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.25),
        text = element_text(size = 7, family = "Arial", face = "plain"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'),
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(0, 0, 2, 0), size = 7))
  #geom_text(data = tStars[Sex == "Female" & variable == "PHQ9"], mapping = aes(x = x, y = y, label = stars))

ggsave("finalPlots/Fig4/F.PHQ.tiff", g1, width = 4.75, height = 1.4, units = "in", dpi = 320)

cairo_pdf(paste0("finalPlots/Fig4/F.PHQ.pdf"), width = 4.4, height = 1.4, family = "Arial")
print(g1)
dev.off()

svglite(paste0("finalPlots/Fig4/F.PHQ.svg"), width = 4.4, height = 1.4, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

# Male PhQ-9
g1 <- ggplot(td[Sex == "Male"], aes(Event, PHQ9, group = PRSHiLo)) +
  geom_errorbar(aes(ymin = PHQ9-phqSE, ymax = PHQ9+phqSE), width = 0.1, linewidth = 0.25) +
  geom_line(aes(color = PRSHiLo, linetype = PRSHiLo), linewidth = 0.35) +
  geom_point(aes(color = PRSHiLo, shape = PRSHiLo), show.legend = FALSE, size = 0.75) +
  facet_wrap(c("Phase"), nrow = 1, labeller = phase_labeller) +
  ylab("PHQ-9") + xlab("") +
  scale_color_manual(values = c("purple", "darkgreen"), labels = c("MDD-PRS > 0", "MDD-PRS ≤ 0")) +
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("MDD-PRS > 0", "MDD-PRS ≤ 0")) +
  scale_shape_manual(values = c(16, 1)) +
  scale_x_discrete(labels = c("Base\nline", "Dec", "Mar", "Jun", "Sep")) +
  scale_y_continuous(expand = expansion(mult = 0), limits = c(0, 10),
                     breaks = c(2.5, 5, 7.5)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.25),
        text = element_text(size = 7, family = "Arial", face = "plain"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'),
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(0, 0, 2, 0), size = 7))
  #geom_text(data = tStars[Sex == "Male" & variable == "PHQ9"], mapping = aes(x = x, y = y, label = stars))
  
ggsave("finalPlots/Fig4/M.PHQ.tiff", g1, width = 4.75, height = 1.4, units = "in", dpi = 320)

cairo_pdf(paste0("finalPlots/Fig4/M.PHQ.pdf"), width = 4.4, height = 1.4, family = "Arial")
print(g1)
dev.off()

svglite(paste0("finalPlots/Fig4/M.PHQ.svg"), width = 4.4, height = 1.4, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

# Female GAD7
g1 <- ggplot(td[Sex == "Female"], aes(Event, GAD7, group = PRSHiLo)) +
  geom_errorbar(aes(ymin = GAD7-gadSE, ymax = GAD7+gadSE), width = 0.1, linewidth = 0.25) +
  geom_line(aes(color = PRSHiLo, linetype = PRSHiLo), linewidth = 0.35) +
  geom_point(aes(color = PRSHiLo, shape = PRSHiLo), show.legend = FALSE, size = 0.75) +
  facet_wrap(c("Phase"), nrow = 1, labeller = phase_labeller) +
  ylab("GAD-7") + xlab("") +
  scale_color_manual(values = c("purple", "darkgreen"), labels = c("MDD-PRS > 0", "MDD-PRS ≤ 0")) +
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("MDD-PRS > 0", "MDD-PRS ≤ 0")) +
  scale_shape_manual(values = c(16, 1)) +
  scale_x_discrete(labels = c("Base\nline", "Dec", "Mar", "Jun", "Sep")) +
  scale_y_continuous(expand = expansion(mult = 0), limits = c(0, 10),
                     breaks = c(2.5, 5, 7.5)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.25),
        text = element_text(size = 7, family = "Arial", face = "plain"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'),
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(0, 0, 2, 0), size = 7))
  #geom_text(data = tStars[Sex == "Female" & variable == "GAD7"], mapping = aes(x = x, y = y, label = stars))

ggsave("finalPlots/Fig4/F.GAD.tiff", g1, width = 4.75, height = 1.4, units = "in", dpi = 320)

cairo_pdf(paste0("finalPlots/Fig4/F.GAD.pdf"), width = 4.4, height = 1.4, family = "Arial")
print(g1)
dev.off()

svglite(paste0("finalPlots/Fig4/F.GAD.svg"), width = 4.4, height = 1.4, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

# Male GAD7
g1 <- ggplot(td[Sex == "Male"], aes(Event, GAD7, group = PRSHiLo)) +
  geom_errorbar(aes(ymin = GAD7-gadSE, ymax = GAD7+gadSE), width = 0.1, linewidth = 0.25) +
  geom_line(aes(color = PRSHiLo, linetype = PRSHiLo), linewidth = 0.35) +
  geom_point(aes(color = PRSHiLo, shape = PRSHiLo), show.legend = FALSE, size = 0.75) +
  facet_wrap(c("Phase"), nrow = 1, labeller = phase_labeller) +
  ylab("GAD-7") + xlab("") +
  scale_color_manual(values = c("purple", "darkgreen"), labels = c("MDD-PRS > 0", "MDD-PRS ≤ 0")) +
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("MDD-PRS > 0", "MDD-PRS ≤ 0")) +
  scale_shape_manual(values = c(16, 1)) +
  scale_x_discrete(labels = c("Base\nline", "Dec", "Mar", "Jun", "Sep")) +
  scale_y_continuous(expand = expansion(mult = 0), limits = c(0, 10),
                     breaks = c(2.5, 5, 7.5)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.25),
        text = element_text(size = 7, family = "Arial", face = "plain"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'),
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(0, 0, 2, 0), size = 7))
  #geom_text(data = tStars[Sex == "Male" & variable == "GAD7"], mapping = aes(x = x, y = y, label = stars))

ggsave("finalPlots/Fig4/M.GAD.tiff", g1, width = 4.75, height = 1.4, units = "in", dpi = 320)

cairo_pdf(paste0("finalPlots/Fig4/M.GAD.pdf"), width = 4.4, height = 1.4, family = "Arial")
print(g1)
dev.off()

svglite(paste0("finalPlots/Fig4/M.GAD.svg"), width = 4.4, height = 1.4, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

##
# Now for Affect Score
phase_labeller <- as_labeller(c(`Phase1` = "Phase 1",
                                `Phase2` = "Phase 2",
                                `Phase3` = "Phase 3"))

cMood <- merge(cMood, cData[, c("ID", "bestZscore.1")], by = "ID")
cMood$ASHiLo <- ifelse(cMood$bestZscore.1 > 0, "High AS", "Low AS")

td <- cMood[, lapply(.SD, mean, na.rm = TRUE),
            by = c("Phase", "Event", "ASHiLo", "Sex"), .SDcols = c("GAD7", "PHQ9")]

td2 <- cMood[, lapply(.SD, std.error, na.rm = TRUE),
             by = c("Phase", "Event", "ASHiLo", "Sex"), .SDcols = c("GAD7", "PHQ9")]

colnames(td2)[c(5, 6)] <- c("gadSE", "phqSE")

td <- merge(td, td2, by = c("Phase", "Event", "ASHiLo", "Sex"))

tStars <- fread("output/PHQ.GAD.Sex.Reg.Contrast.corrected.csv")
tStars$PRSHiLo <- "High PRS"
tStars$stars <- ifelse(tStars$p.value < 0.001, "***",
                       ifelse(tStars$p.value < 0.01, "**",
                              ifelse(tStars$p.value < 0.05, "*", NA)))

# Female PHQ-9
g1 <- ggplot(td[Sex == "Female"], aes(Event, PHQ9, group = ASHiLo)) +
  geom_errorbar(aes(ymin = PHQ9-phqSE, ymax = PHQ9+phqSE), width = 0.1, linewidth = 0.25) +
  geom_line(aes(color = ASHiLo, linetype = ASHiLo), linewidth = 0.35) +
  geom_point(aes(color = ASHiLo, shape = ASHiLo), show.legend = FALSE, size = 0.75) +
  facet_wrap(c("Phase"), nrow = 1, labeller = phase_labeller) +
  ylab("PHQ-9") + xlab("") +
  scale_color_manual(values = c("purple", "darkgreen"), labels = c("AS > 0", "AS ≤ 0")) +
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("AS > 0", "AS ≤ 0")) +
  scale_shape_manual(values = c(16, 1)) +
  scale_x_discrete(labels = c("Base\nline", "Dec", "Mar", "Jun", "Sep")) +
  scale_y_continuous(expand = expansion(mult = 0), limits = c(0, 12),
                     breaks = c(2.5, 5, 7.5, 10)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.25),
        text = element_text(size = 7, family = "Arial", face = "plain"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'),
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(0, 0, 2, 0), size = 7))
#geom_text(data = tStars[Sex == "Female" & variable == "PHQ9"], mapping = aes(x = x, y = y, label = stars))

ggsave("finalPlots/Fig4/F.AS.PHQ.tiff", g1, width = 4.75, height = 1.4, units = "in", dpi = 320)


# Male PhQ-9
g1 <- ggplot(td[Sex == "Male"], aes(Event, PHQ9, group = ASHiLo)) +
  geom_errorbar(aes(ymin = PHQ9-phqSE, ymax = PHQ9+phqSE), width = 0.1, linewidth = 0.25) +
  geom_line(aes(color = ASHiLo, linetype = ASHiLo), linewidth = 0.35) +
  geom_point(aes(color = ASHiLo, shape = ASHiLo), show.legend = FALSE, size = 0.75) +
  facet_wrap(c("Phase"), nrow = 1, labeller = phase_labeller) +
  ylab("PHQ-9") + xlab("") +
  scale_color_manual(values = c("purple", "darkgreen"), labels = c("AS > 0", "AS ≤ 0")) +
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("AS > 0", "AS ≤ 0")) +
  scale_shape_manual(values = c(16, 1)) +
  scale_x_discrete(labels = c("Base\nline", "Dec", "Mar", "Jun", "Sep")) +
  scale_y_continuous(expand = expansion(mult = 0), limits = c(0, 12),
                     breaks = c(2.5, 5, 7.5, 10)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.25),
        text = element_text(size = 7, family = "Arial", face = "plain"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'),
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(0, 0, 2, 0), size = 7))
#geom_text(data = tStars[Sex == "Male" & variable == "PHQ9"], mapping = aes(x = x, y = y, label = stars))

ggsave("finalPlots/Fig4/M.AS.PHQ.tiff", g1, width = 4.75, height = 1.4, units = "in", dpi = 320)



##########################################################################
# Figure 5
##########################################################################

myComparisons <- list(c("Phase 1", "Phase 2"), c("Phase 2", "Phase 3"), c("Phase 1", "Phase 3"))

g1 <- ggboxplot(data = combinedData[Sex == "Female"], x = "Phase", y = "FU4TPPHQ9",
          color = "Phase", palette = "frontiers",
          xlab = FALSE, ylab = "Follow-up PHQ-9", ylim = c(1, 30),
          legend = "none", add = "jitter", add.params = list(size = 0.01),
          facet.by = "PRSHiLo", panel.labs = list(PRSHiLo = c("MDD-PRS > 0", "MDD-PRS ≤ 0"))) +
  stat_compare_means(comparisons = myComparisons,
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = c(21.75, 24, 26.25),
                     tip.length = 0.01, size = 2.1) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0)),
                     breaks = c(0, 10, 20)) +
  theme(text = element_text(size = 7, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'),
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(0,0,1,0), size = 7))

ggsave("finalPlots/Fig5/F.PHQ.tiff", g1, width = 3.54, height = 1.57, units = "in", dpi = 320)

cairo_pdf(paste0("finalPlots/Fig5/F.PHQ.pdf"), width = 3.54, height = 1.57, family = "Arial")
print(g1)
dev.off()

svglite(paste0("finalPlots/Fig5/F.PHQ.svg"), width = 3.54, height = 1.57, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

g1 <- ggboxplot(data = combinedData[Sex == "Male"], x = "Phase", y = "FU4TPPHQ9",
                color = "Phase", palette = "frontiers",
                xlab = FALSE, ylab = "Follow-up PHQ-9", ylim = c(1, 30),
                legend = "none", add = "jitter", add.params = list(size = 0.01),
                facet.by = "PRSHiLo", panel.labs = list(PRSHiLo = c("MDD-PRS > 0", "MDD-PRS ≤ 0"))) +
  stat_compare_means(comparisons = myComparisons,
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = c(21.75, 24, 26.25),
                     tip.length = 0.01, size = 2.1) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0)),
                     breaks = c(0, 10, 20)) +
  theme(text = element_text(size = 7, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'),
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(0,0,1,0), size = 7))

ggsave("finalPlots/Fig5/M.PHQ.tiff", g1, width = 3.54, height = 1.57, units = "in", dpi = 320)

cairo_pdf(paste0("finalPlots/Fig5/M.PHQ.pdf"), width = 3.54, height = 1.57, family = "Arial")
print(g1)
dev.off()

svglite(paste0("finalPlots/Fig5/M.PHQ.svg"), width = 3.54, height = 1.57, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

#Now GAD-7
g1 <- ggboxplot(data = combinedData[Sex == "Female"], x = "Phase", y = "FU4TPGAD7",
                color = "Phase", palette = "frontiers",
                xlab = FALSE, ylab = "Follow-up GAD-7", ylim = c(1, 30),
                legend = "none", add = "jitter", add.params = list(size = 0.01),
                facet.by = "PRSHiLo", panel.labs = list(PRSHiLo = c("MDD-PRS > 0", "MDD-PRS ≤ 0"))) +
  stat_compare_means(comparisons = myComparisons,
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = c(21.75, 24, 26.25),
                     tip.length = 0.01, size = 2.1) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0)),
                     breaks = c(0, 10, 20)) +
  theme(text = element_text(size = 7, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'),
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(0,0,1,0), size = 7))

ggsave("finalPlots/Fig5/F.GAD.tiff", g1, width = 3.54, height = 1.57, units = "in", dpi = 320)

cairo_pdf(paste0("finalPlots/Fig5/F.GAD.pdf"), width = 3.54, height = 1.57, family = "Arial")
print(g1)
dev.off()

svglite(paste0("finalPlots/Fig5/F.GAD.svg"), width = 3.54, height = 1.57, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()


g1 <- ggboxplot(data = combinedData[Sex == "Male"], x = "Phase", y = "FU4TPGAD7",
                color = "Phase", palette = "frontiers",
                xlab = FALSE, ylab = "Follow-up GAD-7", ylim = c(1, 30),
                legend = "none", add = "jitter", add.params = list(size = 0.01),
                facet.by = "PRSHiLo", panel.labs = list(PRSHiLo = c("MDD-PRS > 0", "MDD-PRS ≤ 0"))) +
  stat_compare_means(comparisons = myComparisons,
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = c(21.75, 24, 26.25),
                     tip.length = 0.01, size = 2.1) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0)),
                     breaks = c(0, 10, 20)) +
  theme(text = element_text(size = 7, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'),
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(0,0,1,0), size = 7))

ggsave("finalPlots/Fig5/M.GAD.tiff", g1, width = 3.54, height = 1.57, units = "in", dpi = 320)

cairo_pdf(paste0("finalPlots/Fig5/M.GAD.pdf"), width = 3.54, height = 1.57, family = "Arial")
print(g1)
dev.off()

svglite(paste0("finalPlots/Fig5/M.GAD.svg"), width = 3.54, height = 1.57, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()





##########################################################################
# Figure 6
##########################################################################

# the regression figs
# let's make the data.frame for the propo chart
varProps <- data.frame(pcOutBest$var$contrib)

varProps$Affect <- c("State", "Trait", "Family\nHistory", "State", "State", "Trait",
                      "State", "Family\nHistory", "State", "State")
varProps$subDim <- c("SBState", "SBTrait", "CTQ", "GAD-7", "PANSI+", "Neuroticism", "PSS", "RFQ", "PANSI-", "PHQ-9")

varProps$Affect <- factor(varProps$Affect, ordered = TRUE, levels = c("Trait", "State", "Family\nHistory"))

epitab(cData$FU4TPPHQCutoff, cData$bestZscore.1 > 0)

tiff("finalPlots/Fig6/Pie2.tiff", width = 2.441, height = 2.441, units = "in", res = 320)
PieDonut(varProps, aes(Affect, subDim, count = Dim.1), showPieName = FALSE, family = "Arial",
         donutLabelSize = 2.9, pieLabelSize = 2.9, labelposition = 1, ratioByGroup = FALSE)
dev.off()

cairo_pdf("finalPlots/Fig6/Pie2.pdf", width = 2.441, height = 2.441)
PieDonut(varProps, aes(Affect, subDim, count = Dim.1), showPieName = FALSE, family = "Arial",
         donutLabelSize = 2, pieLabelSize = 2, labelposition = 1, ratioByGroup = FALSE)
dev.off()

svglite(paste0("finalPlots/Fig6/Pie2.svg"), width = 2.441, height = 2.441)
PieDonut(varProps, aes(Affect, subDim, count = Dim.1), showPieName = FALSE, family = "Arial",
         donutLabelSize = 2, pieLabelSize = 2, labelposition = 1, ratioByGroup = FALSE)
dev.off()

g1 <- ggplot(cData, aes(x = dim.1, y = dim.2)) +
  geom_point(aes(shape = FU4TPPHQCutoff), size = 0.5) +
  geom_point(aes(color = FU4TPPHQCutoff), size = 0.5) +
  scale_color_manual(values = c("red", "blue"), labels = c("PHQ-9 < 10", "PHQ-9 ≥ 10")) +
  scale_shape_manual(values = c(16, 0), labels = c("PHQ-9 < 10", "PHQ-9 ≥ 10")) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8), limits = c(-5, 11),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = c(-2, 0, 2, 4),
                     expand = expansion(mult = c(0.025, 0.025))) +
  xlab("PCA Dimension 1 (58.55%)") + 
  ylab("PCA Dimension 2 (12.91%)") +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "longdash", linewidth = 0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(family = "Arial", size = 7),
        legend.spacing = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 1, 1, 1, unit = "pt"),
        legend.key.height = unit(0.075, "cm"),
        legend.key.width = unit(0.09, "cm"),
        legend.justification = "center",
        legend.text.align = 0.5,
        legend.background = element_blank(),
        legend.position = c(0.8672, 0.098),
        legend.box.background = element_rect(colour = "black", linewidth = 0.15),
        text = element_text(family = "Arial", size = 8),
        axis.ticks = element_line(linewidth = 0.25))

ggsave("finalPlots/Fig6/PCAPlot3.tiff", g1, width = 7.3, height = 4.4, units = "cm", dpi = 320)

cairo_pdf(paste0("finalPlots/Fig6/PCAPlot.pdf"), width = 3.75, height = 2, family = "Arial")
print(g1)
dev.off()

svglite(paste0("finalPlots/Fig6/PCAPlot.svg"), width = 3.75, height = 2, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

phase_labeller <- as_labeller(c(`Phase 1` = "Phase 1      r = 0.6, p = 2.2e-11",
                                `Phase 2` = "Phase 2      r = 0.71, p < 2.2e-16",
                                `Phase 3` = "Phase 3      r = 0.64, p = 8e-14"))

g1 <- ggplot(cData, aes(bestZscore.1, FU4TPPHQ9)) +
  geom_point(size = 0.3, shape = 16) +
  geom_smooth(method = "lm", linewidth = 0.3) +
  xlab("Affect Score") + ylab("Follow-up PHQ-9") +
  scale_y_continuous(limits = c(-0.8, 25),
                     expand = expansion(mult = c(0, 0)),
                     breaks = c(0, 10, 20)) +
  scale_x_continuous(limits = c(-2.5, 4.2),
                     expand = expansion(mult = 0),
                     breaks = c(-2, 0, 2, 4)) +
  facet_wrap("Phase", nrow = 1, labeller = phase_labeller) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial", face = "plain", size = 6),
    axis.ticks = element_line(linewidth = 0.2),
    strip.background = element_blank(),
    strip.text = element_text(margin = margin(0, 0, 1, 0, "pt"), hjust = 0,
                              family = "Arial", size = 6)
  )

ggsave("finalPlots/Fig6/AffectScore.PHQ.tiff", g1, width = 4.3, height = 1.4, units = "in", dpi = 320)

cairo_pdf("finalPlots/Fig6/AffectScore.PHQ.pdf", width = 4.3, height = 1.4, family = "Arial")
print(g1)
dev.off()

svglite("finalPlots/Fig6/AffectScore.PHQ.svg", width = 4.3, height = 1.4, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()


cData$affect <- ifelse(cData$bestZscore.1 > 0, "High Affect", "Low Affect")
cData$affect <- relevel(factor(cData$affect), ref = "Low Affect")

affect_labeller <- as_labeller(c(`High Affect` = "AS > 0        r = 0.083, p = 0.32",
                                 `Low Affect` = "AS ≤ 0        r = 0.16, p = 0.028"))

g1 <- ggplot(cData, aes(PRS, FU4TPPHQ9)) +
  geom_point(size = 0.3, shape = 16) +
  geom_smooth(method = "lm", linewidth = 0.3) +
  xlab("PRS") + ylab("Follow-up PHQ-9") +
  scale_y_continuous(limits = c(-0.8, 25),
                     expand = expansion(mult = c(0, 0)),
                     breaks = c(0, 10, 20)) +
  scale_x_continuous(limits = c(-3, 3),
                     expand = expansion(mult = 0),
                     breaks = c(-2, 0, 2)) +
  facet_wrap("affect", nrow = 1, labeller = affect_labeller) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial", face = "plain", size = 6),
    axis.ticks = element_line(linewidth = 0.2),
    strip.background = element_blank(),
    strip.text = element_text(margin = margin(0, 0, 1, 0, "pt"), hjust = 0,
                              family = "Arial", size = 6)
  )

ggsave("finalPlots/Fig6/PRS.AffectScore.PHQ.tiff", g1, width = 2.8667, height = 1.4, units = "in", dpi = 320)

cairo_pdf("finalPlots/Fig6/PRS.AffectScore.PHQ.pdf", width = 2.8667, height = 1.4, family = "Arial")
print(g1)
dev.off()

svglite("finalPlots/Fig6/PRS.AffectScore.PHQ.svg", width = 2.8667, height = 1.4, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

# the regression
cData$Sex <- relevel(factor(cData$Sex), ref = "Male")
bestMod <- lm(FU4TPPHQ9 ~ Sex + PRS + bestZscore.1 + Phase + Sex:Phase, data = cData)
summary(bestMod)

#minMod <- lm(FU4TPPHQ9 ~ Sex + PRS + minZscore + Phase + Sex:Phase, data = cData)

modelsummary(bestMod, statistic = c("Std.error: {std.error}", "P-value: {p.value}"),
             vcov = vcovHC(bestMod, "HC3"))

# Let's do the glm
glm1 <- glm(FU4TPPHQCutoff ~ PRS + Sex * Phase + bestZscore.1, data = cData, family = binomial)
predictedVals <- predict(glm1, type = "response")
g1 <- ggroc(roc(cData$FU4TPPHQCutoff, predictedVals, ci = TRUE)) +
  xlab("Specificity") + ylab("Sensitivity") +
  annotate("text", x = 0.85, y = 0.95, label = "AUROC = 0.835", size = 2.5) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    text = element_text(family = "Arial", size = 9)
  )

ggsave("finalPlots/SuppMisc/roc.tiff", g1, width = 4, height = 4, units = "in", dpi = 320)

cairo_pdf("finalPlots/SuppMisc/roc.pdf", width = 4, height = 4, family = "Arial")
print(g1)
dev.off()

svglite("finalPlots/SuppMisc/roc.svg", width = 4, height = 4, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

modelsummary(glm1, statistic = c("Std.error: {std.error}", "P-value: {p.value}"))

roc(cData$FU4TPPHQCutoff, predictedVals)

# add persistence
moodCombined <- moodData$combined[Event != "Sep BL"]

depTimes <- moodCombined[, list(depTimes = sum(PHQ9 >= 10, na.rm = TRUE),
                                persistence = sum(PHQ9 >= 10, na.rm = TRUE) > 1,
                                anyDep = any(PHQ9 >= 10, na.rm = TRUE)),
                         by = "ID"]

# combine it
tData <- merge(cData, depTimes, all.y = FALSE)

table(tData$persistence, tData$anyDep, tData$Sex, tData$Phase)
table(tData$persistence, tData$Sex, tData$Phase)
table(tData$anyDep, tData$Sex, tData$Phase)

tData$persistence <- ifelse(tData$persistence, "Chronic\nDep", "No Chronic\nDep")

g5 <- ggboxplot(data = tData, x = "persistence", y = "bestZscore.1",
                xlab = FALSE, ylab = "Affect Score", ylim = c(-3.5, 4.2),
                legend = "none", add = "jitter", add.params = list(size = 0.01)) +
  stat_compare_means(comparisons = list(c("No Chronic\nDep", "Chronic\nDep")),
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = 4.1,
                     tip.length = 0.01, size = 1.93) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     breaks = seq(from = -4, to = 4, by = 2)) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'))


ggsave(paste0("finalPlots/Fig6/AffectScore.persistence.Boxplot.tiff"), g5, width = 1.4, height = 1.4, units = "in", dpi=320)

tData$depCutoff <- ifelse(tData$FU4TPPHQCutoff == "Depression", "Depression", "No\nDepression")

# now add High PRS and AS
g6 <- ggboxplot(data = tData[PRS > 0], x = "depCutoff", y = "bestZscore.1",
                xlab = FALSE, ylab = "Affect Score", ylim = c(-3.5, 4.2),
                legend = "none", add = "jitter", add.params = list(size = 0.01)) +
  stat_compare_means(comparisons = list(c("No\nDepression", "Depression")),
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = 4.1,
                     tip.length = 0.01, size = 1.93) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     breaks = seq(from = -4, to = 4, by = 2)) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'))

ggsave(paste0("finalPlots/Fig6/AffectScore.HighPRS.Dep.Boxplot.tiff"), g6, width = 1.4, height = 1.4, units = "in", dpi=320)

g6 <- ggboxplot(data = tData[PRS <= 0], x = "depCutoff", y = "bestZscore.1",
                xlab = FALSE, ylab = "Affect Score", ylim = c(-3.5, 4.2),
                legend = "none", add = "jitter", add.params = list(size = 0.01)) +
  stat_compare_means(comparisons = list(c("No\nDepression", "Depression")),
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = 4.1,
                     tip.length = 0.01, size = 1.93) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     breaks = seq(from = -4, to = 4, by = 2)) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'))

ggsave(paste0("finalPlots/Fig6/AffectScore.LowPRS.Dep.Boxplot.tiff"), g6, width = 1.4, height = 1.4, units = "in", dpi=320)

g7f <- ggboxplot(data = tData[Sex == "Female"], x = "depCutoff", y = "bestZscore.1",
                xlab = FALSE, ylab = "Affect Score", ylim = c(-3.5, 4.2),
                legend = "none", add = "jitter", add.params = list(size = 0.01)) +
  stat_compare_means(comparisons = list(c("No\nDepression", "Depression")),
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = 4.1,
                     tip.length = 0.01, size = 1.93) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     breaks = seq(from = -4, to = 4, by = 2)) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'))

g7m <- ggboxplot(data = tData[Sex == "Male"], x = "depCutoff", y = "bestZscore.1",
                 xlab = FALSE, ylab = "Affect Score", ylim = c(-3.5, 4.2),
                 legend = "none", add = "jitter", add.params = list(size = 0.01)) +
  stat_compare_means(comparisons = list(c("No\nDepression", "Depression")),
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = 4.1,
                     tip.length = 0.01, size = 1.93) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     breaks = seq(from = -4, to = 4, by = 2)) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'))

ggsave(paste0("finalPlots/Fig6/AffectScore.F.Dep.Boxplot.tiff"), g7f, width = 1.4, height = 1.4, units = "in", dpi=320)
ggsave(paste0("finalPlots/Fig6/AffectScore.M.Dep.Boxplot.tiff"), g7m, width = 1.4, height = 1.4, units = "in", dpi=320)

# out <- stargazer(bestMod, header = FALSE, type = "html",
#           title = "Regressing the Follow-up PHQ-9 on Baseline variables",
#           label = "tab:depReg",
#           dep.var.labels = "Highest Follow-up PHQ-9 score",
#           star.cutoffs = c(0.05, 0.01, 0.001),
#           covariate.labels = c("Sex (Female)", "MDD-PRS", "Affect Dimension",
#                                "Phase (Phase 2)", "Phase (Phase 3)",
#                                "Sex(F) x Phase(2)", "Sex(F) x Phase(3)",
#                                "Intercept"))


##########################################################################
# Supplemental Figure
##########################################################################

phase_labeller <- as_labeller(c(`Phase 1` = "Phase 1            r = 0.18, p = 0.073",
                                `Phase 2` = "Phase 2            r = 0.13, p = 0.15",
                                `Phase 3` = "Phase 3            r = 0.0034, p = 0.97"))

g1 <- ggplot(combinedData, aes(PRS, GAD7)) +
  geom_point(size = 0.3, shape = 16) +
  geom_smooth(method = "lm", linewidth = 0.3) +
  xlab("MDD-PRS") + ylab("Baseline GAD-7") +
  scale_y_continuous(limits = c(-0.8, 25),
                     expand = expansion(mult = c(0, 0)),
                     breaks = c(0, 10, 20)) +
  scale_x_continuous(limits = c(-3, 3),
                     expand = expansion(mult = 0),
                     breaks = c(-2:2)) +
  facet_wrap("Phase", nrow = 1, labeller = phase_labeller) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial", face = "plain", size = 6),
    axis.ticks = element_line(linewidth = 0.2),
    strip.background = element_blank(),
    strip.text = element_text(margin = margin(0, 0, 1, 0, "pt"), hjust = 0)
  )

ggsave("finalPlots/Supp6/BL.GAD.PRS.tiff", g1, width = 4.3, height = 1.4, units = "in", dpi = 320)

cairo_pdf("finalPlots/Supp6/BL.GAD.PRS.pdf", width = 4.3, height = 1.4, family = "Arial")
print(g1)
dev.off()

svglite("finalPlots/Supp6/BL.GAD.PRS.svg", width = 4.3, height = 1.4, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

phase_labeller <- as_labeller(c(`Phase 1` = "Phase 1            r = 0.34, p = 0.00045",
                                `Phase 2` = "Phase 2            r = 0.16, p = 0.092",
                                `Phase 3` = "Phase 3            r = -0.0072, p = 0.94"))

g1 <- ggplot(combinedData, aes(PRS, FU4TPGAD7)) +
  geom_point(size = 0.3, shape = 16) +
  geom_smooth(method = "lm", linewidth = 0.3) +
  xlab("MDD-PRS") + ylab("Follow-up GAD-7") +
  scale_y_continuous(limits = c(-0.8, 25),
                     expand = expansion(mult = c(0, 0)),
                     breaks = c(0, 10, 20)) +
  scale_x_continuous(limits = c(-3, 3),
                     expand = expansion(mult = 0),
                     breaks = c(-2:2)) +
  facet_wrap("Phase", nrow = 1, labeller = phase_labeller) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial", face = "plain", size = 6),
    axis.ticks = element_line(linewidth = 0.2),
    strip.background = element_blank(),
    strip.text = element_text(margin = margin(0, 0, 1, 0, "pt"), hjust = 0)
  )

ggsave("finalPlots/Supp6/FU.GAD.PRS.tiff", g1, width = 4.3, height = 1.4, units = "in", dpi = 320)

cairo_pdf("finalPlots/Supp6/FU.GAD.PRS.pdf", width = 4.3, height = 1.4, family = "Arial")
print(g1)
dev.off()

svglite("finalPlots/Supp6/FU.GAD.PRS.svg", width = 4.3, height = 1.4, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()




##########################################################################
# Supplemental Figure 2
##########################################################################

# PHQ-9 over time
p1 <- fread("data/phase1/PHQ9.csv")
p2 <- fread("data/phase2/combined/phq.csv")
p3 <- fread("data/phase3/phq.csv")

p1$ID <- as.integer(substr(p1$USERID, 4, 7))

p1 <- p1[ID %in% combinedData$ID]
p2 <- p2[ID %in% combinedData$ID]
p3 <- p3[ID %in% combinedData$ID]

# remove screening
p2 <- p2[Event != "May 2019"]
p3 <- p3[Event != "May 2020"]

# set the events to correct
p1$time <- ifelse(p1$Event == "BL", "Baseline",
                  ifelse(p1$Event == "3M", "Dec",
                         ifelse(p1$Event == "6M", "Mar",
                                ifelse(p1$Event == "9M", "Jun", "Sep"))))
p1$time <- factor(p1$time, ordered = TRUE,
                  levels = c("Baseline", "Oct", "Nov", "Dec", "Jan", "Feb",
                             "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct.1", "Nov.1", "Dec.1"))

p2$time <- paste0(substr(p2$Event, 1, 3), "-", substr(p2$Event, 7, 8))
p2$time <- factor(p2$time, ordered = TRUE,
                  levels = c("Sep-19", "Oct-19", "Nov-19", "Dec-19", paste0(month.abb, "-20")))

p3$time <- paste0(substr(p3$Event, 1, 3), "-", substr(p3$Event, 7, 8))

p3$time <- factor(p3$time, ordered = TRUE,
                  levels = c("Sep-20", "Oct-20", "Nov-20", "Dec-20", paste0(month.abb, "-21")))

# add sex to the cols
p1 <- merge(p1, combinedData[, c("ID", "Sex")], by = "ID", all.x = TRUE)
p2 <- merge(p2, combinedData[, c("ID", "Sex")], by = "ID", all.x = TRUE)
p3 <- merge(p3, combinedData[, c("ID", "Sex")], by = "ID", all.x = TRUE)

# now get the graph lines
td <- p1[, lapply(.SD, mean, na.rm = TRUE),
         by = c("time", "Sex"), .SDcols = c("PHQ_SUM")]
colnames(td)[3] <- "phq"

td2 <- p1[, lapply(.SD, std.error, na.rm = TRUE),
          by = c("time", "Sex"), .SDcols = c("PHQ_SUM")]
colnames(td2)[3] <- "phqse"

td <- merge(td, td2, by = c("time", "Sex"))

g1 <- ggplot(td, aes(time, phq, group = Sex)) +
  geom_errorbar(aes(ymin = phq-phqse, ymax = phq+phqse), width = 0.1, linewidth = 0.25) +
  geom_line(aes(color = Sex, linetype = Sex), linewidth = 0.5) +
  geom_point(aes(color = Sex, shape = Sex), show.legend = FALSE, size = 0.75) +
  ylab("PHQ-9") + xlab("") + labs(linetype="", color = "") + 
  scale_color_manual(values = c("#F39B7FB2", "#4DBBD5B2"), labels = c("Female", "Male")) +
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("Female", "Male")) +
  scale_shape_manual(values = c(16, 1), labels = c("Female", "Male")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 12),
                     breaks = c(0, 5, 10)) +
  scale_x_discrete(limits = c("Baseline", "Oct", "Nov", "Dec", "Jan", "Feb",
                              "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep")) +
  # annotate("text", x = 7, y = 6, label = "*") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.3),
        panel.grid.major.y = element_line(linewidth = 0.3),
        axis.text.x = element_text(angle = 30, vjust = 0.8, size = 8),
        axis.text = element_text(size = 8),
        legend.position = c(0.9, 0.9),
        legend.background = element_blank(),
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.text = element_text(family = "Arial", face = "plain", size = 8),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'),
        text = element_text(size = 9, family = "Arial", face = "plain"))

ggsave(paste0("finalPlots/Supp2/phase1.tiff"), g1, width = 5.5, height = 3, units = "in", dpi=320)

cairo_pdf(paste0("finalPlots/Supp2/phase1.pdf"), width = 5.5, height = 3, family = "Arial")
print(g1)
dev.off()

svglite(paste0("finalPlots/Supp2/phase1.svg"), width = 5.5, height = 3, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

###
# Phase 2

td <- p2[, lapply(.SD, mean, na.rm = TRUE),
         by = c("time", "Sex"), .SDcols = c("phqSum")]
colnames(td)[3] <- "phq"

td2 <- p2[, lapply(.SD, std.error, na.rm = TRUE),
          by = c("time", "Sex"), .SDcols = c("phqSum")]
colnames(td2)[3] <- "phqse"

td <- merge(td, td2, by = c("time", "Sex"))

g1 <- ggplot(td, aes(time, phq, group = Sex)) +
  geom_errorbar(aes(ymin = phq-phqse, ymax = phq+phqse), width = 0.1, linewidth = 0.25) +
  geom_line(aes(color = Sex, linetype = Sex), linewidth = 0.5) +
  geom_point(aes(color = Sex, shape = Sex), show.legend = FALSE, size = 0.75) +
  ylab("PHQ-9") + xlab("") + labs(linetype="", color = "") + 
  scale_color_manual(values = c("#F39B7FB2", "#4DBBD5B2"), labels = c("Female", "Male")) +
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("Female", "Male")) +
  scale_shape_manual(values = c(16, 1), labels = c("Female", "Male")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 12),
                     breaks = c(0, 5, 10)) +
  scale_x_discrete(limits = c("Sep-19", "Oct-19", "Nov-19", "Dec-19", paste0(month.abb, "-20"))) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.3),
        panel.grid.major.y = element_line(linewidth = 0.3),
        axis.text.x = element_text(angle = 30, vjust = 0.8, size = 8),
        axis.text = element_text(size = 8),
        legend.position = c(0.9, 0.9),
        legend.background = element_blank(),
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.text = element_text(family = "Arial", face = "plain", size = 8),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'),
        text = element_text(size = 9, family = "Arial", face = "plain"))

ggsave(paste0("finalPlots/Supp2/phase2.tiff"), g1, width = 5.5, height = 3, units = "in", dpi=320)

cairo_pdf(paste0("finalPlots/Supp2/phase2.pdf"), width = 5.5, height = 3, family = "Arial")
print(g1)
dev.off()

svglite(paste0("finalPlots/Supp2/phase2.svg"), width = 5.5, height = 3, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

###
# Phase 3
td <- p3[, lapply(.SD, mean, na.rm = TRUE),
         by = c("time", "Sex"), .SDcols = c("phqSum")]
colnames(td)[3] <- "phq"

td2 <- p3[, lapply(.SD, std.error, na.rm = TRUE),
          by = c("time", "Sex"), .SDcols = c("phqSum")]
colnames(td2)[3] <- "phqse"

td <- merge(td, td2, by = c("time", "Sex"))

g1 <- ggplot(td, aes(time, phq, group = Sex)) +
  geom_errorbar(aes(ymin = phq-phqse, ymax = phq+phqse), width = 0.1, linewidth = 0.25) +
  geom_line(aes(color = Sex, linetype = Sex), linewidth = 0.5) +
  geom_point(aes(color = Sex, shape = Sex), show.legend = FALSE, size = 0.75) +
  ylab("PHQ-9") + xlab("") + labs(linetype="", color = "") + 
  scale_color_manual(values = c("#F39B7FB2", "#4DBBD5B2"), labels = c("Female", "Male")) +
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("Female", "Male")) +
  scale_shape_manual(values = c(16, 1), labels = c("Female", "Male")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 12),
                     breaks = c(0, 5, 10)) +
  scale_x_discrete(limits = c("Sep-20", "Oct-20", "Nov-20", "Dec-20", paste0(month.abb, "-21"))) +
  annotate("text", x = 1, y = 7, label = "*") +
  annotate("text", x = 2, y = 7.5, label = "*") +
  annotate("text", x = 3, y = 8.1, label = "*") +
  annotate("text", x = 4, y = 8, label = "**") +
  annotate("text", x = 5, y = 7.67, label = "**") +
  annotate("text", x = 6, y = 7.2, label = "*") +
  annotate("text", x = 7, y = 7.3, label = "*") +
  annotate("text", x = 8, y = 7.1, label = "*") +
  annotate("text", x = 9, y = 6.5, label = "*") +
  annotate("text", x = 10, y = 6.8, label = "*") +
  annotate("text", x = 11, y = 6.85, label = "**") +
  annotate("text", x = 12, y = 6.1, label = "*") +
  annotate("text", x = 13, y = 7, label = "*") +
  annotate("text", x = 14, y = 7.58, label = "*") +
  annotate("text", x = 15, y = 7.5, label = "*") +
  annotate("text", x = 16, y = 7, label = "*") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.3),
        panel.grid.major.y = element_line(linewidth = 0.3),
        axis.text.x = element_text(angle = 30, vjust = 0.8, size = 8),
        axis.text = element_text(size = 8),
        legend.position = c(0.9, 0.9),
        legend.background = element_blank(),
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.text = element_text(family = "Arial", face = "plain", size = 8),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'),
        text = element_text(size = 9, family = "Arial", face = "plain"))

ggsave(paste0("finalPlots/Supp2/phase3.tiff"), g1, width = 5.5, height = 3, units = "in", dpi=320)

cairo_pdf(paste0("finalPlots/Supp2/phase3.pdf"), width = 5.5, height = 3, family = "Arial")
print(g1)
dev.off()

svglite(paste0("finalPlots/Supp2/phase3.svg"), width = 5.5, height = 3, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

# get the regressions
m1 <- lmer(PHQ_SUM ~ Sex * time + (1|ID), data = p1)
m2 <- lmer(phqSum ~ Sex * time + (1|ID), data = p2)
m3 <- lmer(phqSum ~ Sex * time + (1|ID), data = p3)

pairs(emmeans(m1, (~ Sex | time)), simple = list("Sex"), combine = TRUE, adjust = "fdr")
pairs(emmeans(m2, (~ Sex | time)), simple = list("Sex"), combine = TRUE, adjust = "fdr")
pairs(emmeans(m3, (~ Sex | time)), simple = list("Sex"), combine = TRUE, adjust = "fdr")

anova(m1)
anova(m2)
anova(m3)

##########################################################################
# Supplemental Figure 3
##########################################################################

# let's get the averages for each id
dAvgData <- dActData[, list(avgSteps = mean(StepTotal, na.rm = TRUE),
                            avgSleep = mean(TotalMinutesAsleep, na.rm = TRUE)),
                     by = c("ID", "Phase", "Sex", "FU4TPPHQ9")]

phase_labeller <- as_labeller(c(`Phase 1` = "Phase 1          r = -0.2, p = 0.046",
                                `Phase 2` = "Phase 2          r = -0.13, p = 0.24",
                                `Phase 3` = "Phase 3          r = -0.22, p = 0.027"))

g1 <- ggplot(dAvgData[!is.na(avgSteps)], aes(avgSteps, FU4TPPHQ9)) +
  geom_point(size = 0.3, shape = 16) +
  geom_smooth(method = "lm", linewidth = 0.3) +
  xlab("Avg Steps / Day") + ylab("Follow-up PHQ-9") +
  scale_y_continuous(limits = c(0, 25),
                     expand = expansion(mult = c(0.05, 0.05)),
                     breaks = c(0, 10, 20)) +
  scale_x_continuous(limits = c(0, 18000),
                     expand = expansion(mult = c(0.05, 0.05)),
                     breaks = c(5000, 10000, 15000)) +
  facet_wrap("Phase", nrow = 1, labeller = phase_labeller) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial", face = "plain", size = 6),
    axis.ticks = element_line(linewidth = 0.2),
    strip.background = element_blank(),
    strip.text = element_text(margin = margin(0, 0, 1, 0, "pt"), hjust = 0,
                              family = "Arial", size = 6)
  )

ggsave("finalPlots/Supp3/steps.tiff", g1, width = 4.3, height = 1.4, units = "in", dpi = 320)

cairo_pdf("finalPlots/Supp3/steps.pdf", width = 4.3, height = 1.4, family = "Arial")
print(g1)
dev.off()

svglite("finalPlots/Supp3/steps.svg", width = 4.3, height = 1.4, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()


phase_labeller <- as_labeller(c(`Phase 1` = "Phase 1          r = 0.0025, p = 0.98",
                                `Phase 2` = "Phase 2          r = 0.046, p = 0.67",
                                `Phase 3` = "Phase 3          r = 0.18, p = 0.064"))

g1 <- ggplot(dAvgData[!is.na(avgSleep)], aes(avgSleep, FU4TPPHQ9)) +
  geom_point(size = 0.3, shape = 16) +
  geom_smooth(method = "lm", linewidth = 0.3) +
  xlab("Avg Min Asleep / Day") + ylab("Follow-up PHQ-9") +
  scale_y_continuous(limits = c(-3, 25),
                     expand = expansion(mult = c(0, 0.05)),
                     breaks = c(0, 10, 20)) +
  scale_x_continuous(limits = c(0, 800),
                     expand = expansion(mult = c(0.05, 0.05)),
                     breaks = c(200, 400, 600)) +
  facet_wrap("Phase", nrow = 1, labeller = phase_labeller) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial", face = "plain", size = 6),
    axis.ticks = element_line(linewidth = 0.2),
    strip.background = element_blank(),
    strip.text = element_text(margin = margin(0, 0, 1, 0, "pt"), hjust = 0,
                              family = "Arial", size = 6)
  )

ggsave("finalPlots/Supp3/sleep.tiff", g1, width = 4.3, height = 1.4, units = "in", dpi = 320)

cairo_pdf("finalPlots/Supp3/sleep.pdf", width = 4.3, height = 1.4, family = "Arial")
print(g1)
dev.off()

svglite("finalPlots/Supp3/sleep.svg", width = 4.3, height = 1.4, system_fonts = list(sans = "Arial"))
print(g1)
dev.off()

# ggscatter(dAvgData[!is.na(avgSteps)], x = "avgSteps", y = "FU4TPPHQ9", facet.by = "Phase",
#           add = "reg.line", conf.int = TRUE,
#           add.params = list(color = "blue", fill = "lightgray")) +
#   stat_cor(method = "pearson")
# 
# ggscatter(dAvgData[!is.na(avgSleep)], x = "avgSleep", y = "FU4TPPHQ9", facet.by = "Phase",
#           add = "reg.line", conf.int = TRUE,
#           add.params = list(color = "blue", fill = "lightgray")) +
#   stat_cor(method = "pearson")
# 


##########################################################################
# Supplemental Figure SCID
##########################################################################

# let's get the SCID data
scid <- fread("data/SCID.csv")

bestModSurveys <- c("SBState", "SBTrait", "CTQ", "GAD7", "PANSIPos", "NEOPIR", "PSS", "RFQ", "PANSINeg", "PHQ9")

pcOutBest <- PCA(imputePCA(combinedData[, ..bestModSurveys])$completeObs, graph = FALSE)

cData <- cbind(combinedData, data.frame(bestZscore.1 = scale(pcOutBest$ind$coord[, 1]),
                                        bestZscore.2 = scale(pcOutBest$ind$coord[, 2]),
                                        dim.1 = pcOutBest$ind$coord[, 1],
                                        dim.2 = pcOutBest$ind$coord[, 2]))

# now combine it with the data
tData <- merge(cData, scid, by.x = "ID", by.y = "id", all.x = TRUE, all.y = FALSE)

tData$pastAnxDep <- ifelse(tData$pastAnx | tData$pastDep, "Past Anxiety or Depression", "No History of Anxiety or Depression")
tData$pastAnx <- ifelse(tData$pastAnx, "Past\nAnxiety", "No Past\nAnxiety")
tData$currentAnx <- ifelse(tData$currentAnx, "Anxiety", "No\nAnxiety")
tData$pastDep <- ifelse(tData$pastDep, "Past\nDepression", "No Past\nDepression")
tData$currentDep <- ifelse(tData$currentDep, "Depression", "No\nDepression")
tData$pastCurrentAnx <- ifelse(tData$pastCurrentAnx, "Past or\nPresent\nAnxiety", "No Past or\nPresent\nAnxiety")
tData$pastCurrentDep <- ifelse(tData$pastCurrentDep, "Past or\nPresent\nDepression", "No Past or\nPresent\nDepression")

colsToGraph <- c("pastAnx", "currentAnx", "pastDep", "currentDep", "pastCurrentAnx", "pastCurrentDep")

lapply(colsToGraph, function(x) {
  
  myComparisons <- na.omit(unique(tData[[x]]))
  
  phq9ByPhase <- ggboxplot(data = tData[!is.na(tData[[x]])], x = x, y = "FU4TPPHQ9",
                           xlab = FALSE, ylab = "Highest Follow-up PHQ-9", ylim = c(-1, 25),
                           legend = "none", add = "jitter", add.params = list(size = 0.01),
                           facet.by = "Phase") +
    stat_compare_means(comparisons = list(myComparisons),
                       method = "wilcox.test",
                       label = c("p.format"), label.sep = " ",
                       hide.ns = FALSE,
                       label.y = 24,
                       tip.length = 0.01, size = 1.93) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                       breaks = seq(from = 0, to = 25, by = 10)) +
    theme(text = element_text(size = 6, family = "Arial", face = "plain"),
          #axis.text.x = element_text(angle = 30, vjust = 0.8),
          axis.line = element_line(linewidth = 0.3, color = "black"),
          axis.ticks = element_line(linewidth = 0.3, color = "black"),
          plot.margin = unit(c(4, 6, 4, 6), 'pt'))
  
  phq9 <- ggboxplot(data = tData[!is.na(tData[[x]])], x = x, y = "FU4TPPHQ9",
                    xlab = FALSE, ylab = "Highest Follow-up PHQ-9", ylim = c(-1, 25),
                    legend = "none", add = "jitter", add.params = list(size = 0.01)) +
    stat_compare_means(comparisons = list(myComparisons),
                       method = "wilcox.test",
                       label = c("p.format"), label.sep = " ",
                       hide.ns = FALSE,
                       label.y = 24,
                       tip.length = 0.01, size = 1.93) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                       breaks = seq(from = 0, to = 25, by = 10)) +
    theme(text = element_text(size = 6, family = "Arial", face = "plain"),
          #axis.text.x = element_text(angle = 30, vjust = 0.8),
          axis.line = element_line(linewidth = 0.3, color = "black"),
          axis.ticks = element_line(linewidth = 0.3, color = "black"),
          plot.margin = unit(c(4, 6, 4, 6), 'pt'))
  
  gad7ByPhase <- ggboxplot(data = tData[!is.na(tData[[x]])], x = x, y = "FU4TPGAD7",
                           xlab = FALSE, ylab = "Highest Follow-up GAD-7", ylim = c(-1, 25),
                           legend = "none", add = "jitter", add.params = list(size = 0.01),
                           facet.by = "Phase") +
    stat_compare_means(comparisons = list(myComparisons),
                       method = "wilcox.test",
                       label = c("p.format"), label.sep = " ",
                       hide.ns = FALSE,
                       label.y = 24,
                       tip.length = 0.01, size = 1.93) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                       breaks = seq(from = 0, to = 25, by = 10)) +
    theme(text = element_text(size = 6, family = "Arial", face = "plain"),
          #axis.text.x = element_text(angle = 30, vjust = 0.8),
          axis.line = element_line(linewidth = 0.3, color = "black"),
          axis.ticks = element_line(linewidth = 0.3, color = "black"),
          plot.margin = unit(c(4, 6, 4, 6), 'pt'))
  
  gad7 <- ggboxplot(data = tData[!is.na(tData[[x]])], x = x, y = "FU4TPGAD7",
                    xlab = FALSE, ylab = "Highest Follow-up GAD-7", ylim = c(-1, 25),
                    legend = "none", add = "jitter", add.params = list(size = 0.01)) +
    stat_compare_means(comparisons = list(myComparisons),
                       method = "wilcox.test",
                       label = c("p.format"), label.sep = " ",
                       hide.ns = FALSE,
                       label.y = 24,
                       tip.length = 0.01, size = 1.93) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                       breaks = seq(from = 0, to = 25, by = 10)) +
    theme(text = element_text(size = 6, family = "Arial", face = "plain"),
          #axis.text.x = element_text(angle = 30, vjust = 0.8),
          axis.line = element_line(linewidth = 0.3, color = "black"),
          axis.ticks = element_line(linewidth = 0.3, color = "black"),
          plot.margin = unit(c(4, 6, 4, 6), 'pt'))
  
  # save these
  ggsave(paste0("finalPlots/SuppSCID/FUPHQ9.Phase.", x, ".Boxplot.tiff"), phq9ByPhase, width = 4.2, height = 1.4, units = "in", dpi=320)
  ggsave(paste0("finalPlots/SuppSCID/FUPHQ9", x, ".Boxplot.tiff"), phq9, width = 1.4, height = 1.4, units = "in", dpi=320)
  ggsave(paste0("finalPlots/SuppSCID/FUGAD7.Phase.", x, ".Boxplot.tiff"), gad7ByPhase, width = 4.2, height = 1.4, units = "in", dpi=320)
  ggsave(paste0("finalPlots/SuppSCID/FUGAD7", x, ".Boxplot.tiff"), gad7, width = 1.4, height = 1.4, units = "in", dpi=320)
  
  return(NULL)
})

myComparisons <- na.omit(unique(tData$pastDep))
# do a graph of PRS and past depression
prs <- ggboxplot(data = tData[!is.na(pastDep)], x = "pastDep", y = "PRS",
                 xlab = FALSE, ylab = "MDD-PRS", ylim = c(-3.5, 4),
                 legend = "none", add = "jitter", add.params = list(size = 0.01)) +
  stat_compare_means(comparisons = list(myComparisons),
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = 3.75,
                     tip.length = 0.01, size = 1.93) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     breaks = seq(from = -4, to = 4, by = 2)) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'))

prsByPhase <- ggboxplot(data = tData[!is.na(pastDep)], x = "pastDep", y = "PRS",
                        xlab = FALSE, ylab = "MDD-PRS", ylim = c(-3.5, 4),
                        legend = "none", add = "jitter", add.params = list(size = 0.01),
                        facet.by = "Phase") +
  stat_compare_means(comparisons = list(myComparisons),
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = 3.75,
                     tip.length = 0.01, size = 1.93) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     breaks = seq(from = -4, to = 4, by = 2)) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'))

ggsave(paste0("finalPlots/SuppSCID/PRS.Phase.pastDep.Boxplot.tiff"), prsByPhase, width = 4.2, height = 1.4, units = "in", dpi=320)
ggsave(paste0("finalPlots/SuppSCID/PRS.pastDep.Boxplot.tiff"), prs, width = 1.4, height = 1.4, units = "in", dpi=320)

as <- ggboxplot(data = tData[!is.na(pastDep)], x = "pastDep", y = "bestZscore.1",
                 xlab = FALSE, ylab = "Affect Score", ylim = c(-3.5, 4),
                 legend = "none", add = "jitter", add.params = list(size = 0.01)) +
  stat_compare_means(comparisons = list(myComparisons),
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = 3.75,
                     tip.length = 0.01, size = 1.93) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     breaks = seq(from = -4, to = 4, by = 2)) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'))

asByPhase <- ggboxplot(data = tData[!is.na(pastDep)], x = "pastDep", y = "bestZscore.1",
                        xlab = FALSE, ylab = "Affect SCore", ylim = c(-3.5, 4),
                        legend = "none", add = "jitter", add.params = list(size = 0.01),
                        facet.by = "Phase") +
  stat_compare_means(comparisons = list(myComparisons),
                     method = "wilcox.test",
                     label = c("p.format"), label.sep = " ",
                     hide.ns = FALSE,
                     label.y = 3.75,
                     tip.length = 0.01, size = 1.93) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     breaks = seq(from = -4, to = 4, by = 2)) +
  theme(text = element_text(size = 6, family = "Arial", face = "plain"),
        #axis.text.x = element_text(angle = 30, vjust = 0.8),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.3, color = "black"),
        plot.margin = unit(c(4, 6, 4, 6), 'pt'))

ggsave(paste0("finalPlots/SuppSCID/AffectScore.Phase.pastDep.Boxplot.tiff"), asByPhase, width = 4.2, height = 1.4, units = "in", dpi=320)
ggsave(paste0("finalPlots/SuppSCID/AffectScore.pastDep.Boxplot.tiff"), as, width = 1.4, height = 1.4, units = "in", dpi=320)

# add the SCID to the regression
tData$Sex <- relevel(factor(tData$Sex), ref = "Male")
bestMod <- lm(FU4TPPHQ9 ~ Sex + PRS + bestZscore.1 + Phase + Sex:Phase, data = tData[!is.na(pastAnxDep)])
summary(bestMod)

bestMod2 <- lm(FU4TPPHQ9 ~ Sex + PRS + bestZscore.1 + Phase + Sex:Phase + pastAnxDep, data = tData)

bestMod3 <- lm(FU4TPPHQ9 ~ Sex + bestZscore.1 + Phase + Sex:Phase + pastAnxDep, data = tData)

glm1 <- glm(FU4TPPHQCutoff ~ PRS + Sex * Phase + bestZscore.1, data = tData[!is.na(pastAnxDep)], family = binomial)
glm2 <- glm(FU4TPPHQCutoff ~ PRS + Sex * Phase + bestZscore.1 + pastAnxDep, data = tData[!is.na(pastAnxDep)], family = binomial)

predictedVals1 <- predict(glm1, type = "response")
predictedVals2 <- predict(glm2, type = "response")

roc(tData[!is.na(pastAnxDep)]$FU4TPPHQCutoff, predictedVals1)
roc(tData[!is.na(pastAnxDep)]$FU4TPPHQCutoff, predictedVals2)



##########################################################################
# Supplemental Figure: PHQ-9 Sex differences over time
##########################################################################

# first Phase 1

x1 <- moodData[[1]]

# Combine mood and sex
x1 <- merge(x1, combinedData[, c("ID", "Sex")], by = "ID", all.y = FALSE)

td <- x1[, list(meanPHQ = mean(PHQ9, na.rm = TRUE),
               sePHQ = std.error(PHQ9, na.rm = TRUE)),
        by = c("Event", "Sex")]

g1 <- ggplot(td, aes(Event, meanPHQ, group = Sex)) +
  geom_errorbar(aes(ymin = meanPHQ-sePHQ, ymax = meanPHQ+sePHQ), width = 0.05, linewidth = 0.25) +
  geom_line(aes(color = Sex, linetype = Sex), linewidth = 0.4) +
  geom_point(aes(color = Sex, shape = Sex), show.legend = FALSE, size = 0.75) +
  ylab("PHQ-9") + xlab("") + labs(linetype="", color = "") + 
  scale_color_manual(values = c("orange", "cyan"), labels = c("Female", "Male")) +
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("Female", "Male")) +
  scale_shape_manual(values = c(16, 1), labels = c("Female", "Male")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 9),
                     breaks = c(2, 4, 6, 8)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.25),
        text = element_text(size = 6, family = "Arial", face = "plain"),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.text = element_text(size = 6, family = "Arial", face = "plain"),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'))

# now phase 2
x2 <- moodData[[2]]
x2 <- merge(x2, combinedData[, c("ID", "Sex")], by = "ID", all.y = FALSE)

td <- x2[, list(meanPHQ = mean(PHQ9, na.rm = TRUE),
               sePHQ = std.error(PHQ9, na.rm = TRUE)),
        by = c("Event", "Sex")]

td$Event <- my(td$Event)

g2 <- ggplot(td, aes(Event, meanPHQ, group = Sex)) +
  geom_errorbar(aes(ymin = meanPHQ-sePHQ, ymax = meanPHQ+sePHQ), width = 4, linewidth = 0.25) +
  geom_line(aes(color = Sex, linetype = Sex), linewidth = 0.4) +
  geom_point(aes(color = Sex, shape = Sex), show.legend = FALSE, size = 0.75) +
  ylab("PHQ-9") + xlab("") + labs(linetype="", color = "") + 
  scale_color_manual(values = c("orange", "cyan"), labels = c("Female", "Male")) +
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("Female", "Male")) +
  scale_shape_manual(values = c(16, 1), labels = c("Female", "Male")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 9),
                     breaks = c(2, 4, 6, 8)) +
  scale_x_date(date_labels = "%b-%y", breaks = breaks_pretty(n = 14)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.25),
        text = element_text(size = 6, family = "Arial", face = "plain"),
        axis.text.x = element_text(angle = 30, vjust = 0.8,),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.text = element_text(size = 6, family = "Arial", face = "plain"),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'))

# now phase 3
x3 <- moodData[[3]]
x3 <- merge(x3, combinedData[, c("ID", "Sex")], by = "ID", all.y = FALSE)

td <- x3[, list(meanPHQ = mean(PHQ9, na.rm = TRUE),
               sePHQ = std.error(PHQ9, na.rm = TRUE)),
        by = c("Event", "Sex")]

td$Event <- my(td$Event)

g3 <- ggplot(td, aes(Event, meanPHQ, group = Sex)) +
  geom_errorbar(aes(ymin = meanPHQ-sePHQ, ymax = meanPHQ+sePHQ), width = 4, linewidth = 0.25) +
  geom_line(aes(color = Sex, linetype = Sex), linewidth = 0.4) +
  geom_point(aes(color = Sex, shape = Sex), show.legend = FALSE, size = 0.75) +
  ylab("PHQ-9") + xlab("") + labs(linetype="", color = "") + 
  scale_color_manual(values = c("orange", "cyan"), labels = c("Female", "Male")) +
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("Female", "Male")) +
  scale_shape_manual(values = c(16, 1), labels = c("Female", "Male")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 9),
                     breaks = c(2, 4, 6, 8)) +
  scale_x_date(date_labels = "%b-%y", breaks = breaks_pretty(n = 14)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.25),
        text = element_text(size = 6, family = "Arial", face = "plain"),
        axis.text.x = element_text(angle = 30, vjust = 0.8,),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.text = element_text(size = 6, family = "Arial", face = "plain"),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = unit(c(4, 6, 0, 6), 'pt'))

# save the graphs
ggsave(paste0("finalPlots/SuppPHQSexDiff/Phase1.PHQ.SexDiff.tiff"), g1, width = 4.3, height = 3.7, units = "cm", dpi=320)
ggsave(paste0("finalPlots/SuppPHQSexDiff/Phase2.PHQ.SexDiff.tiff"), g2, width = 6.5, height = 5.6, units = "cm", dpi=320)
ggsave(paste0("finalPlots/SuppPHQSexDiff/Phase3.PHQ.SexDiff.tiff"), g3, width = 6.5, height = 5.6, units = "cm", dpi=320)

# stats

# get the significance at each timepoint
m1 <- lmer(PHQ9 ~ Sex * Event + (1|ID), data = x1)
m2 <- lmer(PHQ9 ~ Sex * Event + (1|ID), data = x2)
m3 <- lmer(PHQ9 ~ Sex * Event + (1|ID), data = x3)

contrast(emmeans(m1, ~ Sex | Event), method = "pairwise", adjust = "tukey")
contrast(emmeans(m2, ~ Sex | Event), method = "pairwise", adjust = "tukey")
contrast(emmeans(m3, ~ Sex | Event), method = "pairwise", adjust = "tukey")

anova(m1)
anova(m2)
anova(m3)




