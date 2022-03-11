
###########################################################################
# Load Libraries & Data  ---------------------------------------------
###########################################################################

library(readr)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
library(flextable)
library(gridExtra)
library(survival)
library(survminer)
library(RColorBrewer)
library(readr)

# set path
my_path <- paste0(getwd(), "/")

# read in functions
source(paste0(my_path, "Isolation_Fxns_031122.R"))

# set path to save figures
my_path_figures <- paste0(my_path, "Figures/")

# load data
dat <- read_csv(paste0(my_path, "dat_summary.csv"))

# select colors for plotting
all_colors <- brewer.pal(8, "Dark2")
# display.brewer.pal(8, "Dark2") 
customPalette <- c(all_colors[3], all_colors[1], all_colors[4]) # Cat1 = purple; Cat2 = green; Cat3 = pink
customPalette2 <- c(all_colors[8], all_colors[2], all_colors[6])

###########################################################################
# Table S1 - Vaccine Details  ---------------------------------------------
###########################################################################

# fully vaccinated w/ primary series >= 14 days before dx?
dat$Primary_vx_14d <- rep(NA, nrow(dat))
dat[which(dat$time_since_primary_vx_complete >= 14), "Primary_vx_14d"] <- "Yes"
dat[which(dat$time_since_primary_vx_complete < 14), "Primary_vx_14d"] <- "No"

# fully boosted >= 14 days before dx?
dat$Booster_vx_14d <- rep(NA, nrow(dat))
dat[which(dat$time_since_boost_vx >= 14), "Booster_vx_14d"] <- "Yes"
dat[which(dat$time_since_boost_vx < 14), "Booster_vx_14d"] <- "No"
dat[which((dat$Booster_vx_brand == "Unknown") & (is.na(dat$time_since_boost_vx) == TRUE)), "Booster_vx_14d"] <- "Unknown"

# 2nd booster >= 14 days before dx?
dat$Extra_booster_vx_14d <- rep(NA, nrow(dat))
dat[which(dat$time_since_2nd_boost >= 14), "Extra_booster_vx_14d"] <- "Yes"
dat[which(dat$time_since_2nd_boost < 14), "Extra_booster_vx_14d"] <- "No"
dat[which(dat$Extra_booster_brand == "Unknown"), "Extra_booster_vx_14d"] <- "Unknown" # for UID = 144

# create vaccination history table
dat_vx <- ddply(dat, .(Prior_primary_series, Primary_vx_brand_all, Booster_vx_brand, Extra_booster_brand, Primary_vx_14d, Booster_vx_14d, Extra_booster_vx_14d, Num_vx_doses), nrow)
dat_vx <- dat_vx[order(dat_vx$V1, decreasing = TRUE), ]
dat_vx[is.na(dat_vx)] <- "-"

tab_vx <- flextable(dat_vx)
tab_vx_format <- tab_vx %>% 
  add_header_row(top = TRUE, values = c("Vaccine Brands", 
                                        "Given >= 14 Days Before Diagnosis?", 
                                        "No. Doses", 
                                        "No. Persons"), colwidths = c(4, 3, 1, 1)) %>% 
  set_header_labels(       
    Prior_primary_series = "Prior Primary Series", 
    Primary_vx_brand_all = "Primary Series",                  
    Booster_vx_brand = "Booster",
    Extra_booster_brand = "Second Booster",
    Primary_vx_14d = "Primary Series",
    Booster_vx_14d = "Booster",
    Extra_booster_vx_14d = "Second Booster",
    Num_vx_doses = "",
    V1 = "") %>%
  align(align = "center", part = "header") %>%
  align(align = "center", part = "body") %>%
  bold(i = 1, bold = TRUE, part = "header") %>%
  italic(i = 2, italic = TRUE, part = "all") %>%
  vline(part = "all", j = c(4, 7)) %>% 
  bg(i = c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31), part = "body", bg = "gray95")  
tab_vx_format <- fix_border_issues(tab_vx_format)
tab_vx_format <- fontsize(tab_vx_format, size = 5, part = "all")

save_as_docx(tab_vx_format, path = paste0(my_path_figures, "Table_S1.docx"))

###########################################################################
# Table 1 - Sample Characteristics  ---------------------------------------
###########################################################################

keep_cat <- NULL

# symptom status
unique(dat$Sympt_status)
dat_cat <- ddply(dat, .(time_since_last_neg_cat, Sympt_status), nrow)
dat_cat <- reshape2::dcast(dat_cat, Sympt_status ~ time_since_last_neg_cat, value.var = "V1")
dat_cat$Sympt_status <- as.character(dat_cat$Sympt_status)
dat_cat[which(is.na(dat_cat$Sympt_status) == TRUE), "Sympt_status"] <- "Unknown"
dat_cat[is.na(dat_cat)] <- 0
cat_totals <- colSums(dat_cat[, 2:ncol(dat_cat)])
dat_cat$total <- rowSums(dat_cat[, 2:ncol(dat_cat)])
dat_cat$Perc_Cat1 <- dat_cat$Cat1 / cat_totals[1]
dat_cat$Perc_Cat2 <- dat_cat$Cat2 / cat_totals[2]
dat_cat$Perc_Cat3 <- dat_cat$Cat3 / cat_totals[3]
dat_cat$Perc_Unknown <- dat_cat$Unknown / cat_totals[4]
dat_cat$Perc_Total <- dat_cat$total / sum(dat_cat$total)
colnames(dat_cat) <- c("cat", 
                       "num_cat1", "num_cat2", "num_cat3", "num_catUNK", "num_total",
                       "perc_cat1", "perc_cat2", "perc_cat3", "perc_catUNK", "perc_total")
keep_cat <- rbind.data.frame(keep_cat, dat_cat)

# prior inf > 90 days
dat_cat <- ddply(dat, .(time_since_last_neg_cat, Prior_inf_more90d), nrow)
dat_cat <- reshape2::dcast(dat_cat, Prior_inf_more90d ~ time_since_last_neg_cat, value.var = "V1")
dat_cat[is.na(dat_cat)] <- 0
cat_totals <- colSums(dat_cat[, 2:ncol(dat_cat)])
dat_cat$total <- rowSums(dat_cat[, 2:ncol(dat_cat)])
dat_cat$Perc_Cat1 <- dat_cat$Cat1 / cat_totals[1]
dat_cat$Perc_Cat2 <- dat_cat$Cat2 / cat_totals[2]
dat_cat$Perc_Cat3 <- dat_cat$Cat3 / cat_totals[3]
dat_cat$Perc_Unknown <- dat_cat$Unknown / cat_totals[4]
dat_cat$Perc_Total <- dat_cat$total / sum(dat_cat$total)
colnames(dat_cat) <- c("cat", 
                       "num_cat1", "num_cat2", "num_cat3", "num_catUNK", "num_total",
                       "perc_cat1", "perc_cat2", "perc_cat3", "perc_catUNK", "perc_total")
keep_cat <- rbind.data.frame(keep_cat, dat_cat)

# num vx doses
dat_cat <- ddply(dat, .(time_since_last_neg_cat, Num_vx_doses), nrow)
dat_cat <- reshape2::dcast(dat_cat, Num_vx_doses ~ time_since_last_neg_cat, value.var = "V1")
dat_cat[is.na(dat_cat)] <- 0
cat_totals <- colSums(dat_cat[, 2:ncol(dat_cat)])
dat_cat$total <- rowSums(dat_cat[, 2:ncol(dat_cat)])
dat_cat$Perc_Cat1 <- dat_cat$Cat1 / cat_totals[1]
dat_cat$Perc_Cat2 <- dat_cat$Cat2 / cat_totals[2]
dat_cat$Perc_Cat3 <- dat_cat$Cat3 / cat_totals[3]
dat_cat$Perc_Unknown <- dat_cat$Unknown / cat_totals[4]
dat_cat$Perc_Total <- dat_cat$total / sum(dat_cat$total)
colnames(dat_cat) <- c("cat", 
                       "num_cat1", "num_cat2", "num_cat3", "num_catUNK", "num_total",
                       "perc_cat1", "perc_cat2", "perc_cat3", "perc_catUNK", "perc_total")
keep_cat <- rbind.data.frame(keep_cat, dat_cat)

keep_cat_names <- c(rep("Self-Reported Symptoms Prior to or At Initial Diagnosis", 3), 
                    rep("Prior Infection > 90 Days", 2), 
                    rep("Number of Vaccine Doses", 5))
keep_cat <- cbind.data.frame(keep_cat_names, keep_cat)

# compile into table
library(officer)
big_border = fp_border(color="black", width = 1)


keep_cat$cat <- c("No", "Yes", "Unknown", "No", "Yes", 1, 2, 3, 4, "Unknown")
colnames(keep_cat) <- c("Sample Characteristics", " ", 
                        "Number: <= 4", "Number: 5-9", "Number: >= 10", "Number: Unknown", "Number: Total",
                        "Percent: <=4", "Percent: 5-9", "Percent: >= 10", "Percent: Unknown", "Percent: Total")
keep_cat[, 8:ncol(keep_cat)] <- round(keep_cat[, 8:ncol(keep_cat)], digits = 2) * 100 # WARNING: if change cols, change #
keep_cat <- keep_cat[,  c("Sample Characteristics", " ", 
                          "Number: <= 4", "Percent: <=4",
                          "Number: 5-9", "Percent: 5-9",
                          "Number: >= 10", "Percent: >= 10",
                          "Number: Unknown", "Percent: Unknown",
                          "Number: Total", "Percent: Total")]
tab_cat <- cbind.data.frame(keep_cat$`Sample Characteristics`, keep_cat$` `, 
                            paste0(keep_cat$`Number: <= 4`, " (", keep_cat$`Percent: <=4`, ")"),
                            paste0(keep_cat$`Number: 5-9`, " (", keep_cat$`Percent: 5-9`, ")"),
                            paste0(keep_cat$`Number: >= 10`, " (", keep_cat$`Percent: >= 10`, ")"),
                            paste0(keep_cat$`Number: Unknown`, " (", keep_cat$`Percent: Unknown`, ")"),
                            paste0(keep_cat$`Number: Total`, " (", keep_cat$`Percent: Total`, ")"))

cat_pops <- ddply(dat, .(time_since_last_neg_cat), nrow)
cat_1_name <- paste0("<=4", " (N=", cat_pops[which(cat_pops$time_since_last_neg_cat == "Cat1"), "V1"], ")")
cat_2_name <- paste0("5-9", " (N=", cat_pops[which(cat_pops$time_since_last_neg_cat == "Cat2"), "V1"], ")")
cat_3_name <- paste0(">= 10", " (N=", cat_pops[which(cat_pops$time_since_last_neg_cat == "Cat3"), "V1"], ")")
cat_4_name <- paste0("Unknown", " (N=", cat_pops[which(cat_pops$time_since_last_neg_cat == "Unknown"), "V1"], ")")
total_name <- paste0("Total", " (N=", sum(cat_pops$V1), ")")
colnames(tab_cat) <- c("Characteristic", " ", cat_1_name, cat_2_name, cat_3_name, cat_4_name, total_name)

dat_tab1 <- flextable(tab_cat)
dat_tab1_format <- dat_tab1  %>%
  autofit() %>%
  merge_v(j = 1, part = "body", combine = FALSE) %>%
  bold(j = 1, bold = TRUE, part = "header")  %>%
  bold(i = 1, bold = TRUE, part = "header")  %>%
  flextable::align(align = "left", j = 1, part = "header") %>%
  flextable::align(align = "center", j = 3, part = "header") %>%
  flextable::align(align = "center", j = 2:3, part = "body") %>%
  hline(i = c(3, 5, 10), part = "body") %>%
  vline(part = "all", j = c(2, 6))
dat_tab1_format <- fix_border_issues(dat_tab1_format)

save_as_docx(dat_tab1_format, path = paste0(my_path_figures, "Table_1.docx"))

# check distribution of time from sympt onset to earliest test among symptomatic
table(dat[which(dat$Sympt_status == 1), "time_sympt_to_earliest_test"])

# among initially inconclusives, check number who used inconclusive vs subsequent positive test date as isolation start
time_neg_cat_breakdown <- ddply(dat, .(time_since_last_neg_cat), nrow)
time_neg_cat_breakdown$Percent <- time_neg_cat_breakdown$V1 / sum(time_neg_cat_breakdown$V1)
ddply(dat, .(incon_used_d0), nrow)

# calculate the mean and IQR of number of days symptoms appeared before the initial test, among the symptomatic
dat_sympt <- dat[which(dat$Sympt_status == 1), ]
ddply(dat_sympt, .(time_since_last_neg_cat), summarise, Mean_cat = mean(time_sympt_to_earliest_test, na.rm = TRUE),
      IQR_lo = quantile(time_sympt_to_earliest_test, prob = 0.25, na.rm = TRUE),
      IQR_hi = quantile(time_sympt_to_earliest_test, prob = 0.75, na.rm = TRUE))

###########################################################################
# Figure 1 - Daily Positivity by Time Since Last Negative Category  -----
###########################################################################
# Note: last negative test
## Cat1 = <= 4 days before dx
## Cat2 = 5-9 days before dx
## Cat3 = >= 10 days before dx

# drop persons w/ unknown time since last negative test (N=1)
dat <- dat %>% dplyr::filter(time_since_last_neg_cat != "Unknown")
unique(dat$time_since_last_neg_cat)

# drop initially inclusive persons who used subsequent positive test date as isolation start "pos" (N=7)
dat <- dat %>% dplyr::filter((is.na(incon_used_d0) == TRUE) | (incon_used_d0 == "Inc")) 
unique(dat$incon_used_d0)

# total persons in each time since last negative test category
dat_cat <- dat
total_by_cat <- ddply(dat_cat, .(time_since_last_neg_cat), nrow)
colnames(total_by_cat) <- c("Cat", "Total")
total_cat1 <- total_by_cat[which(total_by_cat$Cat == "Cat1"), "Total"] 
total_cat2 <- total_by_cat[which(total_by_cat$Cat == "Cat2"), "Total"] 
total_cat3 <- total_by_cat[which(total_by_cat$Cat == "Cat3"), "Total"] 

# generate number of isolation days based on data
isolation_days <- data.frame("Isolation_Day" = min(dat$earliest_test_to_release):max(dat$earliest_test_to_release))

# number ending isolation (i.e. testing negative) per day by category
dat_num_neg <- ddply(dat, .(earliest_test_to_release, time_since_last_neg_cat), nrow)
dat_num_neg <- reshape2::dcast(dat_num_neg, earliest_test_to_release ~ time_since_last_neg_cat, value.var = "V1")
colnames(dat_num_neg) <- c("Isolation_Day", "Num_Neg_Cat1", "Num_Neg_Cat2", "Num_Neg_Cat3")
dat_iso <- isolation_days %>% dplyr::left_join(dat_num_neg, by = "Isolation_Day")
dat_iso[is.na(dat_iso)] <- 0 # set days where no one tested negative to 0 from NA

# number positive per day
for (each_day in unique(dat_iso$Isolation_Day)) {
  
  # for day 5
  if (each_day == min(dat_iso$Isolation_Day)) { # num pos this day = total ppl at origin minus number negative on this day
    dat_iso[which(dat_iso$Isolation_Day == each_day), "Num_Pos_Cat1"] <- total_cat1 - dat_iso[which(dat_iso$Isolation_Day == each_day), "Num_Neg_Cat1"]
    dat_iso[which(dat_iso$Isolation_Day == each_day), "Num_Pos_Cat2"] <- total_cat2 - dat_iso[which(dat_iso$Isolation_Day == each_day), "Num_Neg_Cat2"]
    dat_iso[which(dat_iso$Isolation_Day == each_day), "Num_Pos_Cat3"] <- total_cat3 - dat_iso[which(dat_iso$Isolation_Day == each_day), "Num_Neg_Cat3"]
  } 
  # for every day after day 5
  else { # num pos this day = num pos on prior day minus num neg on this day
    dat_iso[which(dat_iso$Isolation_Day == each_day), "Num_Pos_Cat1"] <- dat_iso[which(dat_iso$Isolation_Day == (each_day - 1)), "Num_Pos_Cat1"] - dat_iso[which(dat_iso$Isolation_Day == each_day), "Num_Neg_Cat1"]
    dat_iso[which(dat_iso$Isolation_Day == each_day), "Num_Pos_Cat2"] <- dat_iso[which(dat_iso$Isolation_Day == (each_day - 1)), "Num_Pos_Cat2"] - dat_iso[which(dat_iso$Isolation_Day == each_day), "Num_Neg_Cat2"]
    dat_iso[which(dat_iso$Isolation_Day == each_day), "Num_Pos_Cat3"] <- dat_iso[which(dat_iso$Isolation_Day == (each_day - 1)), "Num_Pos_Cat3"] - dat_iso[which(dat_iso$Isolation_Day == each_day), "Num_Neg_Cat3"]
  }
}

# percent still positive (of those diagnosed at origin (num pos on each day divided by num at origin)
dat_iso$Perc_Still_Pos_Cat1 <- round(dat_iso$Num_Pos_Cat1 / total_cat1, digits = 2)
dat_iso$Perc_Still_Pos_Cat2 <- round(dat_iso$Num_Pos_Cat2 / total_cat2, digits = 2)
dat_iso$Perc_Still_Pos_Cat3 <- round(dat_iso$Num_Pos_Cat3 / total_cat3, digits = 2)

# summarize number positive
dat_iso_num_pos <- dat_iso %>% dplyr::select(Isolation_Day, Num_Pos_Cat1, Num_Pos_Cat2, Num_Pos_Cat3)
dat_iso_num_pos <- reshape2::melt(dat_iso_num_pos, .(Isolation_Day))
colnames(dat_iso_num_pos) <- c("Isolation_Day", "Category", "Number_Positive")
levels(dat_iso_num_pos$Category)
dat_iso_num_pos$Category <- factor(dat_iso_num_pos$Category, labels = c("Cat1", "Cat2", "Cat3"))

# summarize percent still positive
dat_iso_pct_still_pos <- dat_iso %>% dplyr::select(Isolation_Day, Perc_Still_Pos_Cat1, Perc_Still_Pos_Cat2, Perc_Still_Pos_Cat3)
dat_iso_pct_still_pos <- reshape2::melt(dat_iso_pct_still_pos, .(Isolation_Day))
colnames(dat_iso_pct_still_pos) <- c("Isolation_Day", "Category", "Pct_Still_Pos")
levels(dat_iso_pct_still_pos$Category)
dat_iso_pct_still_pos$Category <- factor(dat_iso_pct_still_pos$Category, labels = c("Cat1", "Cat2", "Cat3"))

# combine num positive and pct still positive 
dat_iso_all <- dat_iso_num_pos %>% dplyr::left_join(dat_iso_pct_still_pos, by = c("Isolation_Day", "Category"))

# plot 
sec_axis_scale <- max(dat_iso_all$Pct_Still_Pos, na.rm = TRUE) / max(dat_iso_all$Number_Positive, na.rm = TRUE) 
sec_axis_scale <- 0.005
max_num_pos <- max(dat_iso_all$Number_Positive)
max_num_pos <- 100
max_still_pos <- max(dat_iso_all$Pct_Still_Pos)
max_still_pos <- 0.5

## cat 1
dat_plot <- dat_iso_all %>% dplyr::filter(Category == "Cat1")
dat_plot$Label <- "Last Negative Test <= 4 Days"
p_cat1 <- ggplot(dat_plot) + 
  geom_bar(aes(x = Isolation_Day, y = Number_Positive*sec_axis_scale, fill = Category), alpha = 0.5, stat = "identity", position = "dodge") +
  geom_line(aes(x = Isolation_Day, y = Pct_Still_Pos), color = "black") +
  geom_point(aes(x = Isolation_Day, y = Pct_Still_Pos), color = "black") +
  geom_text(aes(x = Isolation_Day, y = Pct_Still_Pos,
                label = scales::percent(Pct_Still_Pos, accuracy = 1)),
            nudge_x = 0.5, nudge_y = 0.01, show.legend=FALSE, color = "black") +
  theme_bw() +
  facet_wrap(~Label) +
  scale_x_continuous(breaks = seq(1, 20, by = 2)) + 
  scale_y_continuous(name = "Percent Still Positive",
                     labels = function(b) { paste0(round(b * 100, 0), "%")},
                     limits = c(0, max_still_pos),
                     sec.axis = sec_axis(~./sec_axis_scale, name = NULL, breaks = c(0, 20, 40, 60, 80, 100))) + 
  theme(
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black")) + 
  theme(legend.position = "bottom", legend.box = "horizontal") +
  xlab("Isolation Day") +
  scale_fill_discrete(labels = c("<= 4", "5-9", ">=10")) +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.position = "none",
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))

## cat 2
dat_plot <- dat_iso_all %>% dplyr::filter(Category == "Cat2")
dat_plot$Label <- "Last Negative Test 5-9 Days"
p_cat2 <- ggplot(dat_plot) + 
  geom_bar(aes(x = Isolation_Day, y = Number_Positive*sec_axis_scale, fill = Category), alpha = 0.5, stat = "identity", position = "dodge") +
  geom_line(aes(x = Isolation_Day, y = Pct_Still_Pos), color = "black") +
  geom_point(aes(x = Isolation_Day, y = Pct_Still_Pos), color = "black") +
  geom_text(aes(x = Isolation_Day, y = Pct_Still_Pos,
                label = scales::percent(Pct_Still_Pos, accuracy = 1)),
            nudge_x = 0.5, nudge_y = 0.01, show.legend=FALSE, color = "black") +
  theme_bw() +
  facet_wrap(~Label) +
  scale_x_continuous( breaks = seq(1, 20, by = 2)) + 
  scale_y_continuous(name = NULL,
                     labels = function(b) { paste0(round(b * 100, 0), "%")},
                     limits = c(0, max_still_pos),
                     sec.axis = sec_axis(~./sec_axis_scale, name = NULL, breaks = c(0, 20, 40, 60, 80, 100))) + 
  theme(
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black")) + 
  theme(legend.position = "bottom", legend.box = "horizontal") +
  xlab("Isolation Day") +
  scale_fill_discrete(labels = c("<= 4", "5-9", ">=10")) +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.position = "none",
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))

## cat 3
dat_plot <- dat_iso_all %>% dplyr::filter(Category == "Cat3")
dat_plot$Label <- "Last Negative Test >= 10 Days"
p_cat3 <- ggplot(dat_plot) + 
  geom_bar(aes(x = Isolation_Day, y = Number_Positive*sec_axis_scale, fill = Category), alpha = 0.5, stat = "identity", position = "dodge") +
  geom_line(aes(x = Isolation_Day, y = Pct_Still_Pos), color = "black") +
  geom_point(aes(x = Isolation_Day, y = Pct_Still_Pos), color = "black") +
  geom_text(aes(x = Isolation_Day, y = Pct_Still_Pos,
                label = scales::percent(Pct_Still_Pos, accuracy = 1)),
            nudge_x = 0.5, nudge_y = 0.01, show.legend=FALSE, color = "black") +
  theme_bw() +
  facet_wrap(~Label) +
  scale_x_continuous( breaks = seq(1, 20, by = 2)) + 
  scale_y_continuous(name = NULL,
                     limits = c(0, max_still_pos),
                     labels = function(b) { paste0(round(b * 100, 0), "%")},
                     sec.axis = sec_axis(~./sec_axis_scale, name = "Number Positive",  breaks = c(0, 20, 40, 60, 80, 100))) + 
  theme(
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black")) + 
  theme(legend.position = "bottom", legend.box = "horizontal") +
  xlab("Isolation Day") +
  scale_fill_discrete(labels = c("<= 4", "5-9", ">=10")) +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.position = "none",
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))

# isolation duration
max_iso <- max(dat$earliest_test_to_release) 
p_iso_dur <- ggplot(dat, aes(x = earliest_test_to_release, fill = time_since_last_neg_cat)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 5, color = "black", linetype="dashed",) +
  ylab("Density") +
  xlab("Isolation Days: Earliest Test to Release") + 
  labs(fill = "Days Since Last Negative Test") +
  scale_fill_manual(labels = c("<= 4", "5-9", ">=10"), values = customPalette) +
  scale_x_continuous(limits = c(0, max_iso+1), breaks = seq(1, max_iso+1, by = 2)) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))

mylegend <- get_legend(p_iso_dur + theme(legend.position='bottom'))

p <- grid.arrange(ggarrange(p_cat1 + 
                              theme(legend.position="none") + 
                              scale_fill_manual(values = customPalette[1]) + scale_color_manual(values = customPalette[1]),
                            p_cat2 + 
                              theme(legend.position="none") + 
                              scale_fill_manual(values = customPalette[2]) + scale_color_manual(values = customPalette[2]),
                            p_cat3 + 
                              theme(legend.position="none") + 
                              scale_fill_manual(values = customPalette[3]) + scale_color_manual(values = customPalette[3]),
                            nrow = 1,
                            labels = c("A", "B", "C")),
                  ggarrange(p_iso_dur + theme(legend.position="none"),
                            nrow = 1,
                            labels = c("D")),
                  mylegend,
                  nrow = 3,
                  heights = c(5.5, 4, 0.5))

ggsave(paste(my_path_figures, "Fig_1.pdf", sep=""), p, height = 10, width =10)

# count number of persons in each time since the last negative test category
ddply(dat, .(time_since_last_neg_cat), nrow)
      
###########################################################################
# Table 2: Parametric Regression Survival Analysis ------------------------
###########################################################################

dat_surv <- data.frame(dat)

# create interval censoring variable
## time 1 - starting time for interval
dat_surv[, "time_1"] <- dat_surv[, "earliest_test_to_release"]
dat_surv[which(dat_surv$earliest_test_to_release == 5), "time_1"] <- 1 # for anyone testing negative on day 5, interval start is day 1

## time 2 - ending time for interval
dat_surv[, "time_2"] <- dat_surv[, "earliest_test_to_release"]

# drop persons missing CT values (N=27)
ddply(dat_surv, .(CT_at_dx), nrow)
dat_surv <- dat_surv %>% dplyr::filter(is.na(CT_at_dx) == FALSE) 

# drop persons missing symptom status (N=2)
ddply(dat_surv, .(Sympt_status), nrow)
dat_surv <- dat_surv %>% dplyr::filter(Sympt_status != "Unknown") 

# drop persons with 1, 4, of Unknown vaccine doses (N=15)
ddply(dat_surv, .(Num_vx_doses), nrow) 
dat_surv <- dat_surv %>% dplyr::filter(!(Num_vx_doses %in% c("1", "4", "Unknown")))
ddply(dat_surv, .(Num_vx_doses), nrow) # now only includes 2 or 3 vaccine dose persons

# create primary vaccine brand categories
ddply(dat_surv, .(Primary_vx_brand, Other_vx_first), nrow) 
dat_surv[which(dat_surv$Primary_vx_brand %in% c("Pfizer", "Moderna")), "Primary_vx_brand_cat"] <- "mRNA"
dat_surv[which((dat_surv$Primary_vx_brand %in% c("Pfizer", "Moderna") & (is.na(dat_surv$Other_vx_first) == FALSE))), "Primary_vx_brand_cat"] <- "Other_mRNA" # relabel mRNA primary vx series persons who received an intl vx before
dat_surv[which(dat_surv$Primary_vx_brand == "JJ"), "Primary_vx_brand_cat"] <- "JJ"
dat_surv[which(dat_surv$Primary_vx_brand == "Other"), "Primary_vx_brand_cat"] <- "Other" # received intl vx as primary series

# drop persons who received an international vx prior to their primary series or as their primary series (N=8)
ddply(dat_surv, .(Primary_vx_brand_cat), nrow) 
dat_surv <- dat_surv %>% dplyr::filter(Primary_vx_brand_cat %in% c("JJ", "mRNA"))
ddply(dat_surv, .(Primary_vx_brand_cat), nrow) # now only includes mRNA and J&J categories

# create numbers and timing of last vaccine dose categories
num_mth <- 5 # months
time_cutoff <- 31 * num_mth # convert to days
dat_surv[which((dat_surv$Num_vx_doses == 2) & (dat_surv$time_since_last_dose < time_cutoff)), "num_vx_time_dose"] <- paste0("Two_less_", num_mth, "m")
dat_surv[which((dat_surv$Num_vx_doses == 2) & (dat_surv$time_since_last_dose >= time_cutoff)), "num_vx_time_dose"] <- paste0("Two_over_", num_mth, "m")
dat_surv[which((dat_surv$Num_vx_doses == 3) & (dat_surv$time_since_last_dose < time_cutoff)), "num_vx_time_dose"] <- paste0("Three_less_", num_mth, "m")
dat_surv[which((dat_surv$Num_vx_doses == 3) & (dat_surv$time_since_last_dose >= time_cutoff)), "num_vx_time_dose"] <- paste0("Three_over_", num_mth, "m")
ddply(dat_surv, .(num_vx_time_dose), nrow) 
dat_surv[which(dat_surv$Num_vx_doses == 3), "num_vx_time_dose"] <- "Three"

# select columns
dat_surv_sub <- dat_surv %>% dplyr::select("earliest_test_to_release",
                                           "time_1",
                                           "time_2",
                                           "time_since_last_neg_cat",
                                           "Sympt_status",
                                           "CT_at_dx",
                                           "Prior_inf_more90d",
                                           "num_vx_time_dose",
                                           "Primary_vx_brand_cat")
# "earliest_test_date")

# save(dat_surv_sub, file = paste0(my_path, "dat_surv_sub.rda"))

# set reference levels and labels
unique(dat_surv_sub$time_since_last_neg_cat)
dat_surv_sub$time_since_last_neg_cat <- factor(dat_surv_sub$time_since_last_neg_cat,
                                               levels = c("Cat1", "Cat2", "Cat3"),
                                               labels = c("<=4 days", "5-9 days", ">= 10 days"))

unique(dat_surv_sub$Sympt_status)
dat_surv_sub$Sympt_status <- factor(dat_surv_sub$Sympt_status,
                                    levels = c("0", "1"),
                                    labels = c("No", "Yes"))

unique(dat_surv_sub$Prior_inf_more90d)
dat_surv_sub$Prior_inf_more90d <- factor(dat_surv_sub$Prior_inf_more90d,
                                         levels = c("0", "1"),
                                         labels = c("No", "Yes"))

unique(dat_surv_sub$num_vx_time_dose)
dat_surv_sub$num_vx_time_dose <- factor(dat_surv_sub$num_vx_time_dose,
                                        levels = c(paste0("Two_over_", num_mth, "m"), paste0("Two_less_", num_mth, "m"), "Three"),
                                        labels = c(paste0("2 Doses / >=", num_mth, "m"), paste0("2 Doses / <", num_mth, "m"), "3 Doses"))

unique(dat_surv_sub$Primary_vx_brand_cat)
dat_surv_sub$Primary_vx_brand_cat <- factor(dat_surv_sub$Primary_vx_brand_cat,
                                            levels = c("JJ", "mRNA"),
                                            labels = c("JJ", "mRNA"))

# set up survival object
## if type = interval2, value of time2 argument ignored unless event=3
## all day 5 test negative persons have interval [1, 5]; all others have test negative day
surv_obj <- Surv(time = dat_surv_sub$time_1, time2 = dat_surv_sub$time_2, type = "interval2")

# exploratory K-M curves
## overall fit
initial_fit <- survfit(surv_obj ~ 1) # fits single survival curve when ~ 1
plot(initial_fit, xlab = "Isolation Duration", ylab = "Probability of survival (testing pos)")

## stratified by categorical covariates
plot_cov_KM(dat_surv_sub$time_since_last_neg_cat)
plot_cov_KM(dat_surv_sub$Sympt_status)
plot_cov_KM(dat_surv_sub$Prior_inf_more90d)
plot_cov_KM(dat_surv_sub$num_vx_time_dose)
plot_cov_KM(dat_surv_sub$Primary_vx_brand_cat)

# compare model fits using different distributions
dist_list <- c("weibull", "exponential", "lognormal", "loglogistic", "logistic")
dist_results <- NULL
for (each_dist in dist_list) {
  dat_res <- survreg(surv_obj ~
                       time_since_last_neg_cat +
                       Sympt_status +
                       CT_at_dx +
                       Prior_inf_more90d +
                       num_vx_time_dose +
                       Primary_vx_brand_cat,
                     dist = each_dist,
                     data = dat_surv_sub)
  dist_AIC <- extractAIC(dat_res)[2] # computes AIC for a fitted parametric model
  dist_results_single <- cbind.data.frame(each_dist, dist_AIC)
  dist_results <- rbind.data.frame(dist_results, dist_results_single)
}

colnames(dist_results) <- c("Distribution", "AIC")
dist_results <- dist_results[order(dist_results$AIC, decreasing = FALSE), ]

# check how much results change with distribution selection
dist_list <- c("weibull", "exponential", "lognormal", "loglogistic", "logistic", "gaussian")
dist_results <- NULL
dist_results_p <- NULL
for (each_dist in dist_list) {
  dat_res <- survreg(surv_obj ~
                       time_since_last_neg_cat +
                       Sympt_status +
                       CT_at_dx +
                       Prior_inf_more90d +
                       num_vx_time_dose +
                       Primary_vx_brand_cat,
                     dist = each_dist,
                     data = dat_surv_sub)
  dist_coef <- exp(dat_res$coefficients)  # computes AIC for a fitted parametric model
  dist_results <- rbind.data.frame(dist_results, dist_coef)
  dat_res_sum <- summary(dat_res)
  dist_p_single <- dat_res_sum$table[, "p"]
  dist_results_p <- rbind.data.frame(dist_results_p, dist_p_single)
}

## regression coefficients
rownames(dist_results) <- dist_list
colnames(dist_results) <- names(dat_res$coefficients) # doesn't matter which grab names from; all the smae

## p-values
rownames(dist_results_p) <- dist_list
colnames(dist_results_p) <- names(dat_res$coefficients) # NA at end is for log(scale)

# fit accelerated failure time model with selected distribution
## AFT measuring effect of covs to accelerate/decelerate survival time
## regression coefficients = logs of ratios of survival times
## predictors (once exponentiated) act multiplicatively on the failure time 
dat_res <- survreg(surv_obj ~
                     time_since_last_neg_cat +
                     Sympt_status +
                     CT_at_dx +
                     Prior_inf_more90d +
                     num_vx_time_dose +
                     Primary_vx_brand_cat,
                   dist = "lognormal",
                   data = dat_surv_sub)

# check whether the ratio of time-quantile (the acceleration factor) is constant for all fixed values via QQ plot
## QQ plot shows predicted survival time quantitles of cov levels against each other - pairs will be on straight line thru origin if assumption valid

library("assertive")
library("snowfall")
library("interval")
library("Icens")
AFTplot(dat_res)

# convert regression coefficient to ETR by exponentiating
dat_res_conv <- exp(cbind(dat_res$coefficients, confint(dat_res))) 
dat_res_conv <- dat_res_conv[2:nrow(dat_res_conv), ] # drop intercept

# extract p-values
dat_res_summary <- summary(dat_res)
p_vals <- data.frame(dat_res_summary$table[, "p"])
p_vals <- p_vals[2:(nrow(p_vals)-1), ] # drop intercept, log(scale)

# combine ETR and p-values into 1 df
dat_res_ETR <- cbind.data.frame(dat_res_conv, p_vals)
colnames(dat_res_ETR) <- c("ETR", "LB", "UB", "p_vals")
cov_names <- c(rep("Time Since Last Negative Test", 3),
               rep("Symptoms at Dx", 2),
               rep("CT Value at Dx", 1),
               rep("Prior Infection >90 Days", 2),
               rep("No. Dose/Time Since Last", 3),
               rep("Primary Vaccine Brand", 2))

# rename columns
dat_surv_plot <- dat_surv_sub %>% dplyr::rename("Isolation Duration" = "earliest_test_to_release",
                                                "Time Since Last Negative Test" = "time_since_last_neg_cat",
                                                "Symptoms at Dx" = "Sympt_status",
                                                "CT Value at Dx" = "CT_at_dx",
                                                "Prior Infection >90 Days" = "Prior_inf_more90d",
                                                "No. Dose/Time Since Last" = "num_vx_time_dose",
                                                "Primary Vaccine Brand" = "Primary_vx_brand_cat")

# categories for each covariate
cat_names <- c(levels(dat_surv_plot$`Time Since Last Negative Test`),
               levels(dat_surv_plot$`Symptoms at Dx`),
               "-", # for CT at dx; no levels since continuous variable
               levels(dat_surv_plot$`Prior Infection >90 Days`),
               levels(dat_surv_plot$`No. Dose/Time Since Last`),
               levels(dat_surv_plot$`Primary Vaccine Brand`))

# indicate reference group
ref_col <- c("ref", "-", "-",
             "ref", "-",
             "-",
             "ref", "-",
             "ref", "-", "-",
             "ref", "-")

# grab sample sizes for each covariate category
sample_size <- c(ddply(dat_surv_plot, .(`Time Since Last Negative Test`), nrow)$V1,
                 ddply(dat_surv_plot, .(`Symptoms at Dx`), nrow)$V1,
                 nrow(dat_surv_plot),
                 ddply(dat_surv_plot, .(`Prior Infection >90 Days`), nrow)$V1,
                 ddply(dat_surv_plot, .(`No. Dose/Time Since Last`), nrow)$V1,
                 ddply(dat_surv_plot, .(`Primary Vaccine Brand`), nrow)$V1)
                
# combine all columns into 1 df
dat_all_ETR <- cbind.data.frame(cov_names, cat_names, ref_col, sample_size)

# create results table
ETR_output <- build_results_table(dat_res_ETR, dat_all_ETR, c("ETR", "LB", "UB", "p_vals"))
ETR_output$CI <- paste0("(", ETR_output$LB, "-", ETR_output$UB, ")")
ETR_output[ETR_output == "(---)"] <- "-"
ETR_output <- ETR_output[, c("cov_names", "cat_names", "ref_col", "sample_size", "ETR", "CI", "p_vals")]
tab_ETR <- flextable(ETR_output)
tab_ETR_format <- tab_ETR %>% 
  set_header_labels(       
    cov_names = "Covariate", 
    cat_names = "Category",                  
    ref_col = " ",
    sample_size = "Sample Size",
    ETR = "ETR",
    CI = "95% CI",
    p_vals = "P-Value") %>%
  align(align = "left", part = "header") %>%
  align(align = "left", part = "body") %>%
  bold(i = 1, bold = TRUE, part = "header") %>%
  vline(part = "all", j = 4) %>% 
  merge_v(j = 1, part = "body", combine = FALSE) %>%
  bg(i = c(1:3, 6, 9:11), part = "body", bg = "gray95")  
tab_ETR_format <- fix_border_issues(tab_ETR_format)
tab_ETR_format <- fontsize(tab_ETR_format, size = 8, part = "all")

save_as_docx(tab_ETR_format, path = paste0(my_path_figures, "Table_2.docx"))

###########################################################################
# Figure S1: Isolation Duration by Covariate/Category  ---------------------
###########################################################################

p_sympt <- ggplot(dat_surv_plot, aes(x = `Isolation Duration`, fill = `Symptoms at Dx`)) +
  geom_bar(alpha = 0.6, position = "stack") +
  geom_vline(xintercept = 5, color = "black", linetype="dashed",) +
  ylab("No. Persons") +
  xlab("Isolation Duration (Days)") +
  labs(fill = "Symptoms at Diagnosis") +
  scale_fill_manual(values = customPalette2) +
  scale_x_continuous(limits = c(0, max_iso+1), breaks = seq(1, max_iso+1, by = 2)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10))

p_prior <- ggplot(dat_surv_plot, aes(x = `Isolation Duration`, fill = `Prior Infection >90 Days`)) +
  geom_bar(alpha = 0.6, position = "stack") +
  geom_vline(xintercept = 5, color = "black", linetype="dashed",) +
  ylab("No. Persons") +
  xlab("Isolation Duration (Days)") +
  labs(fill = "Prior Infection >90 Days") +
  scale_fill_manual(values = customPalette2) +
  scale_x_continuous(limits = c(0, max_iso+1), breaks = seq(1, max_iso+1, by = 2)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10))

p_num_time <- ggplot(dat_surv_plot, aes(x = `Isolation Duration`, fill = `No. Dose/Time Since Last`)) +
  geom_bar(alpha = 0.6, position = "stack") +
  geom_vline(xintercept = 5, color = "black", linetype="dashed",) +
  ylab("No. Persons") +
  xlab("Isolation Duration (Days)") +
  labs(fill = "No. Doses/Timing") +
  scale_fill_manual(values = customPalette2) +
  scale_x_continuous(limits = c(0, max_iso+1), breaks = seq(1, max_iso+1, by = 2)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10))

p_brand <- ggplot(dat_surv_plot, aes(x = `Isolation Duration`, fill = `Primary Vaccine Brand`)) +
  geom_bar(alpha = 0.6, position = "stack") +
  geom_vline(xintercept = 5, color = "black", linetype="dashed",) +
  ylab("No. Persons") +
  xlab("Isolation Duration (Days)") +
  labs(fill = "Primary Vaccine Brand") +
  scale_fill_manual(values = customPalette2) +
  scale_x_continuous(limits = c(0, max_iso+1), breaks = seq(1, max_iso+1, by = 2)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10))

p_ct <- ggplot(dat_surv_plot, aes(x = `Isolation Duration`, y = `CT Value at Dx`)) +
  geom_point(color = customPalette2[1], alpha = 0.5, position = position_jitter(width = 0.3, height = 0.2), size = 1) +
  stat_summary(color = customPalette2[2], geom = "point", fun = "mean", color = "black", size = 1) +
  geom_errorbar(color = customPalette2[2], stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width = 0.1, color = "blue") +
  scale_y_continuous(trans = "reverse") +
  scale_x_continuous(limits = c(0, max_iso+1), breaks = seq(1, max_iso+1, by = 2)) +
  xlab("Isolation Duration (Days)") +
  ylab("CT Value at Diagnosis") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

p <- grid.arrange(ggarrange(p_sympt + theme(legend.position="bottom", legend.text=element_text(size=8)),
                            p_prior + theme(legend.position="bottom", legend.text=element_text(size=8)),
                            p_num_time + theme(legend.position="bottom", legend.text=element_text(size=8)),
                            p_brand + theme(legend.position="bottom", legend.text=element_text(size=8)),
                            nrow = 2,
                            ncol = 2,
                            labels = c("A", "B", "C", "D")),
                  ggarrange(p_ct + theme(legend.position="none"),
                            nrow = 1,
                            labels = "E"),
                  nrow = 2,
                  heights = c(7, 3))

ggsave(paste(my_path_figures, "Fig_S1.pdf", sep=""), p, height = 10, width =10)


###########################################################################
# Figure S2: Number/Timing of Vaccine Doses and Diagnosis Date  -----------
###########################################################################
# Note: "earliest_test_date" not included in shared dataset; please contact corresponding author

# # num/timing of doses by diagnosis time 
# p_num_time_timeline <- ggplot(dat_surv_plot, aes(x = earliest_test_date, fill = `No. Dose/Time Since Last`)) +
#   geom_bar(alpha = 0.6, position = "stack") +
#   ylab("No. Persons") +
#   xlab("Diagnosis Date") +
#   labs(fill = "No. Doses/Timing") +
#   scale_fill_manual(values = customPalette2) +
#   theme_bw() +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 12),
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=10))
# 
# # smoothed proportion of diagnosed persons over time belonging to each num/timing of doses category
# dat_prop <- ddply(dat_surv_plot, .(earliest_test_date, `No. Dose/Time Since Last`), nrow)
# dat_prop_cast <- reshape2::dcast(dat_prop, earliest_test_date ~ `No. Dose/Time Since Last`, value.var = "V1")
# dat_prop_cast[is.na(dat_prop_cast)] <- 0
# dat_prop_cast$Total <- dat_prop_cast$`2 Doses / <5m` + dat_prop_cast$`2 Doses / >=5m` + dat_prop_cast$`3 Doses`
# dat_prop_cast$Prop_Three <- dat_prop_cast$`3 Doses`/ dat_prop_cast$Total
# dat_prop_cast$Prop_Two_less5 <- dat_prop_cast$`2 Doses / <5m` / dat_prop_cast$Total
# dat_prop_cast$Prop_Two_more5 <- dat_prop_cast$`2 Doses / >=5m`/ dat_prop_cast$Total
# dat_prop_melt <- reshape2::melt(dat_prop_cast, id.vars = c("earliest_test_date"),
#                                 measure.vars = c("Prop_Two_more5", "Prop_Two_less5", "Prop_Three"))
# 
# p_num_time_prop <- ggplot(dat_prop_melt, aes(x = earliest_test_date, y = value, color = variable, fill = variable)) +
#   geom_smooth() +
#   ylab("Proportion of Diagnosed Persons") +
#   xlab("Diagnosis Date") +
#   labs(color = "No. Doses/Timing", fill = "No. Doses/Timing") +
#   scale_color_manual(values = customPalette2, labels = c("2 Doses / >=5m", "2 Doses / <5m", "3 Doses")) +
#   scale_fill_manual(values = customPalette2, labels = c("2 Doses / >=5m", "2 Doses / <5m", "3 Doses")) +
#   theme_bw() +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 12),
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=10))
# 
# mylegend <- get_legend(p_num_time_prop + theme(legend.position='bottom'))
# 
# p <- grid.arrange(ggarrange(p_num_time_timeline + theme(legend.position="none"),
#                             p_num_time_prop + theme(legend.position="none"),
#                             nrow = 1,
#                             labels = c("A", "B")),
#                   mylegend,
#                   nrow = 2,
#                   heights = c(4.5, 0.5))
# 
# ggsave(paste(my_path_figures, "Fig_S2.pdf", sep=""), p, height = 5, width =10)
# 
# # extract smoothed proportions data
# dat_prop_smooth <- ggplot_build(p_num_time_prop)$data[[1]]
# group_1 <- dat_prop_smooth[which(dat_prop_smooth$group == 1), "y"]
# group_1_start <- round(group_1[1], digits = 2)
# group_1_stop <- round(group_1[length(group_1)], digits = 2)
# group_2 <- dat_prop_smooth[which(dat_prop_smooth$group == 2), "y"]
# group_2_start <- round(group_2[1], digits = 2)
# group_2_stop <- round(group_2[length(group_2)], digits = 2)
# group_3 <- dat_prop_smooth[which(dat_prop_smooth$group == 3), "y"]
# group_3_start <- round(group_3[1], digits = 2)
# group_3_stop <- round(group_3[length(group_3)], digits = 2)
