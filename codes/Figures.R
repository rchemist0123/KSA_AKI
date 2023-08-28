

# figure 2 ----------------------------------------------------------------

aki_stage_melt = df[,.SD,.SDcols=patterns("AKI_Day|AKI_initial")] |> 
  melt(variable.name = "Day", value.name="Stage")
aki_stage_melt_prop = aki_stage_melt[!is.na(Stage), .N,by=.(Day, Stage)][order(Day, Stage)][,prop:= round(N/sum(N)*100,1), by=.(Day)]
aki_stage_melt_prop_sum = aki_stage_melt_prop[!is.na(Stage) & Stage!=0,
                                              .(Stage = "Any AKI", prop = sum(prop)),by=Day]
aki_stage_melt_final = rbind(aki_stage_melt_prop,aki_stage_melt_prop_sum, fill=T)
aki_stage_melt_final[,Stage2 := ifelse(substr(Stage,1,1)!="A",paste0("Stage ",Stage," AKI"), Stage) |>
                       factor(levels=c("Stage 1 AKI","Stage 2 AKI", "Stage 3 AKI","Stage 7 AKI",  "Any AKI"))]
aki_stage_melt_final[,Day := fcase(Day=="AKI_initial", 0,
                                   Day =="AKI_Day1", 1,
                                   Day =="AKI_Day2", 2,
                                   Day =="AKI_Day3", 3,
                                   Day =="AKI_Day7", 7,
                                   default = NA
) |> as.factor()]

aki_stage_melt_final|> 
  dplyr::filter(!is.na(Stage2)) |> 
  ggplot(aes(x=Day, y=prop, fill=as.factor(Stage2)))+
  geom_col(position = position_dodge(), color="black") +
  # scale_fill_brewer(palette = "Greys") + 
  scale_fill_grey(start=.9, end=.1) + 
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,60))+
  theme_classic() + 
  labs(x="Days in intensive Care Unit", y= "All Patients (%)") +
  theme(legend.title = element_blank(),
        legend.position = "top")
df[,.(ICUIntvRRT)]


# figure 3 ----------------------------------------------------------------

df[,urine_output_ICUD1 := (Output_ICUD1 / Wt / 24)]
df[,urine_output_ICUD2 := (Output_ICUD2 / Wt / 24)]
df[,urine_output_ICUD3 := (Output_ICUD3 / Wt / 24)]
df[,urine_output_ICUD7 := (Output_ICUD7 / Wt / 24)]

df[,AKI_Day1_urine := ifelse(urine_output_ICUD1 < 0.3, 1,0)]
df[,AKI_Day2_urine := ifelse(urine_output_ICUD2 < 0.3, 1,0)]
df[,AKI_Day3_urine := ifelse(urine_output_ICUD3 < 0.3, 1,0)]
df[,AKI_Day7_urine := ifelse(urine_output_ICUD7 < 0.3, 1,0)]
df[,AKI_stage_urine := do.call(pmax, c(.SD, na.rm=T)),.SDcols=patterns("_urine$")]
urine_cr_stage3_melt = df[!is.na(AKI_stage_urine),.(stage3_yn, AKI_stage_urine)] |> 
  melt() 
urine_cr_stage3_melt_prop = urine_cr_stage3_melt[,.N,.(variable,value)][,prop := round(N/sum(N)*100,1), by=.(variable)]
urine_cr_stage3_melt_prop[,name := fifelse(variable=="stage3_yn","AKI by SCr criteria",
                                           "AKI by UOP criteria")]
urine_cr_stage3_melt_prop[value==1] |> 
  ggplot(aes(x=as.factor(value), y=prop, fill=name)) +
  geom_col(position = position_dodge(), width=.5) +
  scale_x_discrete(label="Stage 3 AKI") +
  scale_y_continuous(expand=c(0,0),
                     limits = c(0,45),
                     breaks = seq(0,45,5))+
  scale_fill_grey() +
  labs(y="Percentile") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15)) +
  geom_label(aes(label="p=0.548",y=43), fill="white", size=8) 
geom_text(aes(label=prop, position=prop),
          position = position_dodge(width = .5),
          vjust=-.5)

mytable(stage3_yn ~ inhos_mortality, df)
mytable(AKI_stage_urine ~ inhos_mortality, df)
with(df, chisq.test(stage3_yn, AKI_stage_urine))
df[,.N,stage3_yn]
df[,.N,AKI_stage_urine]
df[,.N,.(stage3_yn, AKI_stage_urine)][!is.na(AKI_stage_urine)]
matrix(c(1924,2942, 1887,2961),nrow=2) |> chisq.test()



data.frame(
  x = rep("Stage 3 AKI",2),
  variable = c("AKI by SCr criteria","AKI by UOP criteria"),
  value = c(52.9, 48.3)
) |> 
  ggplot(aes(x=x, y=value, fill=variable))+
  geom_col(position = position_dodge(), width=.5) +
  scale_y_continuous(limits = c(0,60), expand = c(0,0)) + 
  theme_classic() + 
  scale_fill_grey() +
  labs(y="Mortality percentile") +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15)) +
  geom_label(aes(label="p=0.005", y=57), fill="white", size=8)

matrix(c(907,1017,976,911),nrow=2) |> chisq.test()


# figure 4 ----------------------------------------------------------------

# Average Cr by ICU days by CKD vs Non-CKD ------------------
## line + ribbon plot
df_aki_melt <- melt(df,
                    id.vars = c("SubjectNo", "AKI_stage","Chronic_kidney_ds"),
                    measure.vars = patterns("Cr_"),
                    variable.name = "Day")[Day !="Cr_baseline_MDRD" &
                                             Day != "Cr_ICUDL"]
df_aki_melt[,Day := fcase(
  Day %like% "baseline", "0",
  Day %like% "Cr_ICUD", sub("Cr_ICUD","",Day)
) |> as.factor()]
df_aki_melt[,.(mean = mean(value, na.rm=T),
               sd = sd(value, na.rm=T)),by=.(Day,AKI_stage,Chronic_kidney_ds)]|> 
  ggplot(aes(x=Day, y=mean,
             group = as.factor(AKI_stage)))+
  geom_line(aes(color=as.factor(AKI_stage))) +
  geom_ribbon(aes(y=mean, ymin=mean-sd ,ymax = mean+sd, fill=as.factor(AKI_stage)),
              alpha=.2) +
  geom_point(aes(color=as.factor(AKI_stage))) +
  scale_y_continuous(expand=c(0,0),
                     limits=c(0,8),
                     breaks=seq(0,8,1))+
  scale_color_discrete(name="AKI Stage") + 
  scale_fill_discrete(name="AKI Stage") + 
  theme_classic() + 
  facet_grid(~Chronic_kidney_ds,
             labeller = labeller(Chronic_kidney_ds=ckd_label))+
  labs(x="ICU Days",
       y= "Average Cr (mg/dL)") +
  theme(legend.position = "top")


# figure 5 ----------------------------------------------------------------

install.packages("forestploter")
library(forestploter)
forest_df <- data.table(
  est = c("Non-AKI","Stage 1","Stage 2", "Stage 3"),
  hr = c(1.00, 1.00, 1.06, 1.32),
  ci_low = c(1.00, 0.80, 0.86, 1.06),
  ci_high = c(1.00, 1.26, 1.32, 1.64)
)

forest_df$`HR (95% CI)` = sprintf("%.2f (%.2f - %.2f)", 
                                  forest_df$hr, forest_df$ci_low, forest_df$ci_high)

forest_df$`    ` = paste0(rep(" ",20), collapse = " ")
p = forest(forest_df[,c(1,6,5)],
           est = forest_df$hr,
           lower = forest_df$ci_low,
           upper = forest_df$ci_high,
           ref_line = 1,
           ci_column = 4,
           xlim = c(0.8, 1.6),
           ticks_at = c(0.8, 1.0, 1.2, 1.4, 1.6),
           theme = forest_theme(base_size=25,
                                ci_lwd = 5))
plot(p)
# AKI group In-hospital mortality (Kaplan-Meier, KM) ---------
require(survival)
require(survminer)
names(df)
load("AKI.RData")

aki = "AKI_stage"
km_inhos_by_aki_stage = function(aki){
  form = paste0("Surv(inhos_duration, inhos_mortality==1) ~ ", aki)
  inhos.aki.group.km_fit = survfit(as.formula(form), data=df)
  inhos_aki_group_km = ggsurvplot(inhos.aki.group.km_fit,
                                  xlim= c(0,100),
                                  break.time.by=10,
                                  censor=F,
                                  pval=T,
                                  risk.table=T,
                                  linetype = c("dotted","dotdash","dashed", "solid"),
                                  legend.labs = c("Non-AKI","Stage 1", "Stage 2", "Stage 3"),
                                  palette="jama",
                                  pval.coord = c(1, 0.2),
                                  xlab="Time (day)"
                                  # surv.median.line="hv",
  )
  inhos_aki_group_km$plot  = inhos_aki_group_km$plot + 
    # scale_color_brewer(palette="Greys") +
    theme(legend.direction = "vertical",
          legend.position = c(0.9, 0.9),
          legend.title = element_blank()) +
    scale_y_continuous(expand=c(0,0))
  
  inhos_aki_group_km$table = inhos_aki_group_km$table +
    theme(axis.title.x= element_blank(),
          axis.title.y = element_blank())
  
  print(inhos_aki_group_km)
}
km_inhos_by_aki_stage(aki="AKI_stage")

require(survminer)
require(survival)

ggsave(paste0("plots/inhos_aki_kmplot.png"),
       width= 150,
       height=120,
       units="mm",
       dpi=300
)

temp = pairwise_survdiff(Surv(inhos_duration, inhos_mortality)~AKI_stage, df, p.adjust.method="bonferroni")
temp$p.value |> as.data.table(keep.rownames = T) |> gt()
pairwise_survdiff(Surv(Hospital_LOS, inhos_mortality)~AKI_stage_lowest, df, p.adjust.method="bonferroni")
