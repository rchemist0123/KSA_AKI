

# figure 2 ----------------------------------------------------------------
aki = "AKI_stage"
km_inhos_by_aki_stage = function(aki){
  form = paste0("Surv(inhos_duration, inhos_mortality==1) ~ ", aki)
  form = paste0("Surv(icu_duration, icu_mortality==1) ~ ", aki)
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
form = paste0("Surv(TZ_to_disch_duration, inhos_mortality==1) ~ ", aki)
form = paste0("Surv(inhos_duration, inhos_mortality==1) ~ ", aki)
pairwise_survdiff(as.formula(form),data=df)

require(survminer)
require(survival)

ggsave(paste0("plots/inhos_aki_kmplot.png"),
       width= 150,
       height=120,
       units="mm",
       dpi=300
)
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


# Figure S5. I/O Restricted cubic spline plot ----------------------------------------------------------------
library(rms)
library(plotRCS)

getRCSplot = function(x, covars, outcome, time=NULL){
  vars = c(x, covars, outcome, time)
  rcsplot(data = df[,..vars] |> na.omit(),
          outcome = outcome,
          time = time,
          exposure = x,
          covariates = covars,
          na.rm=T) +
    labs(x='Net I/O')
}
vars = c('input_output_ICUD1', covars, outcome, time)
getRCSplot(covars = cox_risk_factors, 
           outcome='inhos_mortality',
           time='TZ_to_disch_duration')


getRCSplot(covars = tbl2_vars,
           x='input_output_ICUD1',
           outcome = 'AKI_01_23')

getRCSplot(covars = cox_risk_factors,
           x='input_output_ICUD1',
           outcome = 'inhos_mortality')

getRCSplot(covars = cox_risk_factors,
           x='input_output_ICUD1',
           outcome = 'icu_mortality')

getRCSplot(covars = tbl2_vars,
           x='input_output_ICUD2',
           outcome = 'AKI_01_23')

getRCSplot(covars = cox_risk_factors,
           x='input_output_ICUD2',
           outcome = 'inhos_mortality')

getRCSplot(covars = cox_risk_factors,
           x='input_output_ICUD2',
           outcome = 'icu_mortality')

getRCSplot(covars = tbl2_vars,
           x='input_output_ICUD3',
           outcome = 'AKI_01_23')

getRCSplot(covars = cox_risk_factors,
           x='input_output_ICUD3',
           outcome = 'inhos_mortality')

getRCSplot(covars = cox_risk_factors,
           x='input_output_ICUD3',
           outcome = 'icu_mortality')

df[,input_output_ICUD12 := input_output_ICUD1 + input_output_ICUD2]
getRCSplot(covars = tbl2_vars,
           x='input_output_ICUD12',
           outcome = 'AKI_01_23')

getRCSplot(covars = cox_risk_factors,
           x='input_output_ICUD12',
           outcome = 'inhos_mortality')

getRCSplot(covars = cox_risk_factors,
           x='input_output_ICUD12',
           outcome = 'icu_mortality')

df[,input_output_ICU0 := InputPreICU_ICUD1 - OutputPreICU_ICUD1]

getRCSplot(covars = tbl2_vars,
           x='input_output_ICU0',
           outcome = 'AKI_01_23')

getRCSplot(covars = cox_risk_factors,
           x='input_output_ICU0',
           outcome = 'inhos_mortality')

getRCSplot(covars = cox_risk_factors,
           x='input_output_ICU0',
           outcome = 'icu_mortality')

# Figure S1. Average Cr by ICU days by CKD vs Non-CKD ------------------
## line + ribbon plot
ribbon_line_plot = function(data, ckd){
  require(gridExtra)
  ckd = 0
  df_aki_melt = melt(df[Chronic_kidney_ds==ckd],
                      id.vars = c("SubjectNo", "AKI_stage","Chronic_kidney_ds"),
                      measure.vars = patterns("Cr_new_ICUD\\d|Cr_baseline"),
                      variable.name = "Day")
  df_aki_melt[,Day := fcase(
    Day %like% "baseline", "0",
    Day %like% "Cr_new_ICUD", sub("Cr_new_ICUD","",Day)
  ) |> as.factor()]
  df_aki_melt_aggr = df_aki_melt[,.(mean = mean(value, na.rm=T),
                                    sd = sd(value, na.rm=T),
                                    n=.N),by=.(Day,AKI_stage)]
  df_aki_melt_aggr[, m_sd_n := paste0(format(round(mean,1),nsmall=1)," ± ",
                                      format(round(sd,1),nsmall=1),'\n(n=',scales::comma(n),')')]
  p1 = df_aki_melt_aggr |> 
    ggplot(aes(x=Day, y=mean,
               group = as.factor(AKI_stage)))+
    geom_line(aes(color=as.factor(AKI_stage))) +
    geom_ribbon(aes(y=mean, ymin=mean-sd ,ymax = mean+sd, fill=as.factor(AKI_stage)),
                alpha=.2) +
    geom_point(aes(color=as.factor(AKI_stage))) +
    scale_y_continuous(expand=c(0,0),
                       limits=c(0,6),
                       breaks=seq(0,6,1))+
    scale_color_discrete(name="AKI Stage") + 
    scale_fill_discrete(name="AKI Stage") + 
    theme_classic() + 
    # facet_grid(~Chronic_kidney_ds,
    #            labeller = labeller(Chronic_kidney_ds=ckd_label))+
    labs(x="ICU Days",
         y= "Average Cr (mg/dL)") +
    theme(legend.position = "top",
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          axis.title = element_text(size=15),
          axis.text = element_text(size=15))
  # tbl = dcast(df_aki_melt_aggr, formula = AKI_stage ~ Day, value.var = 'm_sd_n')
  # p2 = ggtexttable(tbl,
  #             theme = ttheme(
  #               rownames.style = rownames_style(
  #                 face = 'plain'),
  #               colnames.style=colnames_style(
  #                 color = "black",
  #                 face = "bold",
  #                 size = 12,
  #                 fill = "white",
  #                 linewidth = 1,
  #                 linecolor = "white"),
  #               tbody.style = tbody_style(fill='white',  linecolor = "white")))
  # p1 + p2 +
  # #   plot_layout(ncol=1,heights = c(2.5,1))
  # p2 = tableGrob(tbl, rows = NULL)
  # p2$widths = unit(rep(1, ncol(p2)), "null")
  # p2$heights = unit(rep(1, nrow(p2)), "null")
  # p3 = ggplot() +
  #   annotation_custom(p2)
  # p1 + p3 + plot_layout(ncol=1, heights = c(2, 1))
  p1
}
ribbon_line_plot(data=df,ckd=0)
ribbon_line_plot(data=df,ckd=1)
dcast(df_aki_melt_aggr,
      formula = AKI_stage ~ Day, value.var = 'n') |> 
  ggtexttable()

df[Chronic_kidney_ds==1 & AKI_stage =='Non-AKI'][!is.na(Cr_new_initial),.N]
df[Chronic_kidney_ds==1 & AKI_stage =='Non-AKI'][!is.na(Cr_new_ICUD1),.N]
df[Chronic_kidney_ds==1 & AKI_stage =='Non-AKI'][!is.na(Cr_new_ICUD2),.N]
df[Chronic_kidney_ds==1 & AKI_stage =='Non-AKI'][!is.na(Cr_new_ICUD3),.N]
df[Chronic_kidney_ds==1 & AKI_stage =='Non-AKI'][!is.na(Cr_new_ICUD7),.N]

df_aki_melt_aggr
df[Chronic_kidney_ds==1 & !is.na(Cr)] |> dim()
df[Chronic_kidney_ds==1 & !is.na(Cr_ICUD1)] |> dim()
df[Chronic_kidney_ds==1 & !is.na(Cr_ICUD2)] |> dim()
df[Chronic_kidney_ds==1 & !is.na(Cr_ICUD3)] |> dim()
df[Chronic_kidney_ds==1 & !is.na(Cr_ICUD7)] |> dim()

ribbon_line_plot(df)
ribbon_line_plot(data = df[inhos_mortality==0])
ribbon_line_plot(data = df[inhos_mortality==1])
ribbon_line_plot(data = df[TypeOfInfection=='community'])
ribbon_line_plot(data = df[TypeOfInfection=='nosocomial'])


# Figure S2. Alluvial -----------------------------------------------------
require(ggalluvial)

df[,.(ICUADMDateTime_ICUD1, ICUDischDateTime, AKI_Day7)]
df[is.na(AKI_Day7),.(ICUADMDateTime_ICUD1, ICUDischDateTime)]

df[,AKI_initial_factor := factor(as.character(AKI_initial), levels=c("3","2","1","0"))]

df[,.SD,.SDcols=patterns("_status")]

df[,.N,AKI_Day1_status]
df[,summary(icu_duration)]
require(alluvial)
df_alluvial = df[AKI_status_day1 != 0 &
                   AKI_status_day2 != 0 &
                   AKI_status_day3 != 0 &
                   AKI_status_day7 != 0
                   ,.N, by=.(AKI_stage, AKI_initial, AKI_status_day1, AKI_status_day2, 
                           AKI_status_day3, AKI_status_day7)]
df_alluvial |> View()
alluvial(TZ=df_alluvial$AKI_initial, 
         `Day 1` = df_alluvial$AKI_status_day1, 
         `Day 2` = df_alluvial$AKI_status_day2, 
         `Day 3` = df_alluvial$AKI_status_day3, 
         `Day 7` = df_alluvial$AKI_status_day7,
         freq=df_alluvial$N,
         # border =ifelse(df_alluvial$AKI_stage == 'Stage 1','darkgreen',
         #                ifelse(df_alluvial$AKI_stage == 'Stage 2','blue','red'))
         col = ifelse(df_alluvial$AKI_stage == 'Stage 1','#00AFBB',
                             ifelse(df_alluvial$AKI_stage == 'Stage 2','#E7B800',"#FC4E07")),
         alpha=.8
         )
df_alluvial[,sum(N)]
df_alluvial[AKI_stage=='Stage 1', sum(N)]
df_alluvial[AKI_stage == 'Stage 1' & AKI_status_day7==1]
df_alluvial[,.N,AKI_stage]
df[,.N,AKI_YN]
# Figure S3A ----------------------------------------------------------------

# Table3 cox_risk_factors 변수 사용
form = paste0("inhos_mortality~",paste0(cox_risk_factors,collapse = "+"),"+(1|Center2)")
lr_fit = glmer(as.formula(form), family=binomial, data=df)
single_rows = c("SEX","s_pulmonary","s_abdominal","s_urinary","s_other",
                "s_unclear","fung_blood_sepsis","AppInitEmpThe_re","GramNegYN")
tbl_regression(
  x = lr_fit,
  exponentiate = T,
  # show_single_row = single_rows,
  pvalue_fun = ~style_pvalue(., digits=3),
  estimate_fun = ~style_ratio(., digits=2)
) |> 
  modify_table_styling(
    column = estimate,
    rows = !is.na(estimate),
    cols_merge_pattern = "{estimate} ({conf.low}-{conf.high})"
  ) |> 
  modify_header(estimate ~ "**OR (95% CI)**") |> 
  modify_column_hide(c(ci))

# Y: inhos_mortality
forest_df <- data.table(
  est = c("No AKI","Stage 1","Stage 2", "Stage 3", rep(1, 20)),
  hr = c(1.00, 0.95, 1.26, 2.49, rep(1,20)),
  ci_low = c(1.00, 0.73, 0.99, 2.05, rep(1,20)),
  ci_high = c(1.00, 1.24, 1.60, 3.03,rep(1,20))
)

# Y: icu_mortality
forest_df <- data.table(
  est = c("No AKI","Stage 1","Stage 2", "Stage 3", rep(1, 20)),
  hr = c(1.00, 0.95, 1.18, 1.59, rep(1,20)),
  ci_low = c(1.00, 0.78, 0.98, 1.38, rep(1,20)),
  ci_high = c(1.00, 1.19, 1.42, 1.83,rep(1,20))
)

require(forester)
forester(
  left_side_data = forest_df[,.(`AKI status` = est)],
  ci_low = forest_df$ci_low,
  ci_high = forest_df$ci_high,
  estimate=forest_df$hr,
  estimate_col_name = "     OR (95% CI)",
  font_family = "Sans",
  ci_sep = " — ",
  # display = F,
  # file_path = "~/Downloads/forestplot.png",
  estimate_precision = 2,
  null_line_at = 1,
  xlim= c(0.5, 3)
)
minmax = function(x){
  return((x-min(x))/(max(x)-min(x)))
}


df[,inhos_duration_minmax := minmax(inhos_duration)]
df[,icu_duration_minmax := minmax(icu_duration)]
df[,TZ_to_disch_duration_minmax := minmax(TZ_to_disch_duration)]

form = paste0("inhos_duration_minmax~",paste0(cox_risk_factors,collapse = "+"))
glm_fit = glm(as.formula(form), family=binomial, data=df)
forest_df <- data.table(
  est = c("No AKI","Stage 1","Stage 2", "Stage 3", rep(1, 20)),
  hr = c(1.00, 1.06, 1.36, 1.88, rep(1,20)),
  ci_low = c(1.00, 0.81, 1.07, 1.58, rep(1,20)),
  ci_high = c(1.00, 1.39, 1.71, 2.25,rep(1,20))
)


# figure S3B ----------------------------------------------------------------
form = paste0("inhos_mortality~",paste0(cox_risk_factors,collapse = "+"))
lr_fit = glm(as.formula(form), family=binomial, data=df)
single_rows = c("SEX","s_pulmonary","s_abdominal","s_urinary","s_other",
                "s_unclear","fung_blood_sepsis","AppInitEmpThe_re","GramNegYN")
tbl_regression(
  x = lr_fit,
  exponentiate = T,
  # show_single_row = single_rows,
  pvalue_fun = ~style_pvalue(., digits=3),
  estimate_fun = ~style_ratio(., digits=2)
) |> 
  modify_table_styling(
    column = estimate,
    rows = !is.na(estimate),
    cols_merge_pattern = "{estimate} ({conf.low}-{conf.high})"
  ) |> 
  modify_header(estimate ~ "**OR (95% CI)**") |> 
  modify_column_hide(c(ci))


forest_df <- data.table(
  est = c("No AKI","Stage 1","Stage 2", "Stage 3", rep(1, 20)),
  hr = c(1.00, 0.94, 1.16, 2.62, rep(1,20)),
  ci_low = c(1.00, 0.72, 0.91, 2.19, rep(1,20)),
  ci_high = c(1.00, 1.22, 1.48, 3.15,rep(1,20))
)
# remotes::install_github("rdboyes/forester")

require(forester)
forester(
  left_side_data = forest_df[,.(`AKI status` = est)],
  ci_low = forest_df$ci_low,
  ci_high = forest_df$ci_high,
  estimate=forest_df$hr,
  estimate_col_name = "     OR (95% CI)",
  font_family = "Sans",
  ci_sep = " — ",
  # display = F,
  # file_path = "~/Downloads/forestplot.png",
  estimate_precision = 2,
  null_line_at = 1,
  xlim= c(0.5, 3)
)


plot(p)
# Figure S4. AKI stage Criteria  ----------------------------------------------------------------

df[,urine_output_ICUD1 := (Output_ICUD1 / Wt / 24)]
df[,urine_output_ICUD2 := (Output_ICUD2 / Wt / 24)]
df[,urine_output_ICUD3 := (Output_ICUD3 / Wt / 24)]
df[,urine_output_ICUD7 := (Output_ICUD7 / Wt / 24)]

df[,AKI_Day1_urine := ifelse(urine_output_ICUD1 < 0.3, 1,0)]
df[,AKI_Day2_urine := ifelse(urine_output_ICUD2 < 0.3, 1,0)]
df[,AKI_Day3_urine := ifelse(urine_output_ICUD3 < 0.3, 1,0)]
df[,AKI_Day7_urine := ifelse(urine_output_ICUD7 < 0.3, 1,0)]

df[,AKI_Day0_cr := ifelse(Cr >= Cr_baseline_MDRD * 3.0 | Cr > 4, 1,0)]
df[,AKI_Day1_cr := ifelse(Cr_ICUD1 >= Cr_baseline_MDRD * 3.0 | Cr_ICUD1 > 4, 1,0)]
df[,AKI_Day2_cr := ifelse(Cr_ICUD2 >= Cr_baseline_MDRD * 3.0 | Cr_ICUD2 > 4, 1,0)]
df[,AKI_Day3_cr := ifelse(Cr_ICUD3 >= Cr_baseline_MDRD * 3.0 | Cr_ICUD3 > 4, 1,0)]
df[,AKI_Day7_cr := ifelse(Cr_ICUD7 >= Cr_baseline_MDRD * 3.0 | Cr_ICUD7 > 4, 1,0)]

df[,AKI_stage3_Cr := do.call(pmax, c(.SD, na.rm=T)),.SDcols=patterns("AKI_Day\\d_cr")]
df[,AKI_stage3_urine := do.call(pmax, c(.SD, na.rm=T)),.SDcols=patterns("_urine$")]
df[,AKI_stage3_Cr_urine := fifelse(AKI_stage=='Stage 3'  & AKI_stage3_urine==1,1,0)]

df[,AKI_stage3_crrt := do.call(pmax, c(.SD, na.rm=T)), .SDcols=patterns("^ICUStayCCRT_ICUD\\d")]
df[,AKI_stage3_yn := fifelse(AKI_stage == 'Stage 3',1,0)]

df[,`:=`(
  AKI_stage3_Cr_new = ifelse(AKI_stage3_Cr ==1 & AKI_stage3_urine ==0,1,0),
  AKI_stage3_urine_new = ifelse(AKI_stage3_Cr == 0 & AKI_stage3_urine ==1,1,0),
  AKI_stage3_Cr_urine_new = ifelse(AKI_stage3_Cr ==1 & AKI_stage3_urine ==1,1,0),
  AKI_stage3_crrt_new = ifelse(AKI_stage3_Cr ==0 & AKI_stage3_urine ==0 & AKI_stage3_crrt == 1,1,0)
)]

df[,.N,AKI_stage3_Cr_new]
df[,.N,AKI_stage3_urine_new]
df[,.N,AKI_stage3_Cr_urine_new]
df[,.N,AKI_stage3_crrt_new]
df[,.N,AKI_stage3_new]
# df[,AKI_stage_Cr_CRRT_UOP := ifelse(AKI_stage3_urine==1,'Stage 3',AKI_stage)]
# 
# df[AKI_stage_Cr_CRRT_UOP=='Stage 3' & (
#   ICUStayCCRT_ICUD1==1 | ICUStayCCRT_ICUD2 == 1|
#     ICUStayCCRT_ICUD3==1 | ICUStayCCRT_ICUD7 ==1
# ),.N]
# 
# 

urine_cr_stage3_melt = df[,.N, .(AKI_stage3_Cr_new, AKI_stage3_urine_new, 
                                 AKI_stage3_Cr_urine_new, 
                                 AKI_stage3_crrt_new)]
# urine_cr_stage3_melt_prop = urine_cr_stage3_melt[,.N,.(variable,value)][,prop := round(N/sum(N)*100,1), by=.(variable)]
# urine_cr_stage3_melt_prop[,name := fcase(variable=="AKI_stage3_yn","AKI by SCr criteria",
#                                          variable == 'AKI_stage3_urine', "AKI by UOP criteria",
#                                          variable == 'AKI_stage3_Cr_urine', "AKI by SCr & UOP criteria",
#                                          default = 'AKI by CRRT criteria') |> 
#                             factor(levels=c('AKI by SCr criteria',
#                                             'AKI by UOP criteria','AKI by CRRT criteria',
                                            # 'AKI by SCr & UOP criteria'))]

stage3_df = data.table(
  criteria = rep(c("Scr",'UOP','Scr & UOP','CRRT'), each=2),
  n=c(4216, 650, 3649, 1217, 4196, 670, 4617, 249),
  yn = rep(c(0,1),4),
  x=rep(1,8)
)
stage3_df[,prop := n/sum(n)*100, by=criteria][]
stage3_df[,criteria2 := factor(criteria, levels=c('UOP','Scr & UOP','Scr', 'CRRT'))]
stage3_df[yn==1]|> 
  ggplot(aes(x=reorder(x,-prop), y=prop, fill=criteria2)) +
  geom_col(position=position_dodge(),width=.5) +
  scale_x_discrete(label="AKI Stage 3") +
  scale_y_continuous(expand=c(0,0),
                     limits = c(0,30),
                     breaks = seq(0,30,10))+
  scale_fill_grey() +
  labs(y="Percentile") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=12))+
  geom_label(aes(label="p<0.001",y=28), fill="white", size=8)
  
matrix(stage3_df$n, nrow=2) |> chisq.test()
mat_test = matrix(stage3_df$n, ncol=2, byrow=T)
rownames(mat_test) = c('Scr','UOP','Scr & UOP','CRRT')
colnames(mat_test) = c(0,1)
chisq.test(mat_test)
require(rcompanion)
pairwiseNominalMatrix(mat_test,compare = "row")

df[AKI_stage3_Cr_new==1,.N,inhos_mortality]
df[AKI_stage3_urine_new==1,.N,inhos_mortality]
df[AKI_stage3_Cr_urine_new==1,.N,inhos_mortality]
df[AKI_stage3_crrt_new==1,.N,inhos_mortality]


# dt_death_by_def_melt = melt(df, id.vars='inhos_mortality',
#                             measure.vars = patterns('AKI_stage3_'))
# mytable(variable ~ inhos_mortality, dt_death_by_def_melt[value==1])
# dt_aki_death = dt_death_by_def_melt[value==1, .N, .(inhos_mortality, variable)][,prop:= N/sum(N)*100, by=variable][]
# dt_aki_death[, name := fcase(variable=="AKI_stage3_yn","AKI by SCr criteria",
#                              variable == 'AKI_stage3_urine', "AKI by UOP criteria",
#                              variable == 'AKI_stage3_Cr_urine', "AKI by SCr & UOP criteria",
#                              default = 'AKI by CRRT criteria') |> 
#                factor(levels=c('AKI by SCr & UOP criteria', 'AKI by CRRT criteria', 'AKI by SCr criteria',
#                                'AKI by UOP criteria'))]

stage3_inhos_df = data.table(
  criteria = rep(c("Scr",'UOP','Scr & UOP','CRRT'), each=2),
  n=c(439, 211, 684, 533, 292, 378, 98, 151),
  yn = rep(c(0,1),4),
  x=rep(1,8)
)
stage3_inhos_df[,prop := n/sum(n)*100, by=criteria][]
stage3_inhos_df[,criteria2 := factor(criteria, 
                                     levels=c('CRRT','Scr & UOP', 'UOP','Scr'))]

stage3_inhos_df[yn==1] |> 
  ggplot(aes(x=as.character(x), y=prop, fill=criteria2))+
  geom_col(position = position_dodge(), width=.5) +
  scale_y_continuous(limits = c(0,70), expand = c(0,0),
                     breaks = seq(0,70,10)) + 
  theme_classic() + 
  scale_x_discrete(labels='AKI Stage 3') +
  scale_fill_grey() +
  labs(y="Mortality percentile") +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.title.x = element_blank(),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=12)) +
  geom_label(aes(label="p<0.001", y=65), fill="white", size=8)

mat_death = matrix(c(439, 211, 684, 533, 292, 378, 98, 151),
                   ncol=2, byrow=T); mat_death 
mat_death |> chisq.test()
rownames(mat_death) = c('CRRT','Scr & UOP', 'UOP','Scr')
colnames(mat_death) = c(0,1)
pairwiseNominalMatrix(mat_death,compare = 'row',
                      method='bonferroni')


# AKI group In-hospital mortality (Kaplan-Meier, KM) ---------
require(survival)
require(survminer)
names(df)
load("AKI.RData")

temp = pairwise_survdiff(Surv(inhos_duration, inhos_mortality)~AKI_stage, df, p.adjust.method="bonferroni")
temp$p.value |> as.data.table(keep.rownames = T) |> gt()
pairwise_survdiff(Surv(Hospital_LOS, inhos_mortality)~AKI_stage_lowest, df, p.adjust.method="bonferroni")





df_aki_melt_lowest[,.(mean = mean(value, na.rm=T),
                sd = sd(value, na.rm=T)),by=.(Day,AKI_stage_lowest,Chronic_kidney_ds)]|> 
  ggplot(aes(x=Day, y=mean,
             group = as.factor(AKI_stage_lowest)))+
  geom_line(aes(color=as.factor(AKI_stage_lowest))) +
  geom_ribbon(aes(y=mean, ymin=mean-sd ,ymax = mean+sd, fill=as.factor(AKI_stage_lowest)),
              alpha=.2) +
  geom_point(aes(color=as.factor(AKI_stage_lowest))) +
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

## boxplot
ggplot(data=df_aki_melt,
       aes(x=as.factor(Day),y=value))+
  geom_boxplot(aes(fill=as.factor(Day))) +
  theme_classic() + 
  scale_fill_brewer(palette = "Pastel1")+
  facet_grid(~Chronic_kidney_ds,
             labeller = labeller(Chronic_kidney_ds=ckd_label)) +
  theme(legend.position = "none",
        strip.text.x = element_text(size=15,
                                    face = "bold")
  )+
  labs(x="Day", y="Value")

# inhos mortality  by AKI stage ----------------------
library(ggsci)
library(ggthemes)
barplot_inhos_by_aki_stage = function(aki){
  temp = df[,.N, keyby=c(aki, "inhos_mortality")][,prob := round(N/sum(N)*100,2), by=aki]
  temp[inhos_mortality==1] |> 
    ggplot(aes_string(x=aki, y= "prob", fill=aki))+
    geom_col() +
    scale_fill_brewer(palette = "Reds") +
    scale_y_continuous(limits=c(0,60), breaks=seq(0,60,10))+
    geom_text(aes(label=prob, y = prob +2), size=5)+
    theme_few() +
    labs(x="AKI stage", y="In-hospital mortality (%)") +
    theme(legend.position="none",
          axis.text = element_text(size=12),
          axis.title = element_text(size=15))
}

barplot_inhos_by_aki_stage("AKI_stage_lowest")

# bundle compliance by AKI stage ------------------------------

require(ggplot2)
require(stringr)
bundles <- df[,.(lac1h, bd1h, anti1h, br1h, app1h,
      lac3h, bd3h, anti3h, br3h, app3h,
      lac6h, bd6h, anti6h, br6h, app6h)] |> names()

# 전반적으로 br의 흐름이 제일 잘 나옴.
lapply(bundles[order(bundles)], \(x) {
  hour <- str_extract(x,"[:digit:]{1}[:alpha:]{1}")
  if(x %like% "lac") title = paste0("Lactate ", hour)
  else if (x %like% "app") title = paste0("Vasopressor ", hour)
  else if (x %like% "bd") title = paste0("Blood culture ", hour)
  else if (x %like% "br") title = paste0("Fluid resuscitation ", hour)
  else if (x %like% "anti") title = paste0("Antibiotic ", hour)
  temp <- df[,.N,by=.(AKI_stage, drug=var), env=list(var=x)][,prop:=N/sum(N)*100,by=AKI_stage]
  all <- temp[, sum(N), by=AKI_stage]$V1
  event <- temp[drug==1,sum(N), by=AKI_stage]$V1
  pval <- prop.trend.test(event, all)$p.value
  pval <- ifelse(pval < 0.001, "<0.001", format(round(pval,4),4))
  trend_label <- paste0("P for trend ",pval)
  temp[drug==1] |> 
    ggplot(aes(x=AKI_stage, y=prop, fill=AKI_stage))+
      geom_col(color="black") +
      ggtitle(title) +
      scale_y_continuous(limits=c(0,105), breaks=seq(0,100,20)) +
      geom_text(aes(y=prop+2, label=format(round(prop,1), nsmall=1)))+
      annotate("text", x=4, y=100, label=trend_label) +
      theme_classic()+
      labs(x="AKI stage", y="%") +
      scale_fill_brewer(palette="Blues") +
      theme(legend.position="none")
    ggsave(paste0("AKI/plots/",title,".png"),
          width= 150,
          height=120,
          units="mm",
          dpi=300
          )
})

df[,.N,by=.(AKI_stage, drug=var), env=list(var=x)][,prop:=N/sum(N)*100,by=AKI_stage][] 

all |> length()
event |> length()
prop.trend.test(all, event)


## input output by AKI stage -------------------------------
df[,.N,by=AKI_stage]
levels(df$AKI_stage)
input_output <- names(df)[names(df) %like% "Input|Output" & !names(df) %like% "ICUDL|PreICU"]
df_input_output_melt <- df[,c(input_output, "AKI_stage"), with=F]  |> 
  melt(id.vars="AKI_stage")

input_output_plot(data=df_input_output_melt,  target="input")
input_output_plot(data=df_input_output_melt,  target="output")



save.image("AKI.Rdata")


# sankey Alluvial AKI vs AKI lowest---------
library(ggalluvial)
dim(df[Chronic_kidney_ds==0])
df[Chronic_kidney_ds==0,.N, by=.(AKI_stage, AKI_stage_lowest)]
df[Chronic_kidney_ds==0, .N, AKI_stage]
df[Chronic_kidney_ds==0, .N, AKI_stage_lowest]
df[Chronic_kidney_ds==0,.N, by=.(AKI_stage, AKI_stage_lowest)][order(AKI_stage)] |>
    ggplot(aes(axis1=AKI_stage,
              weight = N,
               axis2 = AKI_stage_lowest, y = N)) + 
    geom_alluvium(aes(fill=AKI_stage)) +
    geom_stratum()+
    geom_text(stat="stratum",
              aes(label = paste0(after_stat(stratum),"\n n=", scales::comma(after_stat(n)))),
              size = 5) + 
    scale_x_discrete(limits = c("AKI_stage","AKI_stage_lowest"), expand = c(0.15,0.05))+
    annotate("text", label = "MDRD", x=1, y=4300, size=5)+
    annotate("text", label = "Lowest", x=2, y=4300, size=5) +
    theme_void() +
    scale_fill_viridis_d(direction=-1, option="viridis")+
    theme(legend.position = "none")

# TZ AKI stage 기준

# MAX AKI stage 기준
df[,AKI_stage_factor := factor(as.character(AKI_stage), 
                               levels=c('Stage 3','Stage 2', 'Stage 1','Non-AKI'))]
df[!is.na(AKI_Day7),.N, by=.(AKI_stage_factor, AKI_initial_factor, AKI_Day1, AKI_Day2, AKI_Day3, AKI_Day7)] |> 
  ggplot(aes(axis1= AKI_stage_factor,
             weight = N,
             axis2 = AKI_initial_factor,
             axis3 = AKI_Day1,
             axis4 = AKI_Day2,
             axis5 = AKI_Day3,
             axis6 = AKI_Day7,
             y=N)) + 
  geom_alluvium(aes(fill=AKI_stage_factor),
                reverse=T) +
  geom_stratum() + 
  geom_text(stat="stratum",
            aes(label = after_stat(stratum)),
            # aes(label = paste0(after_stat(stratum),"\n n=", scales::comma(after_stat(n)))),
            size = 5) +
  scale_x_discrete(limits = c('AKI max stage', "TZ","ICU Day1", "ICU Day2",
                              'ICU Day3','ICU Day7'), expand = c(0.15,0.05)) +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_fill_manual(values=c('red', 'orange','yellow','lightgreen'))+
  # scale_fill_viridis_d(direction=-1, option="viridis")+
  theme(legend.position = "none")


df[AKI_stage=='Stage 1',.N, by=AKI_Day7==0][!is.na(AKI_Day7),.(AKI_Day7, N,N/sum(N)*100)]
df[AKI_stage=='Stage 2',.N, by=AKI_Day7==0][!is.na(AKI_Day7),.(AKI_Day7, N,N/sum(N)*100)]
df[AKI_stage=='Stage 2',.N, by=AKI_Day7==1][!is.na(AKI_Day7),.(AKI_Day7, N,N/sum(N)*100)]
df[AKI_stage=='Stage 3',.N, by=AKI_Day7==0][!is.na(AKI_Day7),.(AKI_Day7, N,N/sum(N)*100)]
df[AKI_stage=='Stage 3',.N, by=AKI_Day7==1][!is.na(AKI_Day7),.(AKI_Day7, N,N/sum(N)*100)]

