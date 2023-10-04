# Figure 2. KM plot----------------------------------------------------------------
aki = "AKI_stage"
km_inhos_by_aki_stage = function(aki){
  form = paste0("Surv(inhos_duration, inhos_mortality==1) ~ ", aki)
  # form = paste0("Surv(icu_duration, icu_mortality==1) ~ ", aki)
  inhos.aki.group.km_fit = survfit(as.formula(form), data=df)
  inhos_aki_group_km = ggsurvplot(inhos.aki.group.km_fit,
                                  xlim= c(0,100),
                                  break.time.by=10,
                                  censor=F,
                                  risk.table=T,
                                  linetype = c("dotted","dotdash","dashed", "solid"),
                                  pval="P < 0.0001",
                                  legend.labs=c("Without SA-AKI","SA-AKI Stage 1","SA-AKI Stage 2","SA-AKI Stage 3"),
                                  palette="jama",
                                  pval.coord = c(1, 0.2),
                                  xlab="Time (day)"
                                  # surv.median.line="hv",
  )
  inhos_aki_group_km$plot  = inhos_aki_group_km$plot + 
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

# Figure S1. Average Cr by ICU days by CKD vs Non-CKD ------------------
## line + ribbon plot
ribbon_line_plot = function(data, ckd){
  require(gridExtra)
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
  max_y = ifelse(ckd==0,6,8)
  p1 = df_aki_melt_aggr |> 
    ggplot(aes(x=Day, y=mean,
               group = as.factor(AKI_stage))) +
    geom_line(aes(color=as.factor(AKI_stage))) +
    geom_ribbon(aes(y=mean, ymin=mean-sd ,ymax = mean+sd, fill=as.factor(AKI_stage)),
                alpha=.2) +
    geom_point(aes(color=as.factor(AKI_stage))) +
    scale_x_discrete(labels=c("Time zero",'ICU D1', 'ICU D2','ICU D3','ICU D7')) +
    scale_y_continuous(expand=c(0,0),
                       limits=c(0,max_y),
                       breaks=seq(0,max_y,1)) +
    scale_color_discrete(name="SA-AKI Stage",
                         labels=c("Without SA-AKI","SA-AKI Stage 1","SA-AKI Stage 2", "SA-AKI Stage 3")) + 
    scale_fill_discrete(name="SA-AKI Stage",
                        labels=c("Without SA-AKI","SA-AKI Stage 1","SA-AKI Stage 2", "SA-AKI Stage 3")) + 
    theme_classic() + 
    labs(x="Time",
         y= "Mean serum creatinine (mg/dL)") +
    theme(legend.position = "top",
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          axis.title = element_text(size=15),
          axis.text = element_text(size=15))
  p1
}
ribbon_line_plot(data=df,ckd=0)
ribbon_line_plot(data=df,ckd=1)

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
require(alluvial)
df[,.(ICUADMDateTime_ICUD1, ICUDischDateTime, AKI_Day7)]
df[is.na(AKI_Day7),.(ICUADMDateTime_ICUD1, ICUDischDateTime)]

df[,AKI_initial_factor := factor(as.character(AKI_initial), levels=c("3","2","1","0"), labels = c("3","2","1","NoAKI"))]

df[,.SD,.SDcols=patterns("_status")]
AKI_status_day1 != 0 &
  AKI_status_day2 != 0 &
  AKI_status_day3 != 0 &
  AKI_status_day7 != 0
df_alluvial = df[AKI_stage!="Non-AKI",.N, by=.(AKI_stage, AKI_initial, AKI_status_day1, AKI_status_day2, 
                           AKI_status_day3, AKI_status_day7)]

target = df_alluvial[,2:6] |> names()
df_alluvial[,(target):=lapply(.SD, \(x) ifelse(x == 0, 'Without\nAKI', as.character(x))),.SDcols=target]
df_alluvial[,(target):=lapply(.SD, factor, levels=rev(c("Dead", "3","2","1","Without\nAKI","Alive"))),.SDcols = target]
df_alluvial
alluvial(`Time zero` = df_alluvial$AKI_initial, 
         "ICU D1" = df_alluvial$AKI_status_day1, 
         `ICU D2` = df_alluvial$AKI_status_day2, 
         `ICU D3` = df_alluvial$AKI_status_day3, 
         `ICU D7` = df_alluvial$AKI_status_day7,
         freq=df_alluvial$N,
         # border =ifelse(df_alluvial$AKI_stage == 'Stage 1','darkgreen',
         #                ifelse(df_alluvial$AKI_stage == 'Stage 2','blue','red'))
         col = ifelse(df_alluvial$AKI_stage == 'Stage 1','#00AFBB',
                             ifelse(df_alluvial$AKI_stage == 'Stage 2','#E7B800',"#FC4E07")),
         alpha=.8,
         axis_labels = c("Time zero", 'ICU D1','ICU D2','ICU D3','ICU D7')
         )
boxplot(mpg ~ cyl, data=mtcars)
legend(x="top", 
       legend=c("SA-AKI stage 3","SA-AKI stage 2","SA-AKI stage 1"), fill=c("#FC4E07",'#E7B800','#00AFBB'),
       horiz = T)
df_alluvial[,sum(N)]
df_alluvial[AKI_stage=='Stage 1', sum(N)]
df_alluvial[AKI_stage == 'Stage 1' & AKI_status_day7==1]

# Figure S3A: inhos mortality ----------------------------------------------------------------

# Table3 cox_risk_factors 변수 사용
form = paste0("inhos_mortality~",paste0(c("AKI_stage",cox_risk_factors),collapse = "+"),"+(1|Center2)")
form = paste0("icu_mortality~",paste0(c("AKI_stage",cox_risk_factors),collapse = "+"),"+(1|Center2)")
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
  modify_header(estimate ~ "**adjusted OR (95% CI)**") |> 
  modify_column_hide(c(ci))

# Y: inhos_mortality
forest_df = data.table(
  est = c("Without SA-AKI","SA-AKI Stage 1","SA-AKI Stage 2", "SA-AKI Stage 3", rep(1, 20)),
  hr = c(1.00, 0.96, 1.28, 2.52, rep(1,20)),
  ci_low = c(1.00, 0.74, 1.01, 2.07, rep(1,20)),
  ci_high = c(1.00, 1.26, 1.64, 3.07,rep(1,20))
)

# Y: ICU mortality
forest_df = data.table(
  est = c("Without SA-AKI","SA-AKI Stage 1","SA-AKI Stage 2", "SA-AKI Stage 3", rep(1, 20)),
  hr = c(1.00, 1.05, 1.42, 2.98, rep(1,20)),
  ci_low = c(1.00, 0.77, 1.07, 2.38, rep(1,20)),
  ci_high = c(1.00, 1.24, 1.88, 3.73,rep(1,20))
)

require(forester)
forester(
  left_side_data = forest_df[,.(`SA-AKI status` = est)],
  ci_low = forest_df$ci_low,
  ci_high = forest_df$ci_high,
  estimate=forest_df$hr,
  estimate_col_name = "Adjusted OR (95% CI)",
  font_family = "Sans",
  ci_sep = " — ",
  justify = c(0, 0.5),
  # display = F,
  # file_path = "~/Downloads/forestplot.png",
  estimate_precision = 2,
  null_line_at = 1,
  xlim= c(0.5, 3.5)
)

# figure S3B: ICU mortality ----------------------------------------------------------------
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


forest_df = data.table(
  est = c("Without SA-AKI","SA-AKI Stage 1","SA-AKI Stage 2", "SA-AKI Stage 3", rep(1, 20)),
  hr = c(1.00, 1.05, 1.42, 2.98, rep(1,20)),
  ci_low = c(1.00, 0.77, 1.07, 2.38, rep(1,20)),
  ci_high = c(1.00, 1.24, 1.88, 3.73,rep(1,20))
)

forester(
  left_side_data = forest_df[,.(`SA-AKI status` = est)],
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
  xlim= c(0.5, 3.5)
)

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
df[,AKI_stage3_urine := do.call(pmax, c(.SD, na.rm=T)),.SDcols=patterns("AKI_Day\\d_urine")]
df[,AKI_stage3_Cr_urine := fifelse(AKI_stage=='Stage 3'  & AKI_stage3_urine==1,1,0)]

df[,AKI_stage3_crrt := do.call(pmax, c(.SD, na.rm=T)), .SDcols=patterns("^ICUStayCCRT_ICUD\\d")]
df[,AKI_stage3_yn := fifelse(AKI_stage == 'Stage 3',1,0)]

df[,AKI_stage3_crrt_new := ifelse(AKI_stage3_crrt == 1,1,0)]
df[,`:=`(
  AKI_stage3_Cr_new = ifelse(AKI_stage3_crrt_new == 0 & AKI_stage3_Cr ==1 & AKI_stage3_urine ==0,1,0),
  AKI_stage3_urine_new = ifelse(AKI_stage3_crrt_new == 0 & AKI_stage3_Cr == 0 & AKI_stage3_urine ==1,1,0),
  AKI_stage3_Cr_urine_new = ifelse(AKI_stage3_crrt_new == 0 & AKI_stage3_Cr ==1 & AKI_stage3_urine ==1,1,0)
)]

df[,.N,AKI_stage3_Cr_new]
df[,.N,AKI_stage3_urine_new]
df[,.N,AKI_stage3_Cr_urine_new]
df[,.N,AKI_stage3_crrt_new]
df[,.N,AKI_stage3_new]
 

urine_cr_stage3_melt = df[,.N, .(AKI_stage3_Cr_new, AKI_stage3_urine_new, 
                                 AKI_stage3_Cr_urine_new, 
                                 AKI_stage3_crrt_new)]

stage3_df = data.table(
  criteria = rep(c("Scr",'UOP','Scr & UOP','CRRT'), each=2),
  n=c(4541, 321, 3939, 927, 4614, 252, 3576, 1290),
  yn = rep(c(0,1),4),
  x=rep(1,8)
)
stage3_df[,prop := n/sum(n)*100, by=criteria][]
stage3_df[,criteria2 := factor(criteria, levels=c('CRRT', 'UOP','Scr','Scr & UOP'))]
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
  n=c(255, 66, 612, 315, 136, 116, 510, 780),
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
  # geom_label(aes(label="p<0.001", y=65), fill="white", size=8) +
  geom_signif(
    y_position = c(63,49, 38),
    xmin= c(0.795, 0.925, 1.055),
    xmax = c(0.925, 1.055, 1.185),
    annotations = rep('p < 0.001',3),
    tip_length = .05
  )
stage3_inhos_df
fcount(df$AKI_stage)
df |> 
  sbt(AKI_stage3_urine_new==1 & AKI_stage=='Non-AKI') |> 
  get_vars('urine_output_ICUD\\d', regex=T)
  
get_vars(df, 'urine_output_ICUD\\d', regex=T) |> 
  descr()

df[AKI_status_day1!='Dead',urine_output_ICUD1] |> descr()
df[AKI_status_day2!='Dead',urine_output_ICUD2] |> descr()
df[AKI_status_day3!='Dead',urine_output_ICUD3] |> descr()
df[AKI_status_day7!='Dead',urine_output_ICUD7] |> descr()
df[AKI_Day1_urine==1 & AKI_status_day1=='Dead'] |> dim()
df[AKI_Day2_urine==1 & AKI_status_day2=='Dead'] |> dim()
df[AKI_Day3_urine==1 & AKI_status_day3=='Dead'] |> dim()
df[AKI_Day7_urine==1 & AKI_status_day7=='Dead'] |> dim()

df[,summary(.SD), .SDcols=patterns('urine_output_ICUD\\d')]

mat_death = matrix(stage3_inhos_df$n, ncol=2, byrow=T)
mat_death |> chisq.test()
rownames(mat_death) = c('Scr', 'UOP', 'Scr & UOP', 'CRRT')
colnames(mat_death) = c(0,1)
pairwiseNominalMatrix(mat_death,compare = 'row',
                      method='bonferroni')

# Figure S5. I/O Restricted cubic spline plot ----------------------------------------------------------------
library(rms)
library(plotRCS)

# X: input
# covars: ~ + output Day ~
# Outcomes: AKI_01_23, inhos_mortality, icu_mortality
getRCSplot = function(x, target=NA,covars, outcome, time=NULL, xlab){
  vars = c(x, covars, outcome, time)
  if(is.na(target)) data = df[,..vars] |> na.omit()
  else data = df[AKI_01_23 == target, ..vars] |> na.omit()
  rcsplot(data = data ,
          outcome = outcome,
          time = time,
          exposure = x,
          covariates = covars,
          fontfamily = "sans",
          pvalue.label.nonlinear = "P-value for nonlinear spline term",
          na.rm=T) +
    labs(x=xlab, y="adjusted odds ratio")
}

df[,input_ICU0123 := rowSums(.SD, na.rm=T),.SDcols=c('InputPreICU_ICUD1', 'Input_ICUD1','Input_ICUD2','Input_ICUD3')]
df[,output_ICUD0123 := rowSums(.SD, na.rm=T),.SDcols=c('OutputPreICU_ICUD1', 'Output_ICUD1','Output_ICUD2', 'Output_ICUD3')]

getRCSplot(covars = c(cox_risk_factors,'output_ICUD0123'),
           x='input_ICU0123',
           outcome = 'inhos_mortality',
           xlab = 'Cumulative dose of fluid intake')

getRCSplot(covars = c(cox_risk_factors,'output_ICUD0123'),
           x='input_ICU0123',
           outcome = 'icu_mortality',
           xlab = 'Cumulative dose of fluid intake')

getRCSplot(target=1,
  covars = c(cox_risk_factors,'output_ICUD0123'),
           x='input_ICU0123',
           outcome = 'inhos_mortality',
           xlab = 'Cumulative dose of fluid intake')

getRCSplot(target=1,
  covars = c(cox_risk_factors,'output_ICUD0123'),
           x='input_ICU0123',
           outcome = 'icu_mortality',
           xlab = 'Cumulative dose of fluid intake')

# df[,input_ICU012 := rowSums(.SD, na.rm=T),.SDcols=c('InputPreICU_ICUD1', 'Input_ICUD1','Input_ICUD2')]
# df[,output_ICUD012 := rowSums(.SD, na.rm=T),.SDcols=c('OutputPreICU_ICUD1', 'Output_ICUD1','Output_ICUD2')]
# getRCSplot(covars = c(tbl2_vars,'output_ICUD012'),
#            x='input_ICU012',
#            outcome = 'AKI_01_23',
#            xlab='Input before ICU Admission ~ ICU Day2')
# 
# getRCSplot(covars = c(cox_risk_factors,'output_ICUD012'),
#            x='input_ICU012',
#            outcome = 'inhos_mortality',
#            xlab = 'Input before ICU Admission ~ ICU Day2')
# 
# getRCSplot(covars = c(cox_risk_factors,'output_ICUD012'),
#            x='input_ICU012',
#            outcome = 'icu_mortality',
#            xlab = 'Input before ICU Admission ~ ICU Day2')
# 
# getRCSplot(covars = c(tbl2_vars,'output_ICUD0123'),
#            x='input_ICU0123',
#            outcome = 'AKI_01_23',
#            xlab='Input before ICU Admission ~ ICU Day3')




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

# inhos mortality  by AKI stage ----------------------
library(ggsci)
library(ggthemes)
barplot_inhos_by_aki_stage("AKI_stage_lowest")

# bundle compliance by AKI stage ------------------------------

require(ggplot2)
require(stringr)
bundles = df[,.(lac1h, bd1h, anti1h, br1h, app1h,
      lac3h, bd3h, anti3h, br3h, app3h,
      lac6h, bd6h, anti6h, br6h, app6h)] |> names()

# 전반적으로 br의 흐름이 제일 잘 나옴.
lapply(bundles[order(bundles)], \(x) {
  hour = str_extract(x,"[:digit:]{1}[:alpha:]{1}")
  if(x %like% "lac") title = paste0("Lactate ", hour)
  else if (x %like% "app") title = paste0("Vasopressor ", hour)
  else if (x %like% "bd") title = paste0("Blood culture ", hour)
  else if (x %like% "br") title = paste0("Fluid resuscitation ", hour)
  else if (x %like% "anti") title = paste0("Antibiotic ", hour)
  temp = df[,.N,by=.(AKI_stage, drug=var), env=list(var=x)][,prop:=N/sum(N)*100,by=AKI_stage]
  all = temp[, sum(N), by=AKI_stage]$V1
  event = temp[drug==1,sum(N), by=AKI_stage]$V1
  pval = prop.trend.test(event, all)$p.value
  pval = ifelse(pval < 0.001, "<0.001", format(round(pval,4),4))
  trend_label = paste0("P for trend ",pval)
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
input_output = names(df)[names(df) %like% "Input|Output" & !names(df) %like% "ICUDL|PreICU"]
df_input_output_melt = df[,c(input_output, "AKI_stage"), with=F]  |> 
  melt(id.vars="AKI_stage")

input_output_plot(data=df_input_output_melt,  target="input")
input_outrput_plot(data=df_input_output_melt,  target="output")

# MAX AKI stage 기준
# library(ggalluvial)
# df[,AKI_stage_factor := factor(as.character(AKI_stage), 
#                                levels=c('Stage 3','Stage 2', 'Stage 1','Non-AKI'))]
# df[!is.na(AKI_Day7),.N, by=.(AKI_stage_factor, AKI_initial_factor, AKI_Day1, AKI_Day2, AKI_Day3, AKI_Day7)] |> 
#   ggplot(aes(axis1= AKI_stage_factor,
#              weight = N,
#              axis2 = AKI_initial_factor,
#              axis3 = AKI_Day1,
#              axis4 = AKI_Day2,
#              axis5 = AKI_Day3,
#              axis6 = AKI_Day7,
#              y=N)) + 
#   geom_alluvium(aes(fill=AKI_stage_factor),
#                 reverse=T) +
#   geom_stratum() + 
#   geom_text(stat="stratum",
#             aes(label = after_stat(stratum)),
#             # aes(label = paste0(after_stat(stratum),"\n n=", scales::comma(after_stat(n)))),
#             size = 5) +
#   scale_x_discrete(limits = c('AKI max stage', "TZ","ICU Day1", "ICU Day2",
#                               'ICU Day3','ICU Day7'), expand = c(0.15,0.05)) +
#   theme_classic() +
#   theme(
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank()
#   ) +
#   scale_fill_manual(values=c('red', 'orange','yellow','lightgreen'))+
#   # scale_fill_viridis_d(direction=-1, option="viridis")+
#   theme(legend.position = "none")