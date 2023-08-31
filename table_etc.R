# baseline: by AKI stage -------------------------
require(moonBook)
require(gtsummary)
require(gt)

makeBaselineTable(data = df,
                  include=table_vars,
                  dicho_vars = dichotomous_vars,
                  by="AKI_stage")

makeBaselineTable(data = df,
                  include=table_vars,
                  dicho_vars = dichotomous_vars,
                  by="AKI_stage_lowest")

df[,..table_vars]

# baseline: AKI stage PostHoc ---------------------------------------------
df |> dim()

df[,table(AKI_stage, Comorbidity_MOSIAC_DM)]
baseline_pairwise = list() 
for (x in table_vars){
  if(class(df[[x]]) %in% c("numeric","integer") &
     df[,uniqueN(get(x))] > 2){
    baseline_pairwise[[x]] =  pairwise.t.test(x=df[[x]], g=df$AKI_stage, p.adjust.method = "bonferroni")$p.value |> 
      as.data.table(keep.rownames = T)
  } 
}
require(rcompanion)

cat_list = list()
for(x in table_vars){
  if(x %in% dichotomous_vars | uniqueN(df[[x]])==2 | x=="BacMTYN") {
    temp = rcompanion::pairwiseNominalMatrix(df[,table(AKI_stage, get(x))],
                                             fisher = F, gtest = F,
                                             chisq = T, method = "bonferroni",
                                             compare = "row",
                                             digits = 4)$Adjusted
    cat_list[[x]] = temp[upper.tri(temp)]
  }
}
names(cat_list)
do.call(rbind, cat_list) |> as.data.table(keep.rowname=T) |> gt()
temp = rcompanion::pairwiseNominalMatrix(df[,table(AKI_stage, SEX)],
                                         fisher = F, gtest = F,
                                         chisq = T, method = "bonferroni",
                                         compare = "row",
                                         digits = 4)$Adjusted
temp
temp[upper.tri(temp)]

temp2 = list()
for(i in names(baseline_pairwise)){
  temp2[[i]] = lapply(baseline_pairwise[[i]][,2:4], \(x)  na.omit(x)) |> 
    do.call(what="c")
}
do.call(cbind, temp2) |> t() |> as.data.table(keep.rownames=T) |> gt()



# baseline: CKD * AKI -----------------------------------------------------
df$RRT
df |> dim()
library(data.table)
makeBaselineTable(data=df[Chronic_kidney_ds==0],
                  include=table_vars,
                  dicho_vars = dichotomous_vars,
                  by="AKI_YN")

makeBaselineTable(data=df[Chronic_kidney_ds==1],
                  include=table_vars,
                  dicho_vars = dichotomous_vars,
                  by="AKI_YN")

makeBaselineTable(data=df[Chronic_kidney_ds==0],
                  include=table_vars,
                  dicho_vars = dichotomous_vars,
                  by="AKI_YN_lowest")

makeBaselineTable(data=df[Chronic_kidney_ds==1],
                  include=table_vars,
                  dicho_vars = dichotomous_vars,
                  by="AKI_YN")

# baseline: Dead vs Alive in AKI patients ---------------------------------
makeBaselineTable(data=df[AKI_YN==1],
                  by="inhos_mortality")

makeBaselineTable(data=df[AKI_YN==1],
                  by="icu_mortality")

# Early AKI vs Late AKI ---------------------------------------------------
# Early AKI: AKI within 48 hours after ICU admission
df[,AKI_early := fifelse(AKI_initial >= 1 & AKI_Day1>=1 | AKI_Day2 >= 1, 1,0) |> as.factor()]
df[,AKI_late := fifelse((AKI_initial == 0 & AKI_Day1 ==0  & AKI_Day2 ==0) & 
                          (AKI_Day3 >=1 | AKI_Day7 >=1),1,0) |> as.factor()]

df[,AKI_early := fifelse(do.call(pmax,.SD)>0,1,0 ),.SDcols=c("AKI_initial", "AKI_Day1", "AKI_Day2")]
df[AKI_early==1,.(AKI_initial, AKI_Day1, AKI_Day2)]
df[,AKI_late := fifelse(AKI_early==0 & do.call(pmax, c(.SD, na.rm=T))>0,1,0),.SDcols=c("AKI_Day3","AKI_Day7")]
df[AKI_early==0 & AKI_late==0,.N,AKI_initial]
df[,AKI_no_early_late := fcase(AKI_stage =="Non-AKI","Non_AKI",
                               AKI_early ==1 ,"Early_AKI",
                               AKI_late == 1, "Late_AKI")]
makeBaselineTable(data=df[AKI_no_early_late!="Non_AKI"],
                  by="AKI_no_early_late")

# baseline 2 survivor vs non-survivor---------------------------------------------------------------
require(data.table)
require(gtsummary)
getwd()
save.image(".RData")
makeBaselineTable(
  data=df,
  by = "inhos_mortality"
)

df[,.N,inhos_mortality]
df[,.N,icu_mortality]
tbl2 <- gtsummary::tbl_summary(
  data=df,
  by=icu_mortality,
  include = c(Age,SEX, BMI,
              Comorbidity_MOSIAC_DM,
              Comorbidity_MOSIAC_Cardio,
              Comorbidity_MOSIAC_Lung,
              Comorbidity_MOSIAC_CKD,
              Comorbidity_MOSIAC_SMT,
              Comorbidity_MOSIAC_HMM,
              Comorbidity_MOSIAC_IMM,
              Charlson_comorbidity_index_total, sofa_initial,
              IVS_SBP,IVS_DBP, IVS_MBP,IVS_HR,IVS_BT,
              Lactate, Hb, Plt, BUN, Cr,
              AST, ALT, Albumin, INR, CRP, Procalcitonin, ArterialPH,
              PaCO2, PaO2, Bicarbonate,
              s_pulmonary, s_abdominal, s_urinary, 
              s_skinsoft, s_other, s_unclear,
              PosBloodCulture,
              TypeOfInfection,
              GramPosYN, GramNegYN,
              Elig02,# com1h, com3h, com6h,
              AKI_stage,
              lac1h, bd1h, anti1h, br1h, app1h,
              lac3h, bd3h, anti3h, br3h, app3h,
              lac6h, bd6h, anti6h, br6h, app6h,
              ICUStayCCRT_OU, MV, ICU_LOS, Hospital_LOS
  ),
  missing = 'no',
  statistic = list(
    all_continuous() ~ "{mean} ± {sd}",
    all_categorical() ~ "{n} ({p})"
  ),
  type=list(
    c(SEX,s_pulmonary, s_abdominal,
      Elig02, PosBloodCulture,
      lac1h, bd1h, anti1h, br1h, app1h,
      lac3h, bd3h, anti3h, br3h, app3h,
      lac6h, bd6h, anti6h, br6h, app6h,
      s_urinary, s_skinsoft, s_other, 
      s_unclear, GramPosYN, GramNegYN)~ 'dichotomous'
  ),
  value = list(
    SEX ~ 'Male'
  ),
  digits = list(
    all_continuous() ~ 1,
    all_categorical() ~ c(0,1)
  )
) |> 
  add_overall() |> 
  add_p(pvalue_fun = ~style_pvalue(., digits=3)) 
tbl2 |> 
  modify_footnote(everything()~NA) |> 
  modify_table_body(~.x |> dplyr::relocate(stat_2, .after=stat_0))



# require(data.table)



# Cr Day by AKI stage ---------------------------------------------------------

df_aki_melt_lowest <- melt(df,
                           id.vars = c("SubjectNo", "AKI_stage_lowest","Chronic_kidney_ds"),
                           measure.vars = patterns("Cr_"),
                           variable.name = "Day")[Day !="Cr_baseline_MDRD" &
                                                    Day != "Cr_ICUDL"]                                         



df_aki_melt_lowest[,Day := fcase(
  Day %like% "baseline", "0",
  Day %like% "Cr_ICUD", sub("Cr_ICUD","",Day)
) |> as.factor()]


ckd_label <- c("Non-CKD","CKD")
names(ckd_label) <- c(0,1)
require(ggplot2)

## line + ribbon plot
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

require(moonBook)
library(ztable)
options(ztable.type="viewer")
df_aki_melt_lowest$Day
mytable(Day + AKI_stage ~ value, df_aki_melt[Chronic_kidney_ds==0])
mytable(Day + AKI_stage ~ value, df_aki_melt[Chronic_kidney_ds==1])
mytable(Day + AKI_stage_lowest ~ value, df_aki_melt_lowest[Chronic_kidney_ds==0]) |> ztable()




# RMANOVA -----------------------------------------------------------------
# One-way anova: 집단별로 차이 존재?
# 그룹 내 분산과 그룹 간 분산 차이 비교

# 가정
# 이상치 없음
# 정규성
# 구형성 가정 만족
require(ggpubr)
require(rstatix)
require(data.table)
df_aki_melt
res<-anova_test(data=df_aki_melt,
                dv=value, 
                wid=SubjectNo, 
                within=Day) 
get_anova_table(res)


one_way_rmanova <- aov(value~Day+Error(SubjectNo/Day),data=df_aki_melt)
one_way_rmanova |> summary()
# 
# two_way_rmanova <-aov(value~AKI_stage * Day+Error(id/Day),data=df_aki_melt) 
# two_way_rmanova |> summary()
box <- ggboxplot(df_aki_melt,
                 x="Day", y="value")

ggqqplot(df_aki_melt, "value", facet.by = "Day")

one_way_rma <- anova_test(data = df_aki_melt,
                          dv = value, wid = SubjectNo, within = Day,
)
one_way_rma
one_way_rma_post <- pairwise_t_test(
  data = df_aki_melt,
  formula = value ~ Day, paired = TRUE,
  p.adjust.method = "bonferroni"
)
one_way_rma_post <- one_way_rma_post |> 
  add_xy_position(x = "Day")

box + 
  stat_pvalue_manual(one_way_rma_post) +
  labs(
    subtitle = get_test_label(one_way_rma, detailed = TRUE),
    caption = get_pwc_label(one_way_rma_post)
  )

df_aki_melt |> 
  anova_test(value ~ AKI_stage * Day,
             wid= SubjectNo) |> 
  get_anova_table()

with(df_aki_melt,
     aov(value ~ AKI_stage * Day + Error(SubjectNo/AKI_stage * Day)))

# 집단 간: AKI_stage


# 집단 내: 일수
with(df_aki_melt,
     interaction.plot(Day, AKI_stage,value,type="b",
                      pch=c(2,4,6),legend=F, col=c(3,4,6),xlab="Time",ylab="Mean of value") 
)
legend("topleft",legend=c("Group 1","Group 2","Group 3"),pch=c(2,4,6),col=c(3,4,6), bg="gray90")


df_aki_melt |> 
  group_by(Day, Chronic_kidney_ds) |> 
  identify_outliers(value)

summary(fit)
df_aki_melt |> 
  rstatix::pairwise_t_test(formula = value ~ Day, paired = T, p.adjust.method = "bonferroni")
with(df_aki_melt, 
     pairwise.t.test(value,Day,paired=T,p.adjust.method="bonferroni"))

boxplot(value ~ Day*Chronic_kidney_ds,
        col=c("white","grey","grey30","grey10"),
        data=df_aki_melt)




# Linear Mixed Model ------------------------------------------------------

# value: Cr 값
# Day: ICU Day 1,2,3,7
# AKI stage: non,1,2,3

df_aki_melt[,.N,by=SubjectNo]
require(lme4)

fit <- lm(value ~ Day, data=df_aki_melt)
plot(fit,which = 2)


sleepstudy |> names()
lme4::
  # Association between value and day exist afoter controlling for
  # the variation in AKI_stage
  mixed_fit <- lmer(value ~ Day + (1 | AKI_stage:Chronic_kidney_ds), data=df_aki_melt)
mixed_fit |> summary()




# bundle LASSO------------------------------------------------------------------
df
require(glmnet)
require(data.table)
update_dev_pkg()
df_lasso <- df[,.(AKI_YN, 
                  lac1h, # lactate
                  bd1h, # 채혈
                  anti1h, #항생제
                  br1h, # 수액치료
                  app1h, # 승압제
                  com1h, # bundle 성공 여부
                  lac3h, bd3h, anti3h, br3h, app3h,com3h
)]

x <- df_lasso[,!"AKI_YN"] |> as.matrix()
y <- df_lasso[,"AKI_YN"] |> as.matrix()

lambdas <- seq(0, 0.3, by = .05)
glm_lasso <- glmnet(x = x, y = y,
                    family = "binomial",
                    alpha = 1)

cv_lasso <- cv.glmnet(
  x=x, y= y, alpha=1,
  family="binomial",
  nfolds = 10
)

cv_lasso
plot(glm_lasso)
plot(cv_lasso)

coefs <- coef(cv_lasso, s = cv_lasso$lambda.min)
require(data.table)
require(moonBook)

coefs_m <- as.matrix(coefs) |> data.table(keep.rownames = T)
vars <- coefs_m[s1!=0][2:11, rn]
form <- paste0("AKI_YN ~",paste0(vars,collapse = "+"))
glm(as.formula(form),
    family = "binomial",
    data=df_lasso) |> extractOR() |> ztable() |> ztable2flextable()
df_lasso[,.N,by=.(com1h)]
df_lasso[,.N,by=.(com3h)]

# Table 2-2 ---------------------------------------------------------------

cox_dichotomous_vars = df[,.(SEX, 
                             # Comorbidity_MOSIAC_DM,Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
                             # Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,
                             s_pulmonary, s_abdominal, s_urinary, 
                             # s_skinsoft, s_other, s_unclear, 
                             # GramPosYN, GramNegYN, 
                             fung_blood_sepsis, AppInitEmpThe_re,
                             # gp_ag_cl,
                             lac1h, bd1h, anti1h, br1h, app1h)] |> names()
# CRRT, vaso_inotrope_ICUD1

# age, sex, BMI, AKI stage, CCI, SOFA, SAPS3는 p value 에 관계없이 포함
must_include_vars = c("Age","SEX","BMI","Charlson_comorbidity_index_total",
                      "SOFA_ICUD1","SAPS3_ICUD1")
tbl2 = makeTable2(data=df, include = cox_risk_factors, dicho_vars = cox_dichotomous_vars,
                  by = "inhos_mortality")

tbl2_pval_vars = tbl2$table_body |> 
  dplyr::filter(p.value<=0.1) |> 
  select(variable) |> c()

mult_cox_vars = union(must_include_vars, tbl2_pval_vars$variable)
mult_cox_vars2 = setdiff(mult_cox_vars, "AKI_stage")
df_AKI = copy(df[AKI_YN==1])
df_NO_AKI = copy(df[AKI_YN==0])
df_AKI[,AKI_stage := relevel(AKI_stage, ref = "Stage 1")]
form = paste0("Surv(inhos_duration, inhos_mortality==1)~", paste0(cox_risk_factors, collapse = "+"))
cox_fit1 = coxph(as.formula(form),data=df[AKI_stage=="Stage 3"])
cox_fit2 = coxph(as.formula(form),data=df[AKI_stage!="Stage 3"])
cox_fit3 = coxph(as.formula(form),data=df_AKI)
cox_fit4 = coxph(as.formula(form),data=df_NO_AKI)
single_row_vars = c("SEX", "s_abdominal", "s_pulmonary", "s_urinary",
                    # "s_unclear", "s_skinsoft", "s_other", "GramPosYN",  "GramNegYN"
                    "fung_blood_sepsis","AppInitEmpThe_re")
tab1 = tbl2_cox_reg(cox_fit1, single_rows = single_row_vars)
tab2 = tbl2_cox_reg(cox_fit2, single_rows = single_row_vars)
tab3 = tbl2_cox_reg(cox_fit3, single_rows = single_row_vars)
tab4 = tbl2_cox_reg(cox_fit4, single_rows = single_row_vars)
tab3
tbl_merge(tbls = list(
  tab1, tab2, tab3, tab4
))


# Table 2-3 ---------------------------------------------------------------

cox_risk_factors = df[,.( Age,SEX, BMI,
                          Charlson_comorbidity_index_total, 
                          SOFA_ICUD1, SAPS3_ICUD1, 
                          # Lactate, Hb, Plt, Bilirubin, Albumin, CRP_imp, ArterialPH,
                          s_pulmonary, s_abdominal, s_urinary, s_skinsoft, s_other, s_unclear,
                          TypeOfInfection, 
                          bac_blood_sepsis_gp, bac_blood_sepsis_gn, fung_blood_sepsis,
                          AntbBEFCurrent,
                          AppInitEmpThe_re,  GramPosYN, GramNegYN, 
                          gp_ag_cl,
                          anti1h, br1h, app1h,
                          MV
                          # input_output_ICUD1, CRRT, vaso_inotrope_ICUD1 
)] |> names()
cox_dichotomous_vars = df[,.(SEX, 
                             # Comorbidity_MOSIAC_DM,Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
                             # Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,
                             s_pulmonary, s_abdominal, s_urinary, s_skinsoft, s_other, s_unclear, 
                             GramPosYN, GramNegYN, fung_blood_sepsis, AppInitEmpThe_re,
                             # gp_ag_cl,
                             lac1h, bd1h, anti1h, br1h, app1h,MV)] |> names()
# CRRT, vaso_inotrope_ICUD1

# age, sex, BMI, AKI stage, CCI, SOFA, SAPS3는 p value 에 관계없이 포함
must_include_vars = c("Age","SEX","BMI","Charlson_comorbidity_index_total",
                      "SOFA_ICUD1","SAPS3_ICUD1")
tbl2 = makeTable2(data=df, include = cox_risk_factors, dicho_vars = cox_dichotomous_vars,
                  by = "inhos_mortality")

tbl2_pval_vars = tbl2$table_body |> 
  dplyr::filter(p.value<=0.1) |> 
  select(variable) |> c()

mult_cox_vars = union(must_include_vars, tbl2_pval_vars$variable)
mult_cox_vars2 = setdiff(mult_cox_vars, "AKI_stage")
df_AKI = copy(df[AKI_YN==1])
df_NO_AKI = copy(df[AKI_YN==0])
df_AKI[,AKI_stage := relevel(AKI_stage, ref = "Stage 1")]
form = paste0("Surv(inhos_duration, inhos_mortality==1)~", paste0(cox_risk_factors, collapse = "+"))
cox_fit1 = coxph(as.formula(form),data=df[AKI_stage=="Stage 3"])
cox_fit2 = coxph(as.formula(form),data=df[AKI_stage!="Stage 3"])
cox_fit3 = coxph(as.formula(form),data=df_AKI)
cox_fit4 = coxph(as.formula(form),data=df_NO_AKI)
single_row_vars = c("SEX", "s_abdominal", "s_pulmonary", "s_urinary",
                    "s_unclear", "s_skinsoft", "s_other", "GramPosYN",
                    "fung_blood_sepsis","AppInitEmpThe_re",
                    "GramNegYN")
tab1 = tbl2_cox_reg(cox_fit1, single_rows = single_row_vars)
tab2 = tbl2_cox_reg(cox_fit2, single_rows = single_row_vars)
tab3 = tbl2_cox_reg(cox_fit3, single_rows = single_row_vars)
tab4 = tbl2_cox_reg(cox_fit4, single_rows = single_row_vars)

tbl_merge(tbls = list(
  tab1, tab2, tab3, tab4
))

# Table -----------------------------------------------------------------

form = paste0("Surv(inhos_duration, inhos_mortality==1)~",paste0(multi_cox_target$variable,collapse = "+"))
cox_fit = coxph(as.formula(form), data=df)
single_rows = c("SEX","s_pulmonary","s_abdominal","s_urinary","s_other",
                "s_unclear","fung_blood_sepsis","AppInitEmpThe_re","GramNegYN")

tbl_regression(
  x = cox_fit,
  exponentiate = T,
  show_single_row = single_rows,
  pvalue_fun = ~style_pvalue(., digits=3),
  estimate_fun = ~style_ratio(., digits=2)
) |> 
  modify_table_styling(
    column = estimate,
    rows = !is.na(estimate),
    cols_merge_pattern = "{estimate} ({conf.low}-{conf.high})"
  ) |> 
  modify_header(estimate ~ "**HR (95% CI)**") |> 
  modify_column_hide(c(ci))
