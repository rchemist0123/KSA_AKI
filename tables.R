library(gtsummary)
# Table 1 ------------------------------------------------
table_vars = df[,.(Age,SEX, BMI,
                   # Comorbidities
                   Comorbidity_MOSIAC_DM, Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
                   Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,Charlson_comorbidity_index_total, 
                   CFScore,
                   # Source of infection
                   s_pulmonary, s_abdominal, s_urinary, s_skinsoft, s_other, s_unclear,
                   TypeOfInfection, bac_blood_sepsis_gp, bac_blood_sepsis_gn,
                   # Vital
                   IVS_MBP, IVS_HR, IVS_BT, 
                   # Lab
                   Lactate, CRP_imp, Procalcitonin_imp, 
                   # ArterialPH, PaCO2, PaO2, Bicarbonate,
                   sofa_initial, SOFA_Renal, SOFA_except_renal, 
                   hosp_to_icu_duration,
                   TZ_ICUAdm_gap_hour,
                   
                   # Characteristics of ICU D1 
                   SAPS3_ICUD1, Septic_shock_ICUD1,  ICUStayInvasive_ICUD1,
                   ICUStayInvasive_OU,ICUStayECMO_ICUD1,ICUStayECMO_OU)] |> names()

dichotomous_vars = df[,.(SEX, Comorbidity_MOSIAC_DM,Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
                         Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,
                         s_pulmonary, s_abdominal,
                         s_urinary, s_skinsoft, s_other, s_unclear, 
                         bac_blood_sepsis_gp, bac_blood_sepsis_gn, 
                         Septic_shock_ICUD1,  ICUStayInvasive_ICUD1,
                         ICUStayInvasive_OU,ICUStayECMO_ICUD1,ICUStayECMO_OU)] |> names()

makeBaselineTable(data = df,
                  include = table_vars,
                  dicho_vars = dichotomous_vars,
                  by = "AKI_YN")


# Table 2 -----------------------------------------------------------------

tbl2_vars = df[,.(Age,SEX, BMI, CFScore,
             Comorbidity_MOSIAC_DM, Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
             Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM, 
             sofa_initial, s_pulmonary, s_abdominal, s_urinary, 
             TypeOfInfection, bac_blood_sepsis_gp, bac_blood_sepsis_gn,
             AntbBEFCurrent2)] |> names()

dicho_vars = df[,.(SEX,Comorbidity_MOSIAC_DM, Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung, Comorbidity_MOSIAC_CKD,
                   Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,
                   # s_pulmonary, s_abdominal, s_urinary, s_skinsoft, s_other, s_unclear,
                   bac_blood_sepsis_gp, bac_blood_sepsis_gn,
                   AntbBEFCurrent2)] |> names()

library(lme4)
library(car)
library(coxme)
library(finalfit)
library(kableExtra)
tbl2 = finalfit.glm(
  .data=df,
  dependent = "AKI_01_23",
  explanatory = tbl2_vars,
  explanatory_multi = setdiff(tbl2_vars, 'CFScore'),
  keep_models = T,
  random_effect = "Center2"
)

tbl2 |> 
  kbl() |> 
  kable_paper("hover", full_width = F)

# For univariable
temp = lapply(tbl2_vars,\(x){
  form = paste0("AKI_01_23 ~ ",x,"+(1|Center2)")
  fit = glmer(as.formula(form), family=binomial, data=df)
  tab = tbl_regression(
    x=fit,
    exponentiate = T,
    pvalue_fun = ~style_pvalue(., digits=3),
    estimate_fun = ~style_ratio(., digits=2)
  )|> modify_table_styling(
      column = estimate,
      rows = !is.na(estimate),
      cols_merge_pattern = "{estimate} ({conf.low}-{conf.high})"
    ) |> 
    modify_header(estimate ~ "**OR (95% CI)**") |> 
    modify_column_hide(c(ci))
})

tbl_stack(temp)

# For multi
form = paste0("AKI_01_23~",
              paste0(setdiff(tbl2_vars,c("CFScore","Comorbidity_MOSIAC_HMM","Comorbidity_MOSIAC_IMM")),collapse = "+"),'+(1|Center2)')
lr_fit = glmer(as.formula(form),family=binomial, data=df)
summary(lr_fit)
vif(lr_fit)
tbl_regression(
  x = lr_fit,
  exponentiate = T,
  pvalue_fun = ~style_pvalue(., digits=3),
  estimate_fun = ~style_ratio(., digits=2),
  show_single_row = c(
    s_pulmonary, s_abdominal, s_urinary,
    AntbBEFCurrent2
  ),
) |> 
  modify_table_styling(
    column = estimate,
    rows = !is.na(estimate),
    cols_merge_pattern = "{estimate} ({conf.low}-{conf.high})"
  ) |> 
  modify_header(estimate ~ "**OR (95% CI)**") |> 
  modify_column_hide(c(ci))

# Table 3----------------------------------------------------------

# Bivariable

tbl_summary(
  df[AKI_01_23==1],
  by = inhos_mortality,
  include = cox_risk_factors,
  statistic = list(
    all_continuous() ~ "{mean} ± {sd}",
    all_categorical() ~ "{n} ({p})"
  ),
  digits = list(
    all_continuous() ~ 1,
    all_categorical() ~ c(0,1)
  ),
  type = list(
    c(CFScore) ~ 'continuous',
    c(SEX,Septic_shock_ICUD1, s_pulmonary, 
      s_abdominal, s_urinary,AntbBEFCurrent2, AppInitEmpThe_re) ~ "dichotomous"
  ),
  value = list(SEX ~ 'Male'),
  missing="no"
) |> 
  add_p(pvalue_fun = ~ style_pvalue(., digits=3))


cox_risk_factors = df[,.( Age,SEX, BMI, CFScore,
                          # AKI_stage, 
                          Comorbidity_MOSIAC_DM, Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
                          Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM, 
                          SOFA_ICUD1, SAPS3_ICUD1,  Septic_shock_ICUD1,
                          # Lactate, Hb, Plt, Bilirubin, Albumin, CRP_imp, ArterialPH,
                          s_pulmonary, s_abdominal, s_urinary,
                          # s_skinsoft, s_other, s_unclear,
                          TypeOfInfection, 
                          bac_blood_sepsis_gp, bac_blood_sepsis_gn,
                          AntbBEFCurrent2,
                          AppInitEmpThe_re, 
                          lac1h, bd1h, anti1h, br1h, app1h)] |> names()
# For univariable
table3_univ = function(status){
  stopifnot(is.numeric(status))
  temp = lapply(cox_risk_factors,\(x){
    form = paste0("inhos_mortality ~ ",x,"+(1|Center2)")
    fit = glmer(as.formula(form), family=binomial, data=df[AKI_01_23==status])
    tab = tbl_regression(
      x=fit,
      exponentiate = T,
      pvalue_fun = ~style_pvalue(., digits=3),
      estimate_fun = ~style_ratio(., digits=2)
    )|> modify_table_styling(
      column = estimate,
      rows = !is.na(estimate),
      cols_merge_pattern = "{estimate} ({conf.low}-{conf.high})"
    ) |> 
      modify_header(estimate ~ "**OR (95% CI)**") |> 
      modify_column_hide(c(ci))
  })
  
  tbl_stack(temp)
}
table3_univ(0)

# For multivariable
shapiro.test(df[AKI_01_23==0]$shock_index_SBP)
shapiro.test(df[AKI_01_23==1]$shock_index_SBP)
mytable(AKI_01_23 ~ shock_index_SBP + shock_index_DBP, df, digits=3)

form = paste0('inhos_mortality ~ ',paste0(cox_risk_factors, collapse = "+"),"+(1|Center2)")
fit = glmer(as.formula(form), family=binomial, data=df[AKI_01_23==0])
makeTable3(fit)
fit = glmer(as.formula(form), family=binomial, data=df[AKI_01_23==1])
makeTable3Reg(fit)

tbl_summary(
  data=df,
  include = c(shock_index_SBP, shock_index_DBP),
  by=AKI_01_23,
  digits=everything() ~ 3,
  
  missing='no'
) |> add_p(pvalue_fun = ~style_pvalue(., digit=3))


#AKI
a=makeTable3(fit)
a

# check for vif
form = paste0('inhos_mortality ~ ',paste0(setdiff(cox_risk_factors,'AKI_stage'), collapse = "+"),"+(1|Center2)")
fit = glmer(as.formula(form), family=binomial, data=df[AKI_01_23==0])
fit = glmer(as.formula(form), family=binomial, data=df[AKI_01_23==1])
makeSupTbl4(fit)

form = paste0('inhos_mortality ~ ',paste0(setdiff(cox_risk_factors,'AKI_stage'), collapse = "+"),"+(1|Center2)")
fit = glmer(as.formula(form), family=binomial, data=df[AKI_YN==0])
vif(fit)
cox.zph(fit)

# Non AKI
fit = coxph(as.formula(paste0('Surv(inhos_duration, inhos_mortality==1) ~ ',
                              paste0(setdiff(cox_risk_factors,'AKI_stage'), collapse = "+"),"+cluster(Center2)")), data=df[AKI_YN==0])
vif(fit)
b = makeSupTbl4(fit)
tbl_merge(list(a,b))
cox.zph(fit)

# Supple table 1 ----------------------------------------------------------

# missing values frequencies
library(gt)

table_vars2 = c(table_vars, c("Input_ICUD2","Input_ICUD3","Input_ICUD7","Input_ICUDL",
                             "Output_ICUD2","Output_ICUD3","Output_ICUD7","Output_ICUDL"
                             ))
missing_count_prob_tbl <- df[,..table_vars2][,lapply(.SD,\(x) list(count = sum(is.na(x)), prob = sum(is.na(x)/nrow(df))*100 )) |> 
  rbindlist(idcol="Variable")];
missing_count_prob_tbl[order(-count)] |> gt()

missing_count_prob_tbl[order(-count)] |> fwrite("missing_value_count_prob.csv")
  ggplot(aes(x=reorder(Variable,-prob), y=prob))+
  geom_col() +
  theme(axis.text.x = element_text(angle=90))

# Supple table 2 ----------------------------------------------------------

makeBaselineTable(data = df,
                  include = table_vars,
                  dicho_vars = c('s_pulmonary','s_abdominal','s_urinary',
                                 'bac_blood_sepsis_gp','bac_blood_sepsis_gn',
                                 's_skinsoft','s_other','s_unclear','Septic_shock_ICUD1'),
                  by = "AKI_stage")


# Supple table 3 ----------------------------------------------------------

# SA-AKI: Sepsis associated AKI
# Sepsis이면서 동시에 AKI가 있는 경우
# AKI 진단이 Sepsis 진단(TZ)으로부터 7일 이내 이루어진 경우

temp = df[AKI_YN==1, .(SubjectNo, AKI_initial, AKI_Day1, AKI_Day2, AKI_Day3, AKI_Day7, CRRT_within_7days)]
temp[AKI_initial == 0 & 
       AKI_Day1 == 0 & AKI_Day2 == 0 &
       AKI_Day3 == 0 & is.na(AKI_Day7), AKI_Day7 := 3]
temp_melt = melt(temp, id.vars = 'SubjectNo')[order(SubjectNo)]
any_sa_aki_df = temp_melt[temp_melt[,.I[value>0], by=SubjectNo]$V1]
any_sa_aki_df[df, on=.(SubjectNo), AKI_stage := i.AKI_stage]

max_sa_aki_df = temp_melt[temp_melt[,.I[which.max(value)], SubjectNo]$V1]
max_sa_aki_df[df, on=.(SubjectNo), AKI_stage := i.AKI_stage]


for(i in c('','Stage 1','Stage 2','Stage 3')){
  if (i==''){
    cat('AKI stage: All','\n')
    any_sa_aki_df[!is.na(SubjectNo),first(variable), SubjectNo][,.(SubjectNo, day = fcase(V1=='AKI_initial',0,
                                                                                                           V1=='AKI_Day1',1,
                                                                                                           V1=='AKI_Day2',2,
                                                                                                           V1=='AKI_Day3',3,
                                                                                                           V1=='AKI_Day7',7,
                                                                                                           default=NA))][,summary(day)] |> print()
    any_sa_aki_df[!is.na(SubjectNo),first(variable), SubjectNo][,.(SubjectNo, day = fcase(V1=='AKI_initial',0,
                                                                                                           V1=='AKI_Day1',1,
                                                                                                           V1=='AKI_Day2',2,
                                                                                                           V1=='AKI_Day3',3,
                                                                                                           V1=='AKI_Day7',7,
                                                                                                           default=NA))][,.N,day][,prop:=N/sum(N)*100][order(day)] |> print()
    max_sa_aki_df[,.(day = fcase(variable=='AKI_initial',0,
                                                       variable=='AKI_Day1',1,
                                                       variable=='AKI_Day2',2,
                                                       variable=='AKI_Day3',3,
                                                       variable=='AKI_Day7',7,
                                                       default=NA))][,summary(day)] |> print()
    max_sa_aki_df[,.(day = fcase(variable=='AKI_initial',0,
                                                       variable=='AKI_Day1',1,
                                                       variable=='AKI_Day2',2,
                                                       variable=='AKI_Day3',3,
                                                       variable=='AKI_Day7',7,
                                                       default=NA))][,.N,day][,prop:=N/sum(N)*100][order(day)] |> print()
    
  } else {
    cat('\nAKI stage: ',i,'\n\n')
    any_sa_aki_df[!is.na(SubjectNo) & AKI_stage == i,first(variable), SubjectNo][,.(SubjectNo, day = fcase(V1=='AKI_initial',0,
                                                                                                           V1=='AKI_Day1',1,
                                                                                                           V1=='AKI_Day2',2,
                                                                                                           V1=='AKI_Day3',3,
                                                                                                           V1=='AKI_Day7',7,
                                                                                                           default=NA))][,summary(day)] |> print()
    any_sa_aki_df[!is.na(SubjectNo) & AKI_stage == i,first(variable), SubjectNo][,.(SubjectNo, day = fcase(V1=='AKI_initial',0,
                                                                                                           V1=='AKI_Day1',1,
                                                                                                           V1=='AKI_Day2',2,
                                                                                                           V1=='AKI_Day3',3,
                                                                                                           V1=='AKI_Day7',7,
                                                                                                           default=NA))][,.N,day][,prop:=N/sum(N)*100][order(day)] |> print()
    max_sa_aki_df[AKI_stage == i,.(day = fcase(variable=='AKI_initial',0,
                                                       variable=='AKI_Day1',1,
                                                       variable=='AKI_Day2',2,
                                                       variable=='AKI_Day3',3,
                                                       variable=='AKI_Day7',7,
                                                       default=NA))][,summary(day)] |> print()
    max_sa_aki_df[AKI_stage == i,.(day = fcase(variable=='AKI_initial',0,
                                 variable=='AKI_Day1',1,
                                 variable=='AKI_Day2',2,
                                 variable=='AKI_Day3',3,
                                 variable=='AKI_Day7',7,
                                 default=NA))][,.N,day][,prop:=N/sum(N)*100][order(day)] |> print()
  }
  
}

df[!is.na(Cr_ICUD7),.N,.(AKI_stage, AKI_status_day7)]

mytable(AKI_stage ~ AKI_status_day7, df)

# Supple table 4  -----------------------------------------------------------------

vars = df[,.(inhos_mortality, icu_mortality, 
                  icu_duration, inhos_duration,
                  TZ_to_disch_duration,
                  CRRT, RRT_discharge_HD_PD,
                  ICUStayInvasive_OU, invasive_MV_duration,
                  AntbBEFCurrent, AppInitEmpThe_re, gp_ag_cl, AdjCSTx,
                  FirstSourceNonSurg, FirstSourceControl,
                  Vasopressors_ICUD1, Inotropes_ICUD1,
                  Vasopressor_support, Inotrope_support,
                  #Bundlie compliance
                  lac1h, bd1h, anti1h, br1h, app1h,
                  Input_ICUD1, Output_ICUD1,
                  input_output_ICUD1, input_output_ICUD2,
                  input_output_ICUD3, input_output_ICUD7
)] |> names()
df[,.N,by=is.na(input_output_ICUD7)]
tbl = tbl_summary(
  data=df,
  by="AKI_stage",
  include = vars,
  missing = 'no',
  statistic = list(
    all_continuous() ~ "{mean} ± {sd}",
    all_categorical() ~ "{n} ({p})",
    c(inhos_duration,icu_duration, TZ_to_disch_duration,
      invasive_MV_duration) ~ "{median} [{p25}-{p75}]"
  ),
  value = list(
    AntbBEFCurrent ~ 1,
    AppInitEmpThe_re ~ 1
    # SEX ~ 'Male'
  ),
  digits = list(
    all_continuous() ~ 1,
    all_categorical() ~ c(0,1)
  )
) |> 
  add_overall() |> 
  add_p(pvalue_fun = ~style_pvalue(., digits=3),
        test = list(
          all_continuous() ~ "aov"
        )) |> 
  modify_footnote(everything() ~ NA)

tbl_1 = tbl_summary(
  data = df[Hr1SSC_BRInf30min==1],
  include = TZ_to_bolus_dur,
  statistic = list(
    TZ_to_bolus_dur ~ "{median} [{p25}-{p75}]"
  ),
  by = "AKI_stage",
  missing='no',
  digits = list(all_continuous() ~ 1,
                all_categorical() ~ c(0,1)) 
) |> add_overall() |> 
  add_p(pvalue_fun = ~style_pvalue(., digits=3))

tbl_2 = tbl_summary(
  data =df[Hr1SSC_ApplyVP==1],
  include = c("fluid_to_vaso_dur","TZ_to_vaso_dur","vaso_prior_fluid_yn"),
  statistic = list(
    c("fluid_to_vaso_dur","TZ_to_vaso_dur") ~ "{median} [{p25}-{p75}]",
    all_categorical() ~ "{n} ({p})"
  ),
  by = "AKI_stage",
  missing='no',
  digits = list(all_continuous() ~ 1,
                all_categorical() ~ c(0,1)),
  value = vaso_prior_fluid_yn ~ 1
) |> add_overall() |> 
  add_p(pvalue_fun = ~style_pvalue(., digits=3))

tbl_3 = tbl_summary(
  data=df,
  include = LimitLifeSustain,
  by='AKI_stage',
  missing='no',
  digits = list(all_continuous() ~ 1,
                all_categorical() ~ c(0,1)),
  value = LimitLifeSustain ~ 1
)|> add_overall() |> 
  add_p(pvalue_fun = ~style_pvalue(., digits=3))

tbl_stack(
  tbls = list(
    tbl, tbl_1, tbl_2, tbl_3
  )
)

# Supple table 5 ----------------------------------------------------------
require(gt)
makeSupTbl5(75, target = 'br1h', variables = cox_risk_factors)
makeSupTbl5(70, target = 'br1h', variables = cox_risk_factors)
makeSupTbl5(65, target = 'br1h', variables = cox_risk_factors)
makeSupTbl5('both75', target = 'br1h', variables = cox_risk_factors)
makeSupTbl5(60, var = 'app1h')
makeSupTbl5(65, var = 'app1h')
makeSupTbl5(70, var = 'app1h')
makeSupTbl5(75, var = 'app1h')
makeSupTbl5('baseline', var = 'br1h')
makeSupTbl5('baseline', var = 'app1h')
makeSupTbl5('both75', var = 'app1h')

# Supple table 6 ----------------------------------------------------------

df[br1h==0 & SubjectNo %notin% no_sa_aki_id,.N,AKI_stage][,prop:=N/sum(N)*100][]
df[Hr1SSC_BRInf30min==0 & SubjectNo %notin% no_sa_aki_id,.N,AKI_stage][,prop:=N/sum(N)*100][]
mytable(~AKI_stage, df |> dplyr::filter(br1h==0))
df[Hr1SSC_BRInf30min==0 & SubjectNo %notin% no_sa_aki_id] |> 
     mytable(~AKI_stage)
suppleTbl6(data=df[br1h==0])
suppleTbl6(data=df[Hr1SSC_BRInf30min==0])
require(epiDisplay)
df[Hr1SSC_BRInf30min==0 & SubjectNo %notin% no_sa_aki_id,AKI_stage] |> 
  tab1(decimal=1, cum.percent=F)


# Supple table 7 ----------------------------------------------------------

df[,.N,Center]
center = unique(df$Center2)

suppleTbl7(center[2])

tbl_list =lapply(center, suppleTbl7)
tbl_merge(tbl_list, tab_spanner = center)




