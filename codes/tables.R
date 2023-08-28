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

# CRRT, RRT_discharge,

dichotomous_vars = df[,.(SEX, Comorbidity_MOSIAC_DM,Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
                         Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,
                         s_pulmonary, s_abdominal,
                         s_urinary, s_skinsoft, s_other, s_unclear, 
                         bac_blood_sepsis_gp, bac_blood_sepsis_gn, 
                         Septic_shock_ICUD1,  ICUStayInvasive_ICUD1,
                         ICUStayInvasive_OU,ICUStayECMO_ICUD1,ICUStayECMO_OU)] |> names()

setdiff(dichotomous_vars, table_vars)
makeBaselineTable(data = df,
                  include = table_vars,
                  dicho_vars = dichotomous_vars,
                  by = "AKI_YN")

df[AKI_YN==1,.N,AKI_stage][,.(AKI_stage, N, prop=N/sum(N)*100)]



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


# Table 2 -----------------------------------------------------------------

tbl2_vars = df[,.(Age,SEX,BMI,
             Center2, CFScore,
             # Comorbidities
             Comorbidity_MOSIAC_DM, Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
             Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,Charlson_comorbidity_index_total, 
             
             sofa_initial, SAPS3_ICUD1,
             # Source of infection
             s_pulmonary, s_abdominal, s_urinary, s_skinsoft, s_other, s_unclear,
             TypeOfInfection,
             bac_blood_sepsis_gp, bac_blood_sepsis_gn,
             fung_blood_sepsis,
             AntbBEFCurrent, AppInitEmpThe_re
             )] |> names()

dichotomous_vars = df[,.(SEX, Comorbidity_MOSIAC_DM,Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
                         Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,
                         s_pulmonary, s_abdominal,
                         fung_blood_sepsis,AppInitEmpThe_re,
                         s_urinary, s_skinsoft, s_other, s_unclear, 
                         bac_blood_sepsis_gp,bac_blood_sepsis_gn)] |> names()

makeTable3(data = df,
                  include = tbl3_vars,
                  dicho_vars = dichotomous_vars,
                  by = "AKI_01_23")

vars = df[,.(Age,SEX, BMI, CFScore,
             Comorbidity_MOSIAC_DM, Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
             Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM, 
             # Charlson_comorbidity_index_total,
             sofa_initial,
             # Lactate, Hb, Plt, Bilirubin, Albumin, CRP_imp, ArterialPH, PaCO2, PaO2,
             s_pulmonary, s_abdominal, s_urinary, 
             # s_skinsoft, s_other, s_unclear,
             TypeOfInfection,
             bac_blood_sepsis_gp, bac_blood_sepsis_gn,
             AntbBEFCurrent2)] |> names()

dicho_vars = df[,.(SEX,Comorbidity_MOSIAC_DM, Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung, Comorbidity_MOSIAC_CKD,
                   Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,
                   # s_pulmonary, s_abdominal, s_urinary, s_skinsoft, s_other, s_unclear,
                   bac_blood_sepsis_gp, bac_blood_sepsis_gn,
                   AntbBEFCurrent2)] |> names()

library(lme4)
library(car)
form = paste0("AKI_01_23~",paste0(vars,collapse = "+"),'+(1|Center2)')
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
df[,shock_index := HR/SBP]
cox_risk_factors = df[,.( Age,SEX, BMI, CFScore,
                          AKI_stage, 
                          Comorbidity_MOSIAC_DM, Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
                          Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM, 
                          SOFA_ICUD1, SAPS3_ICUD1, 
                          # Lactate, Hb, Plt, Bilirubin, Albumin, CRP_imp, ArterialPH,
                          s_pulmonary, s_abdominal, s_urinary,
                          # s_skinsoft, s_other, s_unclear,
                          TypeOfInfection, 
                          bac_blood_sepsis_gp, bac_blood_sepsis_gn,
                          AntbBEFCurrent2,
                          AppInitEmpThe_re, 
                          lac1h, bd1h, anti1h, br1h, app1h)] |> names()
form = paste0("Surv(inhos_duration, inhos_mortality==1)~", paste0(cox_risk_factors, collapse = "+"),"+cluster(Center2)")
fit = coxph(as.formula(form), data=df[AKI_YN==1])
mytable(AKI_01_23 ~ shock_index, df)
fit
#AKI
a=makeSupTbl4(fit)
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
b
tbl_merge(list(a,b))
cox.zph(fit)
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



# Supple table 1 ----------------------------------------------------------

# missing values frequencies
library(gt)
missing_count_prob_tbl <- df[,..table_vars][,lapply(.SD,\(x) list(count = sum(is.na(x)), prob = sum(is.na(x)/nrow(df))*100 )) |> 
  rbindlist(idcol="Variable")];
missing_count_prob_tbl[order(-count)] |> gt()


missing_count_prob_tbl[order(-count)] |> fwrite("missing_value_count_prob.csv")
  ggplot(aes(x=reorder(Variable,-prob), y=prob))+
  geom_col() +
  theme(axis.text.x = element_text(angle=90))

df[,.N,AppInitEmpThe]
df[]


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

tbl2_vars = df[,.(inhos_mortality, icu_mortality, 
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
library(gtsummary)

tbl2 = tbl_summary(
  data=df,
  by="AKI_stage",
  include = tbl2_vars,
  missing = 'no',
  statistic = list(
    all_continuous() ~ "{mean} ± {sd}",
    all_categorical() ~ "{n} ({p})",
    c(inhos_duration,icu_duration, TZ_to_disch_duration,
      invasive_MV_duration) ~ "{median} [{p25}-{p75}]"
  ),
  # type = list(
  #   # cox_dichotomous_vars ~ 'dichotomous'
  # ),
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

tbl2_1 = tbl_summary(
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

tbl2_2 = tbl_summary(
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

tbl2_3 = tbl_summary(
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
    tbl2, tbl2_1, tbl2_2, tbl2_3
  )
)

df[Hr1SSC_BRInf30min==1, .N, AKI_stage]
df[Hr1SSC_ApplyVP==1, .N, AKI_stage]



# Supple table 5 ----------------------------------------------------------


require(gt)
makeTable5 = function(EGFR, var){
  temp = set_cohort(baseline='Cr_baseline_MDRD', eGFR = EGFR)
  form = paste0("inhos_mortality~", paste0(cox_risk_factors, collapse = "+"),'+(1|Center2)')
  overall = glmer(as.formula(form), family=binomial, data=temp) 
  aki = glmer(as.formula(form), family=binomial, data=temp[AKI_YN==1]) 
  noaki = glmer(as.formula(paste0('inhos_mortality ~ ',paste0(setdiff(cox_risk_factors,'AKI_stage'), collapse = "+"),"+(1|Center2)")), family=binomial, data=temp[AKI_YN==0])
  severAKI = glmer(as.formula(paste0('inhos_mortality ~ ',paste0(setdiff(cox_risk_factors,'AKI_stage'), collapse = "+"),"+(1|Center2)")), family=binomial, data=temp[AKI_01_23==1])
  no_one_aki = glmer(as.formula(paste0('inhos_mortality ~ ',paste0(setdiff(cox_risk_factors,'AKI_stage'), collapse = "+"),"+(1|Center2)")), family=binomial, data=temp[AKI_01_23==0])
  
  tbl1 = tbl_regression(overall, exponentiate = T)$table_body |> dplyr::filter(variable==var) |> dplyr::select(estimate, conf.low, conf.high)
  tbl2 = tbl_regression(aki, exponentiate = T)$table_body |> dplyr::filter(variable==var) |> dplyr::select(estimate, conf.low, conf.high)
  tbl3 = tbl_regression(noaki, exponentiate = T)$table_body |> dplyr::filter(variable==var) |> dplyr::select(estimate, conf.low, conf.high)
  tbl4 = tbl_regression(severAKI, exponentiate = T)$table_body |> dplyr::filter(variable==var) |> dplyr::select(estimate, conf.low, conf.high)
  tbl5 = tbl_regression(no_one_aki, exponentiate = T)$table_body |> dplyr::filter(variable==var) |> dplyr::select(estimate, conf.low, conf.high)
  result = rbind(overall = tbl1,
                 AKI = tbl2,
                 NoAKI = tbl3,
                 severeAKI = tbl4,
                 nooneAKI = tbl5)
  result |> 
    mutate("OR (95% CI)" = paste0(
      format(round(estimate,2),nsmall=2), ' (', 
      format(round(conf.low,2),nsmall=2),'-',
      format(round(conf.high,2), nsmall=2),')' 
    )) |> gt(rownames_to_stub = T)
}

makeTable5(60, var = 'br1h')
makeTable5(60, var = 'app1h')
makeTable5(65, var = 'br1h')
makeTable5(65, var = 'app1h')
makeTable5(70, var = 'br1h')
makeTable5(70, var = 'app1h')
makeTable5(75, var = 'br1h')
makeTable5(75, var = 'app1h')
makeTable5('baseline', var = 'br1h')
makeTable5('baseline', var = 'app1h')
makeTable5('both75', var = 'br1h')
makeTable5('both75', var = 'app1h')

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
df[,.N,AKI_Day2]
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

