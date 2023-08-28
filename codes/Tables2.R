# Table 1 ------------------------------------------------
table_vars = df[,.(Age,SEX, BMI,
                   Comorbidity_MOSIAC_DM,Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
                   Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,Charlson_comorbidity_index_total, 
                   sofa_initial, SOFA_Renal, SOFA_except_renal,  SAPS3_ICUD1,
                   TZ_ICUAdm_gap_hour,
                   IVS_SBP,IVS_DBP, IVS_MBP,IVS_HR,IVS_BT, 
                   Lactate, Hb, Plt, BUN, Cr, Bilirubin, AST, ALT, Albumin, INR, 
                   CRP_imp, Procalcitonin_imp, ArterialPH, PaCO2, PaO2, Bicarbonate,
                   s_pulmonary, s_abdominal, s_urinary, s_skinsoft, s_other, s_unclear,
                   TypeOfInfection, 
                   bac_blood_sepsis, bac_blood_sepsis_gp, bac_blood_sepsis_gn, fung_blood_sepsis,
                   AntbBEFCurrent, AdjCSTx, FirstSourceControl, FirstSourceNonSurg, SurgSourceCont, 
                   AppInitEmpThe_re, glycopeptide, aminoglycoside, colistin, gp_ag_cl, 
                   GramPosYN, GramNegYN,
                   BacMTYN, Septic_shock, 
                   Vasopressors_ICUD1, Inotropes_ICUD1,
                   Vasopressor_support, Inotrope_support,
                   lac1h, bd1h, anti1h, br1h, app1h,
                   lac3h, bd3h, anti3h, br3h, app3h,
                   lac6h, bd6h, anti6h, br6h, app6h,
                   MV, Input_ICUD1, Output_ICUD1,
                   CRRT, RRT_discharge,
                   inhos_mortality, icu_mortality, ICU_LOS, Hospital_LOS)] |> names()

dichotomous_vars = df[,.(SEX, Comorbidity_MOSIAC_DM,Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
                         Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,
                         s_pulmonary, s_abdominal,
                         s_urinary, s_skinsoft, s_other, s_unclear, 
                         GramPosYN, GramNegYN, fung_blood_sepsis, AppInitEmpThe_re,
                         BacMTYN, Septic_shock,
                         Vasopressor_support, Inotrope_support,
                         lac1h, bd1h, anti1h, br1h, app1h,
                         lac3h, bd3h, anti3h, br3h, app3h,
                         lac6h, bd6h, anti6h, br6h, app6h,
                         CRRT,MV)] |> names()

makeBaselineTable(data = df,
                  include = table_vars,
                  dicho_vars = dichotomous_vars,
                  by = "AKI_YN")

# Table 2 -----------------------------------------------------------------
df[,vaso_inotrope_ICUD1 := pmax(Vasopressor_support, Inotrope_support)]
cox_risk_factors = df[,.( Age,SEX, BMI, AKI_stage,
                          Charlson_comorbidity_index_total, 
                          SOFA_ICUD1, SAPS3_ICUD1, 
                          Lactate, Hb, Plt, Bilirubin, Albumin, CRP_imp, ArterialPH, PaCO2, PaO2,
                          s_pulmonary, s_abdominal, s_urinary, s_skinsoft, s_other, s_unclear,
                          TypeOfInfection, 
                          bac_blood_sepsis_gp, bac_blood_sepsis_gn, fung_blood_sepsis,
                          AntbBEFCurrent, AdjCSTx,
                          AppInitEmpThe_re,  GramPosYN, GramNegYN, gp_ag_cl, 
                          lac1h, bd1h, anti1h, br1h, app1h,
                          MV, input_output_ICUD1, CRRT, vaso_inotrope_ICUD1 )] |> names()
cox_dichotomous_vars = df[,.(SEX, 
                             # Comorbidity_MOSIAC_DM,Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
                             # Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,
                             s_pulmonary, s_abdominal, s_urinary, s_skinsoft, s_other, s_unclear, 
                             GramPosYN, GramNegYN, fung_blood_sepsis, AppInitEmpThe_re,
                             gp_ag_cl,
                             lac1h, bd1h, anti1h, br1h, app1h,
                             CRRT, MV, vaso_inotrope_ICUD1)] |> names()

library(gtsummary)
tbl2 = tbl_summary(
  data=df,
  by="inhos_mortality",
  include = cox_risk_factors,
  missing = 'no',
  statistic = list(
    all_continuous() ~ "{mean} ± {sd}",
    all_categorical() ~ "{n} ({p})"
  ),
  type = list(
    cox_dichotomous_vars ~ 'dichotomous'
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
  add_p(pvalue_fun = ~style_pvalue(., digits=3),
        test = list(
          all_continuous() ~ "aov"
        )) |> 
  modify_footnote(everything() ~ NA) |>
  modify_table_body(~.x |> relocate(stat_2, .after=stat_0))
tbl2
tbl2$table_body$p.value
multi_cox_target = tbl2$table_body |> 
  dplyr::filter(p.value<0.2) |> 
  select(variable) |> c()

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



# Supple table 1 ----------------------------------------------------------

# missing values frequencies
missing_count_prob_tbl <- df[,..table_vars][,lapply(.SD,\(x) list(count = sum(is.na(x)), prob = sum(is.na(x)/nrow(df))*100 )) |> 
                                              rbindlist(idcol="Variable")];missing_count_prob_tbl
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
                  dicho_vars = dichotomous_vars,
                  by = "AKI_stage")


# Supple table 3 ----------------------------------------------------------

makeBaselineTable(data = df,
                  include = table_vars,
                  dicho_vars = dichotomous_vars,
                  by = "AKI_stage_lowest")


# Supple table 4 ----------------------------------------------------------

vars = df[,.(Age,SEX, BMI,
             Comorbidity_MOSIAC_DM,Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
             Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,Charlson_comorbidity_index_total, 
             sofa_initial, SAPS3_ICUD1,
             Lactate, Hb, Plt, Bilirubin, Albumin, CRP_imp, ArterialPH, PaCO2, PaO2,
             s_pulmonary, s_abdominal, s_urinary, s_skinsoft, s_other, s_unclear,
             TypeOfInfection, 
             bac_blood_sepsis, bac_blood_sepsis_gp, bac_blood_sepsis_gn, fung_blood_sepsis,
             AntbBEFCurrent, AdjCSTx,AppInitEmpThe_re, GramPosYN, GramNegYN, gp_ag_cl, 
             lac1h, bd1h, anti1h, br1h, app1h,
             MV, input_output_ICUD1, CRRT, vaso_inotrope_ICUD1)] |> names()

dicho_vars = df[,.(SEX,Comorbidity_MOSIAC_DM, Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung, Comorbidity_MOSIAC_CKD,
                   Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,
                   s_pulmonary, s_abdominal, s_urinary, s_skinsoft, s_other, s_unclear,
                   bac_blood_sepsis, bac_blood_sepsis_gp, bac_blood_sepsis_gn, fung_blood_sepsis,
                   AntbBEFCurrent, AdjCSTx, AppInitEmpThe_re, GramPosYN, GramNegYN, gp_ag_cl,
                   lac1h, bd1h, anti1h, br1h, app1h,
                   MV,CRRT, vaso_inotrope_ICUD1)] |> names()

tbl3 = tbl_summary(
  data=df,
  by="stage3_yn",
  include = vars,
  missing = 'no',
  statistic = list(
    all_continuous() ~ "{mean} ± {sd}",
    all_categorical() ~ "{n} ({p})"
  ),
  type = list(
    dicho_vars ~ 'dichotomous'
  ),
  value = list(
    SEX ~ 'Male',
    TypeOfInfection ~ "nosocomial",
    AntbBEFCurrent ~ 1
  ),
  digits = list(
    all_continuous() ~ 1,
    all_categorical() ~ c(0,1)
  )
) |> 
  add_p(pvalue_fun = ~style_pvalue(., digits=3),
        test = list(
          all_continuous() ~ "aov"
        )) |> 
  modify_footnote(everything() ~ NA)

multi_lr_target = tbl3$table_body |> 
  dplyr::filter(p.value<0.2) |>
  dplyr::filter(variable != "CRRT") |> 
  select(variable) |> c()

form = paste0("stage3_yn~",paste0(multi_lr_target$variable,collapse = "+"))
lr_fit = glm(as.formula(form),family=binomial, data=df)
tbl_regression(
  x = lr_fit,
  exponentiate = T,
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


