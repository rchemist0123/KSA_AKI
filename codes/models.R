# 2023-05-19 --------------------
# binomial LR
# 1) AKI vs non-AKI
lr_risk_factors
library(survival)
library(gt)

aki_yn_uni_lr = uni_reg(data=df, vars = lr_risk_factors,
  y="AKI_YN", type="lr", digit=2)

aki_yn_uni_lr |>  
  gt() |> 
  tab_style(
    style = list(
      cell_fill(color="yellow"),
      cell_text(weight = "bold")
    ),
    location = cells_body(
      columns = p.value,
      rows = p.value < 0.05
    )
  )

multi_vars = aki_yn_uni_lr[p.value =="<0.001" | as.numeric(p.value)<0.2, variable]
form = paste0("AKI_YN~", paste0(multi_vars, collapse="+"))
aki_yn_mult_lr_fit = glm(as.formula(form), family=binomial, data=df)
aki_yn_mult_lr$variable
merge(aki_yn_uni_lr, aki_yn_mult_lr, by="variable", all.x=T)

aki_yn_mult_lr = multi_reg(data=df, vars=multi_vars,
          y = "AKI_YN", type = "lr", digit = 2)
require(broom)
library(dplyr)
tidy(aki_yn_mult_lr_fit,  exponentiate = T, conf.int=T)  |> 
  mutate(OR_CI = paste0(
    format(round(estimate,2),2), " (",
    format(round(conf.low,2),2), "—",
    format(round(conf.high,2),2), ")"
      ),
    p.value = ifelse(p.value<0.001,"<0.001",round(p.value,3)))|> 
  select(term, OR_CI, p.value)  |> 
  gt() |> 
  tab_style(
    style = list(
      cell_fill(color="yellow"),
      cell_text(weight = "bold")
    ),
    location = cells_body(
      columns = p.value,
      rows = p.value < 0.05
    )
  )

# 2) Stage 3 vs Non-stage 3

lr_risk_factors
library(survival)
library(gt)
aki_yn_uni_lr = uni_reg(data=df, vars = lr_risk_factors,
  y="stage3_yn_lowest", type="lr", digit=2)
aki_yn_uni_lr |>  
  gt() |> 
  tab_style(
    style = list(
      cell_fill(color="yellow"),
      cell_text(weight = "bold")
    ),
    location = cells_body(
      columns = p.value,
      rows = p.value < 0.05
    )
  )
multi_vars = aki_yn_uni_lr[p.value =="<0.001" | as.numeric(p.value)<0.2, variable]

form = paste0("stage3_yn_lowest~", paste0(multi_vars,collapse="+"))
aki_yn_mult_lr_fit = glm(as.formula(form), family=binomial, data=df)

tidy(aki_yn_mult_lr_fit,  exponentiate = T, conf.int=T)  |> 
  mutate(OR_CI = paste0(
    format(round(estimate,2),2), " (",
    format(round(conf.low,2),2), "—",
    format(round(conf.high,2),2), ")"
      ),
    p.value = ifelse(p.value<0.001,"<0.001",round(p.value,3)))|> 
  select(term, OR_CI, p.value)  |> 
  gt() |> 
  tab_style(
    style = list(
      cell_fill(color="yellow"),
      cell_text(weight = "bold")
    ),
    location = cells_body(
      columns = p.value,
      rows = p.value < 0.05
    )
  )

# 2023-05-11 --------------------------------
# multinomial logistic regression with nnet
require(nnet)
form
nnet_fit = multinom(as.formula(form), data=df)
exp(coef(nnet_fit))
exp(confint(nnet_fit))
library(broom)
tidy(nnet_fit, exponentiate = T, conf.int = T) |>  
  dplyr::filter(y.level =="Stage 3")  |> 
  mutate(or_ci = paste0(format(round(estimate,2),2), " (", 
                        format(round(conf.low,2),2),"—", 
                        format(round(conf.high,2),2),")"),
        pval = ifelse(p.value < 0.001,"<0.001",round(p.value,4))) |> 
  select(term, or_ci, pval) |> 
  gt() |> 
  
library(stargazer)
stargazer(nnet_fit,  type="text", summary=T)

# 2023-05-10 ----------------------------------
# multinomial non, 1&2, 3
# uni:MDRD

lr_risk_factors <- df[,.(Age,SEX, BMI,
                         Comorbidity_MOSIAC_DM,Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
                         Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM,Charlson_comorbidity_index_total,
                         SOFA_except_renal, SAPS3_ICUD1, IVS_MBP,IVS_HR,IVS_BT, Lactate, Hb, Plt, Bilirubin, 
                         Albumin, CRP_imp, ArterialPH, PaCO2, PaO2,
                         s_pulmonary, s_abdominal, s_urinary, s_skinsoft, s_other, s_unclear,
                         TypeOfInfection, 
                         bac_blood_sepsis_gp, bac_blood_sepsis_gn, fung_blood_sepsis,
                         AntbBEFCurrent, AdjCSTx, FirstSourceNonSurg, FirstSourceControl, 
                         AppInitEmpThe_re, glycopeptide, aminoglycoside, colistin, gp_ag_cl, Septic_shock, 
                         lac1h, bd1h, anti1h, br1h, app1h,
                         MV,input_output_ICUD1)] |> names()

require(VGAM)
aki_2_lr_uv = uni_reg(data=df, vars = lr_risk_factors, y="AKI_stage_2", type="multi", digit=3,)
aki_2_lr_uv[rn %like% ":1"] |> gt()
aki_2_lr_uv[rn %like% ":2"] |> gt()


# 2023-06-16 check VIF  & Multivariable
multi_vars = aki_2_lr_uv[rn %like% ":1" & pval < 0.2, covariate]
multi_vars2 = aki_2_lr_uv[rn %like% ":2" & pval < 0.2, covariate]

multi_target_vars = union(multi_vars, multi_vars2)

form = paste0("AKI_stage_2 ~ ", paste0(multi_target_vars,collapse = "+"))
# multi_fit = VGAM::vglm(as.formula(form), data=df, family=multinomial(refLevel="Non-AKI"))

require(nnet)
multi_fit = multinom(as.formula(form), data = df)

tidy(multi_fit, exponentiate = T) |> View()
df[,AKI_stage_2 := factor(AKI_stage_2, levels=c("Non-AKI","Stage 12","Stage 3"))]
car::vif(glm(as.formula(form), data=df, family=binomial)) |> as.data.table(keep.rownames = T) |> gt()
multi_reg(data=df, vars=multi_targeted_vars, y="AKI_stage_2", type="multi")

# 2023-04-21 -----------------------------------
# Multinomial logistic regression
# VGAM:: vglm object : S4 class -> @로 접근.
# install.packages("VGAM")
# require(VGAM)



form <- paste0("AKI_stage_2~", paste0(risk_factors, collapse="+"))
multi_lr_fit <- vglm(as.formula(form), multinomial(refLevel="Non-AKI"), data=df)


rst_tbl  <- as.data.table(
  cbind(
    OR = round(exp(coef(multi_lr_fit)),2),
    ci_low = round(exp(confint(multi_lr_fit))[,1],2),
    ci_high = round(exp(confint(multi_lr_fit))[,2],2),
    pval = ifelse(smry@coef3[,4] < 0.001, "0.001", round(smry@coef3[,4],4))
  ),
  keep.rowname=T
)

rst_tbl[,OR_CI := paste0(format(OR, digit=2, nsmall=2), ' (',format(ci_low, digit=2,nsmall=2),"—",format(ci_high, digit=2,nsmall=2),")" )]

cbind(rst_tbl[rn %like% ":1",.(rn, OR_CI,p1=pval)],
      rst_tbl[rn %like% ":2",.(oc2 = OR_CI, p2=pval)],
      rst_tbl[rn %like% ":3",.(oc3 = OR_CI, p3=pval)]
) |> gt()
rst_tbl$rn
multi_lr_fit <- multinom(as.formula(form), data=df)
multi_lr_fit
smry <- summary(multi_lr_fit)
smry$coefficients

# Cox for in-hospital mortality

require(survival)
cox_risk_factors <- df[,.(AKI_stage_lowest, Age,SEX, BMI, 
      Charlson_comorbidity_index_total, 
      SAPS3_ICUD1, IVS_MBP, IVS_HR,IVS_BT, Lactate, Hb, Plt, Bilirubin, 
      Albumin, CRP_imp, ArterialPH, PaCO2, PaO2,
      s_pulmonary, s_abdominal, s_urinary, s_skinsoft, s_other, s_unclear,
      TypeOfInfection, 
      bac_blood_sepsis_gp, bac_blood_sepsis_gn, fung_blood_sepsis,
      AntbBEFCurrent, AdjCSTx, FirstSourceNonSurg, FirstSourceControl, 
      AppInitEmpThe_re, glycopeptide, aminoglycoside, colistin, gp_ag_cl, Septic_shock, 
      lac1h, bd1h, anti1h, br1h, app1h,
      MV,input_output_ICUD1, CRRT )] |> names()
df[,.N,AKI_stage_lowest]
df[,.N,AKI_stage]
form <- paste0("Surv(Hospital_LOS, inhos_mortality==1) ~", paste0(cox_risk_factors, collapse="+"))
cox_fit <- coxph(as.formula(form),data=df)
tidy(cox_fit, exponentiate = T, conf.int = T) |> 
  mutate(HR_CI = paste0(
    format(round(estimate,2),2), " (",
    format(round(conf.low,2),2), "—",
    format(round(conf.high,2),2), ")"
  ),
  p.value = ifelse(p.value<0.001,"<0.001",round(p.value,3)))|> 
  select(term, HR_CI, p.value)  |> 
  gt() |> 
  tab_style(
    style = list(
      cell_fill(color="yellow"),
      cell_text(weight = "bold")
    ),
    location = cells_body(
      columns = p.value,
      rows = p.value < 0.05
    )
  )


tbl_regression(
  x= cox_fit,
  exponentiate = TRUE,
  pvalue_fun = ~style_pvalue(., digits=3),
  show_single_row = c(s_pulmonary, s_abdominal,  s_skinsoft,s_urinary, s_other, s_unclear,
                      fung_blood_sepsis, AppInitEmpThe_re)
) |> 
  bold_p() |> 
  modify_column_merge(pattern = "{estimate} ({conf.low}-{conf.high})",
                    rows = !is.na(estimate)) |> 
  modify_column_hide(columns=ci) |> 
#  add_significance_stars(hide_ci = T, hide_se = T) |> 
  modify_header(
    label = "**Variable**",
    estimate = "**HR (95% CI)**",
    p.value = "**P value**"
  ) |> 
  modify_footnote(everything() ~ NA) |> 
  as_gt() |> 
  tab_style(
    style = list(cell_fill(color="yellow")),
    locations = cells_body(columns = p.value,
                           rows = p.value < 0.05 )
  )




# logistic regression, risk factors, Y: AKI_YN ------------------------------------------------------------

risk_factors <- df[,.(Age,SEX, BMI,
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
                      TypeOfInfection,
                      PosBloodCulture, 
                      GramPosYN, GramNegYN,
                      ICUStayCCRT_OU, MV)] |> names()

library(nnet)
multinom(AKI_stage ~ Age + SEX + BMI + 
                    ICUStayCCRT_OU + MV, data=df)
vars <- c("AKI_YN", risk_factors)
single_row_vars <- c('SEX', 's_pulmonary', 's_abdominal', 's_urinary', 's_skinsoft',
                     's_other', 's_unclear', 'PosBloodCulture', 'GramPosYN',
                     'GramNegYN')
aki_uvr_lr_tbl <- tbl_uvregression(
  data = df[,..vars],
  method = glm,
  method.args = list(family=binomial),
  y = AKI_YN,
  exponentiate = T,
  hide_n = T,
  pvalue_fun = ~style_pvalue(., digits=3),
  show_single_row = c(SEX, s_pulmonary, s_abdominal, s_urinary, s_skinsoft,
                      s_other, s_unclear, PosBloodCulture, GramPosYN,
                      GramNegYN),
  estimate_fun = ~style_ratio(.,digits = 2)
) |> modify_table_styling(
  column = estimate,
  rows = !is.na(estimate),
  cols_merge_pattern = "{estimate} ({conf.low} - {conf.high})"
) |> 
  modify_header(estimate ~ "**OR (95% CI)**") |> 
  modify_column_hide(c(ci))

aki_uv_lr_rst <- aki_uvr_lr_tbl$table_body |> setDT()
aki_mv_lr_vars <- aki_uv_lr_rst[p.value < 0.2,variable] 
aki_lr_form <- paste0("AKI_YN","~ ", paste0(aki_mv_lr_vars, collapse = " + "))
single_rows <- intersect(single_row_vars, aki_mv_lr_vars)
aki_mult_lr_fit <- glm(as.formula(aki_lr_form), family=binomial, data=df)

aki_mult_lr_tbl <- tbl_regression(
  x = aki_mult_lr_fit,
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
  modify_header(estimate ~ "**OR (95% CI)**") |> 
  modify_column_hide(c(ci))

tbl_merge(list(aki_uvr_lr_tbl, aki_mult_lr_tbl),
          tab_spanner = c(
            "**Univariable**",
            "**Multivariable**"
          )) |> 
  # gtsummary::show_header_names()
  as_gt() |>
  gt::tab_style(
    style= list(
      cell_fill(color="yellow"),
      cell_text(color="red", weight = "bold")
    ),
    locations = cells_body(
      columns = p.value_2,
      rows = p.value_2 <0.05
    )
  ) 
  # tab_row_group(
  #   label = "Bacteria",
  #   rows = 43:44
  # ) |> 
  # tab_row_group(
  #   label = "Source of infection",
  #   rows = 33:38
  # ) |> 
  # tab_row_group(
  #   label = "Laboratory finding",
  #   rows = 18:32
  # ) |> 
  # tab_row_group(
  #   label = "Vital sign",
  #   rows = 13:17
  # ) |> 
  # tab_row_group(
  #   label="Comorbidity",
  #   rows = 4:12
  # )

  # risk factor 2, Y: inhos-mortality ------------------------------------------
## 1) logistic regression ------------------

risk_factors2 <- df[,.(Age,SEX, BMI,
                       Comorbidity_MOSIAC_DM, Comorbidity_MOSIAC_Cardio,
                       Comorbidity_MOSIAC_Lung, Comorbidity_MOSIAC_CKD,
                       Comorbidity_MOSIAC_SMT,Comorbidity_MOSIAC_HMM,
                       Comorbidity_MOSIAC_IMM,Charlson_comorbidity_index_total, sofa_initial,
                       IVS_SBP,IVS_DBP, IVS_MBP,IVS_HR,IVS_BT,
                       Lactate, Hb, Plt, BUN, Cr,
                       AST, ALT, Albumin, INR, CRP, Procalcitonin, ArterialPH,
                       PaCO2, PaO2, Bicarbonate,
                       s_pulmonary, s_abdominal, s_urinary, 
                       s_skinsoft, s_other, s_unclear,
                       TypeOfInfection,
                       PosBloodCulture,
                       GramPosYN, GramNegYN,
                       AKI_stage, ICUStayCCRT_OU, MV)] |> names()

require(gtsummary)
require(data.table)
require(gt)
if(!is.data.table(data)) setDT(data)
vars <- c("inhos_mortality", risk_factors2)
single_row_vars <- c('SEX', 's_pulmonary', 's_abdominal', 's_urinary', 's_skinsoft',
                     's_other', 's_unclear', 'PosBloodCulture', 'GramPosYN',
                     'GramNegYN')
inhos_uvr_lr_tbl <- tbl_uvregression(
  data = df[,..vars],
  method = glm,
  method.args = list(family=binomial),
  y = inhos_mortality,
  exponentiate = T,
  hide_n = T,
  pvalue_fun = ~style_pvalue(., digits=3),
  show_single_row = c(SEX, s_pulmonary, s_abdominal, s_urinary, s_skinsoft,
                      s_other, s_unclear, PosBloodCulture, GramPosYN,
                      GramNegYN),
  estimate_fun = ~style_ratio(.,digits = 2)
) |> modify_table_styling(
  column = estimate,
  rows = !is.na(estimate),
  cols_merge_pattern = "{estimate} ({conf.low} - {conf.high})"
) |> 
  modify_header(estimate ~ "**OR (95% CI)**") |> 
  modify_column_hide(c(ci))

inhos_uv_lr_rst <- inhos_uvr_lr_tbl$table_body |> setDT()
inhos_mv_lr_vars <- inhos_uv_lr_rst[p.value < 0.2,variable] 
inhos_lr_form <- paste0("inhos_mortality","~ ", paste0(inhos_mv_lr_vars, collapse = " + "))
single_rows <- intersect(single_row_vars, inhos_mv_lr_vars)
inhos_mult_lr_fit <- glm(as.formula(inhos_lr_form), family=binomial, data=df)

inhos_mult_lr_tbl <- tbl_regression(
  x = inhos_mult_lr_fit,
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
  modify_header(estimate ~ "**OR (95% CI)**") |> 
  modify_column_hide(c(ci))

tbl_merge(list(inhos_uvr_lr_tbl, inhos_mult_lr_tbl),
          tab_spanner = c(
            "**Univariable**",
            "**Multivariable**"
          )) |> 
  # gtsummary::show_header_names()
  as_gt() |>
  gt::tab_style(
    style= list(
      cell_fill(color="yellow"),
      cell_text(color="red", weight = "bold")
    ),
    locations = cells_body(
      columns = p.value_2,
      rows = p.value_2 <0.05
    )
  ) |> 
  tab_row_group(
    label = "Bacteria",
    rows = 43:44
  ) |> 
  tab_row_group(
    label = "Source of infection",
    rows = 33:38
  ) |> 
  tab_row_group(
    label = "Laboratory finding",
    rows = 18:32
  ) |> 
  tab_row_group(
    label = "Vital sign",
    rows = 13:17
  ) |> 
  tab_row_group(
    label="Comorbidity",
    rows = 4:12
  )


## 2) cox regression --------------------------------------------------------
vars <- c("inhos_mortality","inhos_duration", risk_factors2)
single_row_vars <- c('SEX', 's_pulmonary', 's_abdominal', 's_urinary', 's_skinsoft',
                     's_other', 's_unclear', 'PosBloodCulture', 'GramPosYN',
                     'GramNegYN')
inhos_uvr_cox_tbl <- tbl_uvregression(
  data = df[,..vars],
  method = coxph,
  # method.args = list(family=binomial),
  y = Surv(inhos_duration, inhos_mortality),
  exponentiate = T,
  hide_n = T,
  pvalue_fun = ~style_pvalue(., digits=3),
  show_single_row = c(SEX, s_pulmonary, s_abdominal, s_urinary, s_skinsoft,
                      s_other, s_unclear, PosBloodCulture, GramPosYN,
                      GramNegYN),
  estimate_fun = ~style_ratio(.,digits = 2)
) |> modify_table_styling(
  column = estimate,
  rows = !is.na(estimate),
  cols_merge_pattern = "{estimate} ({conf.low} - {conf.high})"
) |> 
  modify_header(estimate ~ "**HR (95% CI)**") |> 
  modify_column_hide(c(ci))

inhos_uv_cox_rst <- inhos_uvr_cox_tbl$table_body |> setDT()
inhos_mv_cox_vars <- inhos_uv_cox_rst[p.value < 0.2,variable] 
inhos_cox_form <- paste0("Surv(inhos_duration, inhos_mortality)","~ ", paste0(inhos_mv_cox_vars, collapse = " + "))
single_rows <- intersect(single_row_vars, inhos_mv_cox_vars)
inhos_cox_mult_fit <- coxph(as.formula(inhos_cox_form), data=df)

# require(gdata)
# mv("uvr_tbl","uvr_cox_tbl")
inhos_mult_cox_tbl <- tbl_regression(
  x = inhos_cox_mult_fit,
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

tbl_merge(list(inhos_uvr_cox_tbl, inhos_mult_cox_tbl),
          tab_spanner = c(
            "**Univariable**",
            "**Multivariable**"
          )) |> 
  # gtsummary::show_header_names()
  as_gt() |>
  gt::tab_style(
    style= list(
      cell_fill(color="yellow"),
      cell_text(color="red", weight = "bold")
    ),
    locations = cells_body(
      columns = p.value_2,
      rows = p.value_2 <0.05
    )
  ) |> 
  tab_row_group(
    label = "Bacteria",
    rows = 43:44
  ) |> 
  tab_row_group(
    label = "Source of infection",
    rows = 33:38
  ) |> 
  tab_row_group(
    label = "Laboratory finding",
    rows = 18:32
  ) |> 
  tab_row_group(
    label = "Vital sign",
    rows = 13:17
  ) |> 
  tab_row_group(
    label="Comorbidity",
    rows = 4:12
  )


## 3) only early vs late --------------------------------

risk_factors3 <- df[,.(Age,SEX, BMI,
                       Comorbidity_MOSIAC_DM, Comorbidity_MOSIAC_Cardio,
                       Comorbidity_MOSIAC_Lung, Comorbidity_MOSIAC_CKD,
                       Comorbidity_MOSIAC_SMT,Comorbidity_MOSIAC_HMM,
                       Comorbidity_MOSIAC_IMM,Charlson_comorbidity_index_total, sofa_initial,
                       IVS_SBP,IVS_DBP, IVS_MBP,IVS_HR,IVS_BT,
                       Lactate, Hb, Plt, BUN, Cr,
                       AST, ALT, Albumin, INR, CRP, Procalcitonin, ArterialPH,
                       PaCO2, PaO2, Bicarbonate,
                       s_pulmonary, s_abdominal, s_urinary, 
                       s_skinsoft, s_other, s_unclear,
                       TypeOfInfection,
                       PosBloodCulture,
                       GramPosYN, GramNegYN,
                       AKI_no_early_late, ICUIntvRRT, MV)] |> names()
vars <- c("inhos_mortality","inhos_duration", risk_factors3)
single_row_vars <- c('SEX', 's_pulmonary', 's_abdominal', 's_urinary', 's_skinsoft',
                     's_other', 's_unclear', 'PosBloodCulture', 'GramPosYN',
                     'GramNegYN')

df_early_late <- df[AKI_no_early_late != "Non_AKI"]
inhos_uvr_cox_tbl <- tbl_uvregression(
  data =df_early_late[,..vars],
  method = coxph,
  # method.args = list(family=binomial),
  y = Surv(inhos_duration, inhos_mortality),
  exponentiate = T,
  hide_n = T,
  pvalue_fun = ~style_pvalue(., digits=3),
  show_single_row = c(SEX, s_pulmonary, s_abdominal, s_urinary, s_skinsoft,
                      s_other, s_unclear, PosBloodCulture, GramPosYN,
                      GramNegYN),
  estimate_fun = ~style_ratio(.,digits = 2)
) |> modify_table_styling(
  column = estimate,
  rows = !is.na(estimate),
  cols_merge_pattern = "{estimate} ({conf.low} - {conf.high})"
) |> 
  modify_header(estimate ~ "**HR (95% CI)**") |> 
  modify_column_hide(c(ci))

inhos_uv_cox_rst <- inhos_uvr_cox_tbl$table_body |> setDT()
inhos_mv_cox_vars <- inhos_uv_cox_rst[p.value < 0.2,variable] 
inhos_cox_form <- paste0("Surv(inhos_duration, inhos_mortality)","~ ", paste0(inhos_mv_cox_vars, collapse = " + "))
single_rows <- intersect(single_row_vars, inhos_mv_cox_vars)
inhos_cox_mult_fit <- coxph(as.formula(inhos_cox_form), data= df_early_late)

inhos_mult_cox_tbl <- tbl_regression(
  x = inhos_cox_mult_fit,
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

tbl_merge(list(inhos_uvr_cox_tbl, inhos_mult_cox_tbl),
          tab_spanner = c(
            "**Univariable**",
            "**Multivariable**"
          )) |> 
  # gtsummary::show_header_names()
  as_gt() |>
  gt::tab_style(
    style= list(
      cell_fill(color="yellow"),
      cell_text(color="red", weight = "bold")
    ),
    locations = cells_body(
      columns = p.value_2,
      rows = p.value_2 <0.05
    )
  ) |> 
  tab_row_group(
    label = "Bacteria",
    rows = 43:44
  ) |> 
  tab_row_group(
    label = "Source of infection",
    rows = 33:38
  ) |> 
  tab_row_group(
    label = "Laboratory finding",
    rows = 18:32
  ) |> 
  tab_row_group(
    label = "Vital sign",
    rows = 13:17
  ) |> 
  tab_row_group(
    label="Comorbidity",
    rows = 4:12
  )


# risk factor 2, Y: icu-mortality -----------------------------

## 1) logistic regression ------------------

require(gtsummary)
require(data.table)
require(gt)
risk_factors2 <- df[,.(Age,SEX, BMI,
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
                       TypeOfInfection,
                       PosBloodCulture,
                       GramPosYN, GramNegYN,
                       AKI_stage, MV)] |> names()

vars <- c("icu_mortality", risk_factors2)
single_row_vars <- c('SEX', 's_pulmonary', 's_abdominal', 's_urinary', 's_skinsoft',
                     's_other', 's_unclear', 'PosBloodCulture', 'GramPosYN',
                     'GramNegYN')
icu_uvr_lr_tbl <- tbl_uvregression(
  data = df[,..vars],
  method = glm,
  method.args = list(family=binomial),
  y = icu_mortality,
  exponentiate = T,
  hide_n = T,
  pvalue_fun = ~style_pvalue(., digits=3),
  show_single_row = c(SEX, s_pulmonary, s_abdominal, s_urinary, s_skinsoft,
                      s_other, s_unclear, PosBloodCulture, GramPosYN,
                      GramNegYN),
  estimate_fun = ~style_ratio(.,digits = 2)
) |> modify_table_styling(
  column = estimate,
  rows = !is.na(estimate),
  cols_merge_pattern = "{estimate} ({conf.low} - {conf.high})"
) |> 
  modify_header(estimate ~ "**OR (95% CI)**") |> 
  modify_column_hide(c(ci))

icu_uv_lr_rst <- icu_uvr_lr_tbl$table_body |> setDT()
icu_mv_lr_vars <- icu_uv_lr_rst[p.value < 0.2,variable] 
form <- paste0("icu_mortality","~ ", paste0(icu_mv_lr_vars, collapse = " + "))
single_rows <- intersect(single_row_vars, icu_mv_lr_vars)
multi_fit <- glm(as.formula(form), family=binomial, data=df)

mult_tbl <- tbl_regression(
  x = multi_fit,
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
  modify_header(estimate ~ "**OR (95% CI)**") |> 
  modify_column_hide(c(ci))

tbl_merge(list(icu_uvr_lr_tbl, mult_tbl),
          tab_spanner = c(
            "**Univariable**",
            "**Multivariable**"
          )
) |> 
  # gtsummary::show_header_names()
  as_gt() |>
  gt::tab_style(
    style= list(
      cell_fill(color="yellow"),
      cell_text(color="red", weight = "bold")
    ),
    locations = cells_body(
      columns = p.value_2,
      rows = p.value_2 <0.05
    )
  ) |> 
  tab_row_group(
    label = "Bacteria",
    rows = 43:44
  ) |> 
  tab_row_group(
    label = "Source of infection",
    rows = 33:38
  ) |> 
  tab_row_group(
    label = "Laboratory finding",
    rows = 18:32
  ) |> 
  tab_row_group(
    label = "Vital sign",
    rows = 13:17
  ) |> 
  tab_row_group(
    label="Comorbidity",
    rows = 4:12
  )


## 2) cox regression --------------------------------------------------------
vars <- c("icu_mortality","ICU_LOS", risk_factors2)
single_row_vars <- c('SEX', 's_pulmonary', 's_abdominal', 's_urinary', 's_skinsoft',
                     's_other', 's_unclear', 'PosBloodCulture', 'GramPosYN',
                     'GramNegYN')
icu_uvr_cox_tbl <- tbl_uvregression(
  data = df[,..vars],
  method = coxph,
  # method.args = list(family=binomial),
  y = Surv(ICU_LOS, icu_mortality),
  exponentiate = T,
  hide_n = T,
  pvalue_fun = ~style_pvalue(., digits=3),
  show_single_row = c(SEX, s_pulmonary, s_abdominal, s_urinary, s_skinsoft,
                      s_other, s_unclear, PosBloodCulture, GramPosYN,
                      GramNegYN),
  estimate_fun = ~style_ratio(.,digits = 2)
) |> modify_table_styling(
  column = estimate,
  rows = !is.na(estimate),
  cols_merge_pattern = "{estimate} ({conf.low} - {conf.high})"
) |> 
  modify_header(estimate ~ "**HR (95% CI)**") |> 
  modify_column_hide(c(ci))

icu_uv_cox_rst <- icu_uvr_cox_tbl$table_body |> setDT()
icu_mv_cox_vars <- icu_uv_cox_rst[p.value < 0.2,variable] 
form <- paste0("Surv(icu_duration, icu_mortality)","~ ", paste0(icu_mv_cox_vars, collapse = " + "))
single_rows <- intersect(single_row_vars, icu_mv_cox_vars)
icu_cox_multi_fit <- coxph(as.formula(form), data=df)

# require(gdata)
# mv("uvr_tbl","uvr_cox_tbl")
icu_mult_cox_tbl <- tbl_regression(
  x = multi_fit,
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

tbl_merge(list(uvr_cox_tbl, mult_cox_tbl),
          tab_spanner = c(
            "**Univariable**",
            "**Multivariable**"
          )
) |> 
  # gtsummary::show_header_names()
  as_gt() |>
  gt::tab_style(
    style= list(
      cell_fill(color="yellow"),
      cell_text(color="red", weight = "bold")
    ),
    locations = cells_body(
      columns = p.value_2,
      rows = p.value_2 <0.05
    )
  ) |> 
  tab_row_group(
    label = "Bacteria",
    rows = 43:44
  ) |> 
  tab_row_group(
    label = "Source of infection",
    rows = 33:38
  ) |> 
  tab_row_group(
    label = "Laboratory finding",
    rows = 18:32
  ) |> 
  tab_row_group(
    label = "Vital sign",
    rows = 13:17
  ) |> 
  tab_row_group(
    label="Comorbidity",
    rows = 4:12
  )
