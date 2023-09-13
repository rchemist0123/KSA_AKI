require(gtsummary)
require(gt)
require(data.table)

makeBaselineTable = function(data, include, dicho_vars, by){
  test = ifelse(uniqueN(df[[by]]) > 2,"aov","t.test")
  tbl1 = gtsummary::tbl_summary(
    data=data,
    by=by,
    include = include,
    missing = 'no',
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p})",
      # c(AST, ALT) ~ "{median} [{p25}–{p75}] ({min},{max})",
      c(TZ_ICUAdm_gap_hour,hosp_to_icu_duration) ~ "{median} ({p25}—{p75})"
    ),
    type = list(
    c(CFScore, SOFA_Renal) ~ "continuous",
    dicho_vars ~ 'dichotomous'
    ),
    value = list(
      SEX ~ 'Male'
      # AntbBEFCurrent ~ 1
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
  
  # if(length(unique(data[[by]]))==2) tbl1 = tbl1 |> modify_table_body(~.x |> relocate(stat_2, .after=stat_0))
  
  return(tbl1)
}

makeTable3 = function(fit){
  tbl_regression(
    x = fit,
    exponentiate = T,
    show_single_row = c(s_pulmonary, s_abdominal, s_urinary,
                        AntbBEFCurrent2, SEX,
                        Comorbidity_MOSIAC_DM, Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
                        Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, 
                        AppInitEmpThe_re),
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
}

makeSupTbl5 = function(EGFR, var){
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

suppleTbl6 = function(data){
  vars = df[,.(Age,SEX, BMI, CFScore,
               # Comorbidities
               Comorbidity_MOSIAC_DM, Comorbidity_MOSIAC_Cardio, Comorbidity_MOSIAC_Lung,Comorbidity_MOSIAC_CKD,
               Comorbidity_MOSIAC_SMT, Comorbidity_MOSIAC_HMM, Comorbidity_MOSIAC_IMM, CHF, Charlson_comorbidity_index_total, 
               AKI_stage,
               # Source of infection
               s_pulmonary, s_abdominal, s_urinary, s_skinsoft, s_other, s_unclear,
               TypeOfInfection, GramPosYN, GramNegYN,
               # Vital
               IVS_MBP, IVS_HR, IVS_BT, 
               # Lab
               Lactate, CRP_imp, Procalcitonin_imp, 
               # ArterialPH, PaCO2, PaO2, Bicarbonate,
               sofa_initial, SOFA_Renal, SOFA_except_renal, TZ_ICUAdm_gap_hour,
               
               # Characteristics of ICU D1 
               SAPS3_ICUD1, Septic_shock_ICUD1,  ICUStayInvasive_ICUD1,
               ICUStayInvasive_OU,ICUStayECMO_ICUD1,ICUStayECMO_OU)] |> names()
  tbl6 = tbl_summary(
    data= data,
    include = c(vars, tbl2_vars),
    missing='no',
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p})",
      c(TZ_ICUAdm_gap_hour,
        inhos_duration, icu_duration, invasive_MV_duration,
      ) ~ "{median} ({p25}—{p75})"
    ),
    type = list(
      c(CFScore,SOFA_Renal) ~ "continuous",
      c('s_pulmonary','s_abdominal','s_urinary',
        's_skinsoft','s_other','s_unclear',
        'GramPosYN','GramNegYN','Septic_shock_ICUD1','AntbBEFCurrent',
        'AppInitEmpThe_re') ~ "dichotomous"
    ),
    value = list(
      SEX ~ 'Male',
      c(AppInitEmpThe_re,AntbBEFCurrent) ~ 1
    ),
    
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ c(0,1)
    )
  )
  if(nrow(data[Hr1SSC_BRInf30min==1])>0){
    tbl6_1 = tbl_summary(
      data = data[Hr1SSC_BRInf30min==1],
      include = TZ_to_bolus_dur,
      statistic = list(
        TZ_to_bolus_dur ~ "{median} [{p25}-{p75}]"
      ),
      missing='no',
      digits = list(all_continuous() ~ 1,
                    all_categorical() ~ c(0,1)) 
    )
  }
 
  tbl6_2 = tbl_summary(
    data = data[Hr1SSC_ApplyVP==1],
    include = c("fluid_to_vaso_dur","TZ_to_vaso_dur","vaso_prior_fluid_yn"),
    statistic = list(
      c("fluid_to_vaso_dur","TZ_to_vaso_dur") ~ "{median} [{p25}-{p75}]",
      all_categorical() ~ "{n} ({p})"
    ),
    missing='no',
    digits = list(all_continuous() ~ 1,
                  all_categorical() ~ c(0,1)),
    value = vaso_prior_fluid_yn ~ 1
  )
  
  tbl6_3 = tbl_summary(
    data=data,
    include = LimitLifeSustain,
    missing='no',
    digits = list(all_continuous() ~ 1,
                  all_categorical() ~ c(0,1)),
    value = LimitLifeSustain ~ 1
  )
  
  if(nrow(data[Hr1SSC_BRInf30min==1])>0){
    tbl_stack(
      tbls = list(
        tbl6, tbl6_1, tbl6_2, tbl6_3
      )
    )
  } else {
    tbl_stack(
      tbls = list(
        tbl6, tbl6_2, tbl6_3
      )
    )
  }
}
suppleTbl7 = function(hospital){
  cat(hospital,'\n')
  tbl_summary(
    data = df[Center %like% substr(hospital,1,5)],
    include = c(sofa_initial, SAPS3_ICUD1,
                AKI_stage,
                inhos_mortality,
                icu_mortality,
                icu_duration, inhos_duration,
                CRRT,
                RRT_discharge_HD_PD,
                lac1h, bd1h, anti1h, br1h, app1h),
    statistic = list(
      all_categorical() ~ "{n} ({p})",
      all_continuous() ~ "{mean} ± {sd}",
      c(icu_duration, inhos_duration) ~ "{median} [{p25}-{p75}]"
    ),
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ c(0,1)
    )
  )
}

# univariable regression
# OR (HR), CI, p-val
uni_reg = function(data, vars, y, type, digit=3, time=NULL) {
  tbls = lapply(vars, \(var) {
    if(type %in% c("glm","lr","LR","logistic")){
      form = paste0(y, "~", var)
      fit = glm(as.formula(form), family=binomial, data=data)
      est = coef(fit) |> exp()
      ci = suppressMessages(confint(fit)) |> exp()  |> sapply(`[[`,1)
      p = summary(fit)$coefficients[,4]
      data.table(
        variable = var,
        `OR (95% CI)` = paste0(format(round(est, digit), nsmall=2), " (",
                                paste0(format(round(ci[2],digit),nsmall=2),"—",
                                       format(round(ci[4], digit), nsmall=2)),")"),
        p.value = ifelse(p <0.001, "<0.001", round(p,4))
      )[2,]
    } else if (type %in% c("multi","multiLR","multinomial")) {
      form = paste0(y, "~", var)
      fit = vglm(as.formula(form), multinomial(refLevel="Non-AKI"), data=data)
      smry = summaryvglm(fit)
      rst_tbl  = as.data.table(
        cbind(
          covariate = var,
          OR_CI = paste0(
                       format(round(exp(coef(fit)),2),2), " (",
                       format(round(exp(confint(fit))[,1],2),2), "—",
                       format(round(exp(confint(fit))[,2],2),2), ")"),
          pval = ifelse(smry@coef3[,4] < 0.001,"0.001", round(smry@coef3[,4],4))
        ),
        keep.rowname=T
      )[3:4,]
    } else if (type %in% c("cox","coxph")){
      if(is.null(time)) stop ("No time variable!")
      form = paste0("Surv(",time,",",y,")~",var)
      fit = coxph(as.formula(form), data=data) |> summary()
      est = fit$conf.int[1]
      ci = fit$conf.int[3:4]
      p = fit$coefficients[5]
      data.table(
        variable = var,
        `HR (95% CI)` = paste0(round(est,3), " (",paste0(round(ci[1],digit),"—",round(ci[2], digit)),")"),
        p.value = ifelse(p <0.001, "<0.001", format(round(p,4),4))
      )
    }
  })
  tbls |> setNames(vars) |> 
    do.call(what="rbind")
}

uni_reg(data=df, vars=vars, y="inhos_mortality", type="glm")
uni_reg(data=df, vars=c("SEX","Age"), y="inhos_mortality", type="coxph", time="Hospital_LOS")

# multivariable regression
multi_reg = function(data, vars, y, type, digit=2, time = NULL) {
  if (type %in%  c("glm","lr","logistic")) {
    form = paste0(y,"~",paste0(vars, collapse=" + "))
    fit = glm(as.formula(form), family=binomial, data)
    est = coef(fit) |> exp()
    ci = suppressMessages(confint(fit)) |> exp()
    p = summary(fit)$coefficients[,4]
    rst_tbl = data.table(
      variable = c("Intercept", vars),
      `OR (95% CI)` = paste0(format(round(est, digit),nsmall=2), " (",
                             paste0(format(round(ci[,1], digit),nsmall=2),"—",
                                    format(round(ci[,2], digit),nsmall=2)),")"),
      p.value = ifelse(p <0.001, "<0.001", round(p, 4))
    )
  } else if (type=="coxph") {
    form = paste0("Surv(",time,",",y,")~",paste0(vars))
    fit = coxph(as.formula(form), data)
  } else if (type %in% c("multinomial","multi")){
    form = paste0(y, "~", paste0(vars, collapse="+"))
    fit = vglm(as.formula(form), multinomial(refLevel="Non-AKI"), data=df)
    smry = summaryvglm(fit)
    coef = coef(fit)
    ci = confint(fit)
    p = smry@coef3[,4]
    coefs_no_const = coef[!names(coef) %like% "Intercept"]; print(length(coefs_no_const))
    cis_no_const = ci[!rownames(ci) %like% "Intercept",]; print(length(cis_no_const[,1]))
    ps_no_const = p[!names(p) %like% "Intercept"]; print(length(ps_no_const))
    rst_tbl  = as.data.table(
      cbind(
        OR_CI = paste0(
                       format(round(exp(coefs_no_const),2),2), " (",
                       format(round(exp(cis_no_const)[,1],2),2), "—",
                       format(round(exp(cis_no_const)[,2],2),2), ")"),
        pval = ifelse(ps_no_const < 0.001, "<0.001", round(ps_no_const,4))
      ),
      keep.rowname=TRUE
    )
  }
  return(rst_tbl)
}

makeTable2 = function(data, include, dicho_vars, by) {
  tbl2 = tbl_summary(
    data=data,
    by=by,
    include = include,
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
    modify_footnote(everything() ~ NA)
    # modify_table_body(~.x |> relocate(stat_2, .after=stat_0))
  return(tbl2)
}

tbl2_cox_reg = function(fit, single_rows=NULL) {
  tbl_regression(
    x = fit,
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
}

input_output_plot <- function(data, target=c("input","output")) {
  data[tolower(variable) %like% target, .(mean=mean(value,na.rm=T)),by=.(AKI_stage, variable)]  |> 
    ggplot(aes(x=variable, y=mean)) +
    geom_col(aes(fill=AKI_stage), position=position_dodge()) +
    scale_y_continuous(labels=scales::comma) + 
    scale_x_discrete(labels=c(1,2,3,7)) + 
    scale_fill_brewer(palette="Spectral") + 
    labs(y=paste0(tools::toTitleCase(target), " (mL)"), x="Day") +
    theme_classic() + 
    theme(legend.position = "top")
}

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
# 

