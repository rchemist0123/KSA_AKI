library(data.table)
require(moonBook)
require(stringr)
options(datatable.print.class = T)
load("AKI/AKI.RData")
dt <- fread('ksa5_cleaned_update.csv')

# AKI Definition ----------------------------------------------------------

# 1. baseline Cr 
# CKD : 입원 중 가장 낮은 Cr (투석을 하지 않은 기간)

# Non-CKD
## 1) 입원 중 가장 낮은 Cr (투석을 하지 않은 기간)
## 2) eGFR 75로 설정하여 계산한 SCr

set_cohort = function(baseline, eGFR=75){
  #' baseline: Cr_baseline_MDRD or Cr_baseline
  dt[,Center2 := fifelse(Center %like% '서울아산병원','서울아산병원',Center)]
  dt_icu = dt[ICUADMDateTime_ICUD1!=""] # ICU 환자 대상
  dt_icu[,ICUDischTime := fifelse(nchar(ICUDischTime)==4,
                                  paste0("0",ICUDischTime),
                                  ICUDischTime)]
  dt_icu[,ICUDischDateTime := paste(ICUDischDate, ICUDischTime)]
  
  # 1. CKD가 있는 경우 -> Cr timezero ~ Cr ICUDay7 최소값(nadir).
  # 2. CKD가 없는 경우 - >
  # - 1) Cr timezero ~ Cr IcuDay7 최솟값.
  # - 2) MDRD 
  # - 3) CKD-EPI
  dt_icu[,Cr_baseline := pmin(Cr, Cr_ICUD1, Cr_ICUD2, Cr_ICUD3, Cr_ICUD7, Cr_ICUDL, na.rm=T)]
  if(eGFR %in% c('nadir','baseline')){
    dt_icu[,Cr_baseline_MDRD := fifelse(
      Chronic_kidney_ds==1, # CKD가 있는 경우
      pmin(Cr, Cr_ICUD1, Cr_ICUD2, Cr_ICUD3, Cr_ICUD7, na.rm=T),
      # non-CKD --> MDRD식으로 계산.
      fifelse(SEX==1, 
              ((175*Age^(-0.203))/Cr_baseline)^(1/1.154),
              ((0.742*175*Age^(-0.203))/Cr_baseline)^(1/1.154)
      )
    )]
  } else if (eGFR == 'both75'){
    dt_icu[,Cr_baseline_MDRD := fifelse(SEX==1, 
                                        ((175*Age^(-0.203))/75)^(1/1.154),
                                        ((0.742*175*Age^(-0.203))/75)^(1/1.154))]
  }
  else {
    dt_icu[,Cr_baseline_MDRD := fifelse(
      Chronic_kidney_ds==1, # CKD가 있는 경우
      pmin(Cr, Cr_ICUD1, Cr_ICUD2, Cr_ICUD3, Cr_ICUD7, na.rm=T),
      # non-CKD --> MDRD식으로 계산.
      fifelse(SEX==1, 
              ((175*Age^(-0.203))/eGFR)^(1/1.154),
              ((0.742*175*Age^(-0.203))/eGFR)^(1/1.154)
      )
    )]
  }
  # baseline Cr: 10 이상 제외.
  # ICU Cr : 15 이상 제외.
  cr_except_id = dt_icu[Cr_baseline >= 10 |
                     Cr_ICUD1 >= 15 |
                     Cr_ICUD2 >= 15 |
                     Cr_ICUD3 >= 15 |
                     Cr_ICUD7 >= 15 , SubjectNo]

  data = dt_icu[!SubjectNo %in% cr_except_id]
  
  data[,inhos_mortality := AliveDead-1]
  data[,icu_mortality := AliveDeadOutcome -1]
  setnames(data, c("Elig02","ICUStayCCRT_OU"), c("Septic_shock","CRRT"))
  data[,icu_duration := difftime(ICUDischDatetime, ICUADMDateTime_ICUD1, units='days') |> as.numeric()]
  data[,inhos_duration := difftime(HosDiscDateTime, HosAdmDateTime, units='days') |> as.numeric()]
  data[,hosp_to_icu_duration := difftime(ICUADMDateTime_ICUD1, HosAdmDateTime, units='days') |> as.numeric()]
  data[,TZ_to_disch_duration := difftime(HosDiscDateTime, TZDT, units='days') |> as.numeric()]
  # 우선 AKI가 발생하지 않은 사람들은 다 제외하고 AKI 계산
  # baseline Cr과 max Cr 비교했을 때, 한 번이라도 1.5배 이상 증가하지 않은 사람들 제외

  data[,`:=`(
    AKI_initial = fcase( # fast case when
      Cr >= base * 3.0 | Cr > 4, 3,
      Cr >= base * 2.0 & Cr <= base * 2.9, 2,
      Cr %between% list(base * 1.5, base*1.9), 1,
      default = 0
    ),
    AKI_Day1 = fcase(
      Cr_ICUD1 > base * 3.0 | Cr_ICUD1 > 4 | ICUStayCCRT_ICUD1 == 1, 3,
      Cr_ICUD1 %between% list(base * 2.0, base * 2.9), 2,
      Cr_ICUD1 %between% list(base * 1.5, base * 1.9), 1,
      default = 0
    ),
    AKI_Day2 = fcase(
      Cr_ICUD2 > base * 3.0 | Cr_ICUD2 > 4 | ICUStayCCRT_ICUD2 == 1, 3,
      Cr_ICUD2 %between% list(base * 2.0, base * 2.9), 2,
      Cr_ICUD2 %between% list(base * 1.5, base * 1.9), 1,
      default = 0 # 퇴원
    ),
    AKI_Day3 = fcase(
      Cr_ICUD3 > base * 3.0 | Cr_ICUD3 > 4 | ICUStayCCRT_ICUD3 == 1, 3,
      Cr_ICUD3 %between% list(base * 2.0, base * 2.9), 2,
      Cr_ICUD3 %between% list(base * 1.5, base * 1.9), 1,
      default = 0
    )
  ), env=list(base = baseline)]

  target = data[,.SD,.SDcols=patterns('DateTime|Datetime')] |> names(); target

# AKI Day 7 정의
  data[,(target):= lapply(.SD, \(x)
                        ifelse(nchar(x) == 15,
                              paste0(substr(x, 1, 10), ' 0', substr(x, 12, 15), ':00'),
                              ifelse(nchar(x) == 16,
                                      paste0(x, ':00'),
                                      NA
                                ))),.SDcols=target]
  # difftime(ICUADMDateTime_ICUD1, HosAdmDateTime,units='days') |> as.numeric()<1
  data[,AKI_Day7 := fcase(
      Cr_ICUD7 > base * 3.0 | Cr_ICUD7 > 4 | ICUStayCCRT_ICUD7 == 1, 3,
      Cr_ICUD7 %between% list(base * 2.0, base * 2.9), 2,
      Cr_ICUD7 %between% list(base * 1.5, base * 1.9), 1,
      default = 0
    ), env=list(base = baseline)]

  target = c('AKI_initial','AKI_Day1','AKI_Day2','AKI_Day3','AKI_Day7')
  data[,AKI_stage_temp := do.call(pmax, c(.SD, na.rm=T)),.SDcols=target]
  data[,CRRT_within_7days := do.call(pmax, c(.SD, na.rm=T)),.SDcols=patterns('ICUStayCCRT_ICUD\\d')]
  data[,AKI_stage := fifelse(CRRT_within_7days==1 | AKI_stage_temp == 3, "Stage 3", 
                           fifelse(is.na(AKI_stage_temp) | AKI_stage_temp==0, "Non-AKI", 
                                   paste0("Stage ", AKI_stage_temp)))]
  
  # data[,AKI_stage := do.call(pmax, c(.SD, na.rm=T)),.SDcols=target]
  # data[,AKI_stage := fifelse(ICUStayCCRT_OU==1, "Stage 3", 
  #                      fifelse(is.na(AKI_stage)|AKI_stage==0, "Non-AKI", 
  #                            paste0("Stage ", AKI_stage))) |> as.factor()]

  # Exclude TZ minus
  data[,TZ_DateTime := fifelse(nchar(TZTM)==4, paste0(TZDT," 0",TZTM,":00"), paste0(TZDT," ",TZTM,":00"))]
  data[,TZ_ICUAdm_gap_hour := difftime(ICUADMDateTime_ICUD1, TZ_DateTime, units="hours") |> as.numeric()]

  tz_minus_id = data[TZ_ICUAdm_gap_hour < 0, SubjectNo]

  data <- data[!SubjectNo %in% tz_minus_id]
    # "Add new variables"
  data[,SEX := fifelse(SEX==1,'Male','Female')]
  data[,TypeOfInfection := ifelse(TypeInfection_MOSAICS==1,'community','nosocomial')]
  data[,`:=`(
    s_pulmonary = ifelse(SiteInfection_MOSAICS %like% '1',1,0) |> as.factor(),
    s_abdominal = ifelse(SiteInfection_MOSAICS %like% '2',1,0) |> as.factor(),
    s_urinary = ifelse(SiteInfection_MOSAICS %like% '3', 1,0) |> as.factor(),
    s_skinsoft = ifelse(SiteInfection_MOSAICS %like% '4',1,0) |> as.factor(),
    s_other = ifelse(SiteInfection_MOSAICS %like% '5|7',1,0) |> as.factor(),
    # s_neuro = ifelse(SiteInfection_MOSAICS %like% '7',1,0) |> as.factor(),
    s_unclear = ifelse(SiteInfection_MOSAICS %like% '6',1,0) |> as.factor(),
    m_bacteria = ifelse(MicrobiologyType %like% '1',1,0) |> as.factor(),
    m_virus = ifelse(MicrobiologyType %like% '2',1,0) |> as.factor(),
    m_fungus = ifelse(MicrobiologyType %like% '3',1,0) |> as.factor(),
    m_tuber = ifelse(MicrobiologyType %like% '4',1,0) |> as.factor(),
    m_other = ifelse(MicrobiologyType %like% '9',1,0) |> as.factor(),
    GramPosYN = fifelse(is.na(GramPosYN),0,GramPosYN) |> as.factor(),
    GramNegYN = fifelse(is.na(GramNegYN),0,GramNegYN) |> as.factor()
  )]
  
  data[,PosBloodCulture := fcase(is.na(BacPosSpec),0,
                              BacPosSpec %like% '1',1,
                              default=0) |> as.factor()]
  data[,AntbBEFCurrent2 := fifelse(AntbBEFCurrent==2,0,AntbBEFCurrent) |> as.factor()]
  data[,SOFA_except_renal := rowSums(.SD, na.rm=T) |> as.numeric(), .SDcols=c("SOFA_PF","SOFA_PLT","SOFA_BIL","SOFA_Cardio","SOFA_GCS") ]
  data[,AKI_YN := fifelse(AKI_stage == "Non-AKI",0,1)]
  data[,AKI_01_23 := fifelse(AKI_stage %like% "2|3",1,0)]
  # 2023-04-14 procalcitonin, CRP:  결측 시 TZ와 차이가 48시간 이내면 ICU Day1 수치 사용
  data[, Procalcitonin_imp := fifelse(is.na(Procalcitonin) & TZ_ICUAdm_gap_hour < 48, Procalcitonin_ICUD1, Procalcitonin)]
  data[, CRP_imp := fifelse(is.na(CRP) & TZ_ICUAdm_gap_hour < 48, CRP_ICUD1, CRP)]
  
  # BacPosSpec=1 이 들어간것 ("1,2" 등 포함)  & BacTestMet=1 이 들어간것
  data[,bac_blood_sepsis := fifelse(BacPosSpec %like% "1|2" & BacTestMet==1, 1,0)]
  data[,bac_blood_sepsis_gp := fifelse(bac_blood_sepsis ==1 & GramPosYN ==1,1,0)]
  data[,bac_blood_sepsis_gn := fifelse(bac_blood_sepsis == 1 & GramNegYN ==1,1,0)]
  #  FungPosSpec=1이 들어간것 ("1,2" 등 포함) & FungTestMet=1이 들어간것 
  data[,fung_blood_sepsis := fifelse(FungPosSpec  %like% "1|2" & FungTestMet==1,1,0) |> as.factor()]
  data[,BacMTYN := fifelse(is.na(BacMTYN)|BacMTYN!=1,0,1) |> as.factor()]
  # 초기 사용 항생제
  data[,glycopeptide := fifelse(InitialEmpAntibiot  %like%  "GP",1,0)]
  data[,aminoglycoside := fifelse(InitialEmpAntibiot %like% "AG",1,0)]
  data[,colistin := fifelse(InitialEmpAntibiot %like% "CL",1,0)]
  data[,gp_ag_cl := fifelse(do.call(pmax, c(.SD, na.rm=T)) == 1,1,0) ,.SDcols=c("glycopeptide","aminoglycoside","colistin")]

  # ICU 퇴실 이후 RRT 
  data[,RRT_discharge_HD_PD := fcase(ICUIntvRRTType %in% c(2,3), "Yes",
                                    default = "No")]
  data[,AppInitEmpThe_re := fifelse(AppInitEmpThe ==3,NA, 
                            fifelse(AppInitEmpThe == 2, 0,1)) |> as.factor()]
  
# Input - ouptut
  data[, input_output_ICUD1 := Input_ICUD1 - Output_ICUD1]
  data[, input_output_ICUD2 := Input_ICUD2 - Output_ICUD2]
  data[, input_output_ICUD3 := Input_ICUD3 - Output_ICUD3]
  data[, input_output_ICUD7 := Input_ICUD7 - Output_ICUD7]
  
  # Vasopressors & Inotropes Support
  data[,Vasopressor_support := do.call(pmax, c(.SD, na.rm=T)),.SDcols=patterns("Vasopressors_ICUD")]
  data[,Inotrope_support := do.call(pmax, c(.SD, na.rm=T)),.SDcols=patterns("Inotropes_ICU")]
  data[,vaso_inotrope_ICUD1 := pmax(Vasopressor_support, Inotrope_support)]
  data[,Septic_shock_ICUD1 := fifelse(Vasopressors_ICUD1==1 & Lactate_ICUD1 >=2 ,1,0) |> as.factor()]
  
  # Ventilation duration
  data[,invasive_MV_duration := difftime(
    as.Date(ICUStayInvasiveEDAT_OU,'%Y/%m/%d'),
    as.Date(ICUStayInvasiveSDAT_OU,'%Y/%m/%d'), units='days') |> as.numeric()]
  
  # TZ to treatment time
  data[Hr1SSC_ApplyVP==1 & Hr1SSC_BRInf30minDate != "" &
       Hr1SSC_BRInf30minTime != "",Fluid_init_datetime := fifelse(nchar(Hr1SSC_BRInf30minTime)==4 , 
                                                                  paste0(as.Date(Hr1SSC_BRInf30minDate,'%Y/%m/%d')," 0",Hr1SSC_BRInf30minTime,":00"), 
                                                                  paste0(as.Date(Hr1SSC_BRInf30minDate,'%Y/%m/%d')," ",Hr1SSC_BRInf30minTime,":00"))]
  data[Hr1SSC_ApplyVP==1 & Hr1SSC_ApplyVPDate != "" &
       Hr1SSC_ApplyVPTime != "",Vasopressor_inti_datetime := fifelse(nchar(Hr1SSC_ApplyVPTime)==4, 
                                                                     paste0(as.Date(Hr1SSC_ApplyVPDate, '%Y/%m/%d')," 0",Hr1SSC_ApplyVPTime,":00"), 
                                                                     paste0(as.Date(Hr1SSC_ApplyVPDate, '%Y/%m/%d')," ",Hr1SSC_ApplyVPTime,":00"))]
  data[,.(Hr1SSC_BRInf30minDate, Hr1SSC_BRInf30minTime,Fluid_init_datetime, Vasopressor_inti_datetime)]
  
  data[,fluid_to_vaso_dur := difftime(Vasopressor_inti_datetime,Fluid_init_datetime, units="hour") |> as.numeric()]
  data[,TZ_to_vaso_dur := difftime(Vasopressor_inti_datetime,TZ_DateTime, units="hour") |> as.numeric()]
  data[,TZ_to_bolus_dur := difftime(Hr1SSC_BRInf30minDatetime,TZ_DateTime, units="hour") |> as.numeric()]
  data[,vaso_prior_fluid_yn := fifelse(fluid_to_vaso_dur < 0, 1, 0) |> as.factor()]
  
  return(data)
}

df = set_cohort(baseline="Cr_baseline_MDRD") # MDRD
df2 = set_cohort(baseline="Cr_baseline") # lowest
df2[Chronic_kidney_ds==0,.N, by=AKI_stage]

df[Hr1SSC_BRInf30min==1,.(
  Hr1SSC_BRInf30minDatetime
)]

df[df2, on=.(SubjectNo), `:=` (AKI_stage_lowest = i.AKI_stage,
                              AKI_YN_lowest = i.AKI_YN)]
target = c("AKI_stage","AKI_stage_lowest")
target_new = paste0(c("AKI_stage","AKI_stage_lowest"),"_num")
df[,(target_new):=lapply(.SD,\(x) fcase(x == "Non-AKI", 0,
                                    x == "Stage 1",1,
                                    x == "Stage 2",2,
                                    default = 3)),.SDcols=target]

df[,.(gaps = mean(abs(AKI_stage_num - AKI_stage_lowest_num))), by=.(AKI_stage)]

# Non, 1-2, 3: for regression
target_new = paste0(c("AKI_stage","AKI_stage_lowest"),"_2")
df[,stage3_yn := fifelse(AKI_stage == "Stage 3",1,0)]
df[,stage3_yn_lowest := fifelse(AKI_stage_lowest == "Stage 3",1,0)]

df[,(target_new):=lapply(.SD,\(x) fifelse(x %in% c("Stage 1","Stage 2"), "Stage 12", as.character(x))),.SDcols=target]


# datetime 전처리 ------------------------------------------------------------
# 시간앞에 0 없으면 붙여주기
save.image("AKI/AKI.RData")

save.image("AKI.RData")
load("AKI.RData")


