df[,SEX:= fifelse(SEX==1,'Male','Female')]
df[,TypeOfInfection := ifelse(TypeInfection_MOSAICS==1,'community','nosocomial')]
df[,`:=`(
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
  GramPosYN = ifelse(is.na(GramPosYN),0, GramPosYN) |> as.factor(),
  GramNegYN = ifelse(is.na(GramNegYN),0, GramNegYN) |> as.factor()
)]
df[,inhos_mortality := AliveDead-1]
# df[,inhos_duration := fifelse(inhos_mortality == 1)]
df[,icu_mortality := AliveDeadOutcome -1]
setnames(df, c("Elig02","ICUStayCCRT_OU"),
              c("Septic_shock","CRRT"))

df[TZ_ICUAdm_gap_hour<0,.(SubjectNo, TZ_DateTime, ICUADMDateTime_ICUD1)]
df[,ICU_LOS := difftime(ICUDischDatetime, ICUADMDateTime_ICUD1, units='days') |> as.numeric()]
df[,Hospital_LOS := difftime(HosDiscDateTime, HosAdmDateTime, units='days') |> as.numeric()]
df[,median(ICU_LOS), icu_mortality]
df[ICU_LOS > 100,.(ICUADMDateTime_ICUD1, ICUDischDateTime)]
df[,PosBloodCulture := fcase(is.na(BacPosSpec),0,
                              BacPosSpec %like% '1',1,
                              default=0) |> as.factor()]
df[,SOFA_except_renal := rowSums(.SD, na.rm=T) |> as.numeric, .SDcols=c("SOFA_PF","SOFA_PLT","SOFA_BIL","SOFA_Cardio","SOFA_GCS") ]
df[,AKI_YN := fifelse(AKI_stage == "Non-AKI",0,1)]
df[,inhos_duration := difftime(HosDiscDateTime, ICUADMDateTime_ICUD1, units="days") |> as.numeric()]

df[PosBloodCulture==1, .N, keyby=.(AKI_stage, GramPosYN)][,prob := round(N/nrow(df)*100,2)][GramPosYN==1]
df[PosBloodCulture==1, .N, keyby=.(AKI_stage, GramNegYN)][,prob := round(N/nrow(df)*100,2)][GramNegYN==1]
df[PosBloodCulture==1, .N, GramPosYN]
df[PosBloodCulture==1, .N, GramNegYN]

df[PosBloodCulture==1, chisq.test(GramNegYN, AKI_stage)]
mytable(AKI_stage ~ GramPosYN, df[PosBloodCulture==1])
mytable(AKI_stage ~ GramNegYN, df[PosBloodCulture==1])

# 2023-04-14 procalcitonin, CRP:  결측 시 TZ와 차이가 48시간 이내면 ICU Day1 수치 사용
df[, Procalcitonin_imp := fifelse(is.na(Procalcitonin) & TZ_ICUAdm_gap_hour < 48, Procalcitonin_ICUD1, Procalcitonin)]
df[, CRP_imp := fifelse(is.na(CRP) & TZ_ICUAdm_gap_hour < 48, CRP_ICUD1, CRP)]

# BacPosSpec=1 이 들어간것 ("1,2" 등 포함)  & BacTestMet=1 이 들어간것
df[,bac_blood_sepsis := fifelse(BacPosSpec %like% "1|2" & BacTestMet==1, 1,0)]
df[,bac_blood_sepsis_gp := fifelse(bac_blood_sepsis ==1 & GramPosYN ==1,1,0)]
df[,bac_blood_sepsis_gn := fifelse(bac_blood_sepsis == 1 & GramNegYN ==1,1,0)]
#  FungPosSpec=1이 들어간것 ("1,2" 등 포함) & FungTestMet=1이 들어간것 
df[,fung_blood_sepsis := fifelse(FungPosSpec  %like% "1|2" & FungTestMet==1,1,0) |> as.factor()]

# 초기 사용 항생제
df[,glycopeptide := fifelse(InitialEmpAntibiot  %like%  "GP",1,0)]
df[,aminoglycoside := fifelse(InitialEmpAntibiot %like% "AG",1,0)]
df[,colistin := fifelse(InitialEmpAntibiot %like% "CL",1,0)]
df[,gp_ag_cl := fifelse(do.call(pmax, c(.SD, na.rm=T)) == 1,1,0) ,.SDcols=c("glycopeptide","aminoglycoside","colistin")]

# ICU 퇴실 이후 RRT 
df[,ICUIntvRRTType_re := fcase(ICUIntvRRTType %in% c("1,","1"),"CRRT",
                                  ICUIntvRRTType == 2, "HD", ICUIntvRRTType == 3, "PD",
                                  default = "No")]
df[!is.na(ICUIntvRRTType_re),.N,ICUIntvRRTType_re]

# with(df[!is.na(ICUIntvRRTType_re)],
#   table(ICUIntvRRTType_re,AKI_stage)
# )
# tbl_summary(
#   data =df[!is.na(ICUIntvRRTType_re)],
#   by=AKI_stage,
#   include=ICUIntvRRTType_re
# ) |> 
# add_overall()
# table(df$ICUIntvRRTType_re, df$AKI_stage)
# chisq.test(table(df$ICUIntvRRTType_re, df$AKI_stage))

# 약제 감수성
df[,AppInitEmpThe_re := ifelse(AppInitEmpThe ==3,NA, 
                          ifelse(AppInitEmpThe == 2, 0,1)) |> as.factor()]
df[, .N, AppInitEmpThe_re]
# Input - Output day 1
df[, input_output_ICUD1 := Input_ICUD1- Output_ICUD1]

# Vasopressor & Inotrope support
df[,summary(.SD),.SDcols=patterns("Vasopressors")]

df[,.N,Inotrope_support]

