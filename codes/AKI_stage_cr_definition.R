df[,.SD,.SDcols=patterns('Cr_ICUD')]

df[,.SD,.SDcols=patterns('AKI_initial|AKI_Day')]
df[,AKI_status_day1 := fcase(
  is.na(Cr_ICUD1) & is.na(Cr_ICUD2) & icu_duration < 1 & icu_mortality == 1, 'Dead',
  is.na(Cr_ICUD1) & is.na(Cr_ICUD2) & icu_duration < 1 & icu_mortality == 0, 'Alive',
  is.na(Cr_ICUD1) & icu_duration >=1, as.character(AKI_Day1),
  is.na(Cr_ICUD1) & !is.na(Cr_ICUD2), as.character(AKI_Day1),
  !is.na(Cr_ICUD1) & icu_duration >= 0, as.character(AKI_Day1),
  default = NA
)|> factor(levels = c('Dead','Alive','0','1','2','3'))]
df[,.N,AKI_status_day1]
df[is.na(AKI_status_day1),.(Cr_ICUD1,Cr_ICUD2, icu_duration, icu_mortality)]
df[AKI_status_day1 == 'Dead', .(Cr_ICUD1,Cr_ICUD2, ICUDurationDays, icu_duration, icu_mortality)]

df[,AKI_status_day2 := fcase(
  AKI_status_day1 == 'Dead','Dead',
  AKI_status_day1 == 'Alive','Alive',
  is.na(Cr_ICUD2) & icu_duration < 2 & icu_mortality == 1, 'Dead',
  is.na(Cr_ICUD2) & icu_duration < 2 & icu_mortality == 0, 'Alive',
  !is.na(Cr_ICUD2) & is.na(Cr_ICUD3) & is.na(Cr_ICUD7) & icu_mortality == 1 ,'Dead',
  !is.na(Cr_ICUD2) & is.na(Cr_ICUD3) & is.na(Cr_ICUD7) & icu_mortality == 0 ,'Alive',
  is.na(Cr_ICUD2) & icu_duration >=2, as.character(AKI_Day2),
  is.na(Cr_ICUD2) & !is.na(Cr_ICUD3) & icu_duration >=2, as.character(AKI_Day2),
  !is.na(Cr_ICUD2) & icu_duration >= 1, as.character(AKI_Day2),
  default = NA
) |> factor(levels = c('Dead','Alive','0','1','2','3'))]
df[,.N,AKI_status_day2]
df[AKI_status_day2=='Dead',.(Cr_ICUD2, Cr_ICUD3, icu_duration, icu_mortality)] |> View()
df[is.na(AKI_status_day2) ,.SD, .SDcols=patterns('Cr_ICUD')]

df[,AKI_status_day3 := fcase(
  AKI_status_day2 == 'Dead','Dead',
  AKI_status_day2 == 'Alive','Alive',
  is.na(Cr_ICUD3) & is.na(Cr_ICUD7) & icu_mortality == 1 ,'Dead',
  is.na(Cr_ICUD3) & is.na(Cr_ICUD7) & icu_mortality == 0 ,'Alive',
  is.na(Cr_ICUD3) & icu_duration < 3 & icu_mortality == 1, 'Dead',
  is.na(Cr_ICUD3) & icu_duration < 3 & icu_mortality == 0, 'Alive',
  is.na(Cr_ICUD3) & icu_duration >=3, as.character(AKI_Day3),
  is.na(Cr_ICUD3) & !is.na(Cr_ICUD7) & icu_duration >=3, as.character(AKI_Day3),
  !is.na(Cr_ICUD3) & is.na(Cr_ICUD7) & icu_mortality==0,'Alive',
  !is.na(Cr_ICUD3) & is.na(Cr_ICUD7) & icu_mortality==1,'Dead',
  !is.na(Cr_ICUD3) & icu_duration >= 2, as.character(AKI_Day3),
  default = NA
)|> factor(levels = c('Dead','Alive','0','1','2','3'))]

df[,AKI_status_day7 := fcase(
  AKI_status_day3 == 'Dead','Dead',
  AKI_status_day3 == 'Alive','Alive',
  icu_duration<7 & icu_mortality ==0, 'Alive',
  icu_duration<7 & icu_mortality == 1, 'Dead',
  !is.na(Cr_ICUD7) & icu_duration >= 7, as.character(AKI_Day7),
  default = NA
)|> factor(levels = c('Dead','Alive','0','1','2','3'))]
df[is.na(AKI_status_day7),.(AKI_Day7, Cr_ICUD3, Cr_ICUD7, icu_duration, icu_mortality)]
df[,.N,AKI_status_day7]
df[,.N,AKI_status_day3]

df[,summary(Cr_ICUD1)]
df[,Cr_new_ICUD1 := Cr_ICUD1]
df[,Cr_new_initial := Cr_baseline]
df[,Cr_new_ICUD2 := fcase(
  AKI_status_day1 %in% c('Alive','Dead'), Cr_ICUDL,
  AKI_status_day1 %in% c('Alive','Dead') & is.na(Cr_ICUDL), Cr_ICUD1,
  AKI_status_day1 %in% c('0','1','2','3'), Cr_ICUD2,
  default= NA
)]

df[,Cr_new_ICUD3 := fcase(
  AKI_status_day2 %in% c('Alive','Dead'), Cr_ICUDL,
  AKI_status_day2 %in% c('Alive','Dead') & is.na(Cr_ICUDL), Cr_ICUD2,
  AKI_status_day2 %in% c('0','1','2','3'), Cr_ICUD3,
  default= NA
)]

df[,Cr_new_ICUD7 := fcase(
  AKI_status_day3 %in% c('Alive','Dead'), Cr_ICUDL,
  AKI_status_day3 %in% c('Alive','Dead') & is.na(Cr_ICUDL), Cr_ICUD3,
  AKI_status_day3 %in% c('0','1','2','3'), Cr_ICUD7,
  default= NA
)]
df[,.SD,.SDcols=patterns('AKI_status_day\\d')]
df[,summary(Cr_new_ICUD7)]

df[AKI_status_day1 %in%  c('Alive','Dead'), .(Cr_ICUD2)]