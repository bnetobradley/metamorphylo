## recalulate licor computed values

library("tidyverse")

# calculate extra leaf related metrics which scale from basic set of measurements
recalc_licor68 <- . %>%
  select(-(A:TleafCnd)) %>%
  mutate(
    E = Flow*CorrFact*(H2O_s-H2O_r)/(100*S*(1000-CorrFact*H2O_s)),
    A = Flow*CorrFact*(CO2_r-CO2_s*(1000-CorrFact*H2O_r)/(1000-CorrFact*H2O_s))/(100*S),
    Ca = CO2_s - ifelse(CorrFact>1, A*S*100/(Fan*Fan_speed), 0),
    Ci = CO2_s - ifelse(CorrFact>1, A*S*100/(Fan*Fan_speed), 0)
    #df$Pci <- df$Ci*(df$Pa+df$ΔPcham)/1000
    #df$Pca <- (df$CO2_s - IF(df$CorrFact>1, df$A*df$S*100/(df$Fan*df$Fan_speed), 0))*(df$Pa+df$ΔPcham)/1000
    #df$gsw <- 2/((1/Q20-1/P20)+SIGN(Q20)*SQRT((1/Q20-1/P20)*(1/Q20-1/P20) + 4*AR20/((AR20+1)*(AR20+1))*(2*1/Q20*1/P20-1/P20*1/P20)))
    #df$gbw <- AG20+AF20*df$S+AE20*df$S*df$S
    #df$gtw <- I20*(1000-(1000*0.61365*EXP(17.502*U20/(240.97+U20))/(df$Pa+df$ΔPcham)+df$H2O_s)/2)/(1000*0.61365*EXP(17.502*U20/(240.97+U20))/(df$Pa+df$ΔPcham)-df$H2O_s)
    #df$gtc <- 1/((AR20+1)/(O20/1.6)+1/(P20/1.37)) + AR20/((AR20+1)/(O20/1.6) + AR20/(P20/1.37))
    #df$Rabs <- (AN20*AP20)
    #df$TleafEB <- (BB20+(S20+2*0.95*0.0000000567*(((BB20+$B$7)+273)^4-(BB20+273)^4)-44100*I20)/(1.84*29.3*P20+8*0.95*0.0000000567*(BB20+273)^3))
    #df$TleafCnd <- ($C$7*BC20+$D$7*BD20+$E$7*T20)
    #df$SVPleaf <- 0.61365*EXP(17.502*U20/(240.97+U20))
    #df$RHcham <- (X20/Y20*100)
    #df$VPcham <- df$H2O_s*(df$Pa+df$ΔPcham)/1000
    #df$SVPcham <- 0.61365*EXP(17.502*BB20/(240.97+BB20))
    #df$VPDleaf <- (V20-df$H2O_s*(df$Pa+df$ΔPcham)/1000)
    #df$LatHFlux <- (-I20*44100)
    #df$SenHFlux <- 2*29.3*P20*0.92*(BB20-U20)
    #df$NetTherm <- 2*0.95*0.0000000567*(((BB20+$B$7)+273)^4-(U20+273)^4)
    #df$EBSum <- S20+AC20+AA20+AB20
  )