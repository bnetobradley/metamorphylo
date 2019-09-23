## recalulate licor computed values

library("tidyverse")

# calculate extra leaf related metrics which scale from basic set of measurements
recalc_licor68 <- . %>%
  select(-(A:EBSum)) %>%
  mutate(
    E = Flow*CorrFact*(H2O_s-H2O_r)/(100*S*(1000-CorrFact*H2O_s)),
    A = Flow*CorrFact*(CO2_r-CO2_s*(1000-CorrFact*H2O_r)/(1000-CorrFact*H2O_s))/(100*S),
    Ca = CO2_s - ifelse(CorrFact>1, A*S*100/(Fan*Fan_speed), 0),
    gbw = blfa_3+blfa_2*S+blfa_1*S*S,
    Rabs = (Qabs*convert),
    TleafEB = (Tair+(Rabs+2*0.95*0.0000000567*(((Tair+0)+273)^4-(Tair+273)^4)-44100*E)/(1.84*29.3*gbw+8*0.95*0.0000000567*(Tair+273)^3)),
    ## n.b. deltaTw set at 0
    TleafCnd = (1*Tleaf+0*Tleaf2+0*TleafEB),
    ## n.b. fT1, fT2 and fTeb set at 1, 0, 0
    gtw = E*(1000-(1000*0.61365*exp(17.502*TleafCnd/(240.97+TleafCnd))/(Pa+ΔPcham)+H2O_s)/2)/(1000*0.61365*exp(17.502*TleafCnd/(240.97+TleafCnd))/(Pa+ΔPcham)-H2O_s),
    gsw = 2/((1/gtw-1/gbw)+sign(gtw)*sqrt((1/gtw-1/gbw)*(1/gtw-1/gbw) + 4*K/((K+1)*(K+1))*(2*1/gtw*1/gbw-1/gbw*1/gbw))),
     gtc = 1/((K+1)/(gsw/1.6)+1/(gbw/1.37)) + K/((K+1)/(gsw/1.6) + K/(gbw/1.37)),
     Ci = ((gtc-E/2)*Ca-A)/(gtc+E/2),
    Pci = Ci*(Pa+ΔPcham)/1000,
    Pca = (CO2_s - ifelse(CorrFact>1, A*S*100/(Fan*Fan_speed), 0))*(Pa+ΔPcham)/1000,
     SVPleaf = 0.61365*exp(17.502*TleafCnd/(240.97+TleafCnd)),
   VPcham = H2O_s*(Pa+ΔPcham)/1000,
   SVPcham = 0.61365*exp(17.502*Tair/(240.97+Tair)),
   RHcham = (VPcham/SVPcham*100),
    VPDleaf = (SVPleaf-H2O_s*(Pa+ΔPcham)/1000),
    LatHFlux = (-E*44100),
    SenHFlux = 2*29.3*gbw*0.92*(Tair-TleafCnd),
    NetTherm = 2*0.95*0.0000000567*(((Tair+0)+273)^4-(TleafCnd+273)^4),
    ## n.b. deltaTw set at 0
    EBSum = Rabs+NetTherm+LatHFlux+SenHFlux,
   
   message("Constant values printed in the header of the licor output file are ignored in this recalculation.
This means that the values for deltaTw (0), fT1 (1), fT2 (0) or fTeb (0) are set as denoted in brackets.")
  )
