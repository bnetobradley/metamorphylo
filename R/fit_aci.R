#install.packages("plantecophys")
library(plantecophys)
library(bayCi)
library(units)
library(magrittr)
source(file = "R/correct_racir.R")

## Fit an A-Ci curve with plantecophy package
## racir data needs to be pre-corrected for this

# read in area corrected data file
areacorracir <- read_csv("data/areacorrectracirfiles.csv")
areacorracir <- areacorracir %>% select(-X1)

## read in dataframe with data curves and empty curves matched by time
empdatmatch <- read.csv("data/emptydatamatch.csv")
empdatmatch <- empdatmatch %>% select(numb, unique_id)

## get rid of empty and data curves that have less than 330 observations 
## start with 170349 observations
unq <- unique(areacorracir$id)
removefile <- c()
for (i in 1:length(unq)) {
  data <- areacorracir %>% filter(id == unq[i])
  ifelse((nrow(data) >= 300),
         removefile <- NA, 
         removefile <- print(unq[i]))
  ifelse((is.na(removefile)),
         areacorracir <- areacorracir , 
    areacorracir <- areacorracir[-which(areacorracir$id == unq[i]),])
}
## end up with 169654 observations


## Fit A-Ci curves using ecophys package
pep_fit_bi_cap <- data.frame(matrix(nrow = 430, ncol = 4))
colnames(pep_fit_bi_cap) <- c("Vcmax", "Jmax", "Rd", "UniqueID")
for (i in 1:nrow(empdatmatch)){
  tryCatch({
    empty <- areacorracir %>% 
      filter(id == empdatmatch$numb[i]) %>% 
      mutate(
        A  = set_units( A,  umol / m^2 / s ),
        Cr = set_units( CO2_r, umol / mol     ),
        time = set_units(time, s)
      ) %>% 
      bayCi:::prepare_empty() %>%
      dplyr::mutate_if(~inherits(.x, "units"), units::drop_units) 
    empty %<>% filter(empty$use == TRUE) %>% filter(obs %in% 45:285)
    notempty <- areacorracir %>% filter(id == empdatmatch$unique_id[i]) %>% filter(obs %in% empty$obs)
    corrected <- correct_racir(empty = empty, file = notempty)
    fitecophys <- fitaci(data = corrected, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", PPFD = "Q"), fitmethod = "bilinear")
    # open file for saving plot
    mypath <- paste("data/ecophysplots/", empdatmatch$unique_id[i], ".jpg", sep = "")
    jpeg(file=mypath)
    plot(fitecophys, sub = empdatmatch$unique_id[i])
    dev.off()
    pep_fit_bi_cap[i, 1] <- (fitecophys$pars)[1,1]
    pep_fit_bi_cap[i, 2] <- (fitecophys$pars)[2,1]
    pep_fit_bi_cap[i, 3] <- (fitecophys$pars)[3,1]
    pep_fit_bi_cap[i, 4] <- paste(empdatmatch$unique_id[i])
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#write_csv(pep_fit, "data/ecophysfit.csv")

######
###### check data quality

funky <- c("2019-05-08-0842_luzula_sylvatica",
"2019-05-08-1030_helwingia_chinensis",
"2019-05-09-0740_dryopteris_dilatata",
"2019-05-09-0958_dryopteris_carthusiana",
"2019-05-09-1029_dryopteris_goldiana",
"2019-05-10-0825_iris_wilsonii",
"2019-05-10-0915_stachyurus_chinensis",
"2019-05-10-0951_trochodendron_aralioides",
"2019-05-10-1009_stewartia_serrata",
"2019-05-11-0717_dryopteris_erythrosora",
"2019-05-11-0752_drimys_winteri",
"2019-05-11-0837_sciadopitys_verticillata",
"2019-05-11-1007_eurya_japonica",
"2019-05-11-1043_rhododendron_hyperythrum",
"2019-05-13-0914_narcissus_minor",
"2019-05-22-0802_rhododendron_morii",
"2019-05-23-1031_gymnocarpium_dryopteris",
"2019-05-26-0739_gentiana_acaulis",
"2019-05-26-0756_hypericum_densiflorum",
"2019-05-26-1041_equisetum_scirpoides",
"2019-05-27-1005_dryopteris_expansa",
"2019-06-01-1017_polypodium_vulgare",
"2019-06-02-0934_cornus_kousa",
"2019-06-03-0722_helwingia_chinensis",
"2019-06-03-0929_cercidiphyllum_japonicum",
"2019-06-06-0724_daphniphyllum_macropodum",
"2019-06-06-0746_stewartia_serrata",
"2019-06-06-0819-trochodendron_aralioides",
"2019-06-06-0850_matteuccia_struthiopteris",
"2019-06-06-0909_stachyurus_chinensis",
"2019-06-06-0927_iris_wilsonii",
"2019-06-08-0723_rhododendron_hyperythrum",
"2019-06-09-0728_quercus_macranthera",
"2019-06-09-0946_cephalotaxus_harringtonia",
"2019-06-09-1031_dryopteris_erythrosora",
"2019-06-11-0810_molinia_caerulea",
"2019-06-12-0955_halocarpus_bidwilii",
"2019-06-13-0737_pinus_cembra",
"2019-06-13-0828_chamaedaphne_calyculata",
"2019-06-13-0906_dianthus_microlepis",
"2019-06-13-0922_gentiana_acaulis",
"2019-06-13-1008_abies_pinsapo",
"2019-06-14-0718_arachniodes_simplicior",
"2019-06-16-1026_tsuga_sieboldii",
"2019-06-17-0824_scirpus_microcarpus",
"2019-06-17-0954_syringa_oblata",
"2019-06-18-0704_liquidambar_formosana",
"2019-06-18-0746_polystichum_setiferum",
"2019-06-18-0838_microbiota_decussata",
"2019-06-18-0902_tsuga_mertensiana",
"2019-06-18-0955_acorus_americanus",
"2019-06-20-0836_rhododendron_barbatum",
"2019-06-22-0758_iris_wilsonii",
"2019-06-22-0900_trochodendron_aralioides",
"2019-06-25-0704_glandora_diffusa",
"2019-06-25-0721_gentiana_acaulis",
"2019-06-25-0721_hypericum_densiflorum",
"2019-06-25-0942_halocarpus)bidwilii",
"2019-07-01-1030_dianthus_microlepis",
"2019-07-03-0853_dryopteris_carthusiana",
"2019-07-03-0931_arachniodes_simplicior",
"2019-07-04-0955_sciadopitys_verticillata",
"2019-07-08-0802_oxalis_oregana",
"2019-07-08-0943_ephedra_frustillata",
"2019-07-09-0834_adiantum_pedatum",
"2019-07-09-0952_molinia_caerulea",
"2019-07-16-0720_picea_alcoquiana",
"2019-07-16-1042_equisetum_scirpoides",
"2019-07-18-0927_sarracenia_minor",
### less funky
"2019-05-08-0903_dryopteris_wallichiana",
"2019-05-09-0935_dicksonia_antarctica",
"2019-05-09-1054_cyperipedium_japonicum",
"2019-05-11-1101_clethra_delavayi",
"2019-05-23-0824_rhododendron_campylocarpum",
"2019-05-26-0954_tolmeia_menziesii",
"2019-05-30-0659_adiantum_pedatum",
"2019-05-31-0801_quercus_macranthera",
"2019-06-03-0758_tapiscia_sinensis",
"2019-06-03-1031_rhododendron_morii",
"2019-06-08-0851_actinidia_kolomikta",
"2019-06-09-0807_ampelopsis_megalophylla",
"2019-06-11-0937_symphoricarpus_mollis",
"2019-06-11-1004_paulownia_tormentosa",
"2019-06-11-1045-adiantum_pedatum",
"2019-06-13-1056_fuchsia_magellanica",
"2019-06-15-0720_dryopteris_carthusiana",
"2019-06-15-0903_rhododendron_campylocarpum",
"2019-06-17-0755_sarracenia_minor",
"2019-06-18-0731_dryopteris_expansa",
"2019-06-18-0806_picea_alcoquiana",
"2019-06-18-0919_tolmeia_menziesii",
"2019-06-18-0939_myrica_gale",
"2019-06-19-1054_pinus_sylvestris",
"2019-06-20-1045_staphylea_bumalda",
"2019-06-22-0827_stachyurus_chinensis",
"2019-06-22-0929_stewartia_serrata",
"2019-06-22-1015_drimys_winteri",
"2019-06-25-0800_abies_pinsapo",
"2019-06-25-1037_berberis_darwinii",
"2019-07-01-0803_tsuga_sieboldii",
"2019-07-03-0821_dryopteris_goldiana",
"2019-07-04-0736_eurya_japonica",
"2019-07-04-1044_quercus_macranthera",
"2019-07-16-0929_glandora_diffusa")

## filter out RACiR curves that appear mis-calibrated or display non-unidirectional responses 
funkfiltered <- pep_fit_bi_cap[which(!(pep_fit_bi_cap$UniqueID %in% funky)),]
funkfiltered <- funkfiltered %>% filter(Vcmax > 0) %>% filter(Jmax > 0) %>% filter(Jmax < 500)

#write_csv(funkfiltered, "data/ecophysgenerousfit.csv")
