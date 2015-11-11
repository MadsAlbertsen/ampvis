#' Make raw MiDAS function data into ampvis format
#'
#' Make raw MiDAS function data into ampvis format
#'
#' @usage amp_cleanMiF(data)
#'
#' @param data (required) A data frame with MiDAS functions.
#' 
#' @export
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_cleanMiF <- function(data){
  
  MiF <- mutate(data,
                FIL = Filamentous.morphology,
                AOB = paste(AOB.Other,AOB.In.situ),
                AOB = ifelse(AOB %in% c("POS POS","NEG POS","NT POS","POS NT", "POS VAR", "VAR POS"), "POS", AOB),
                AOB = ifelse(AOB %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", AOB),
                AOB = ifelse(AOB %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", AOB),
                AOB = ifelse(AOB == "NT NT", "NT", AOB),
                NOB = paste(NOB.Other,NOB.In.situ),
                NOB = ifelse(NOB %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", NOB),
                NOB = ifelse(NOB %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", NOB),
                NOB = ifelse(NOB %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", NOB),               
                NOB = ifelse(NOB == "NT NT", "NT", NOB),
                Anammox = paste(Anammox.Other,Anammox.In.situ),
                Anammox = ifelse(Anammox %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", Anammox),
                Anammox = ifelse(Anammox %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", Anammox),
                Anammox = ifelse(Anammox %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", Anammox),               
                Anammox = ifelse(Anammox == "NT NT", "NT", Anammox),
                AU.MIX = paste(Autotroph.mixotroph.Other,Autotroph.mixotroph.In.situ),
                AU.MIX = ifelse(AU.MIX %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", AU.MIX),
                AU.MIX = ifelse(AU.MIX %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", AU.MIX),
                AU.MIX = ifelse(AU.MIX %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", AU.MIX),               
                AU.MIX = ifelse(AU.MIX == "NT NT", "NT", AU.MIX),
                PAO = paste(PAO.Other,PAO.In.situ),
                PAO = ifelse(PAO %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", PAO),
                PAO = ifelse(PAO %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", PAO),
                PAO = ifelse(PAO %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", PAO),               
                PAO = ifelse(PAO == "NT NT", "NT", PAO),
                GAO = paste(GAO.Other,GAO.In.situ),
                GAO = ifelse(GAO %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", GAO),
                GAO = ifelse(GAO %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", GAO),
                GAO = ifelse(GAO %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", GAO),               
                GAO = ifelse(GAO == "NT NT", "NT", GAO),
                HET = paste(Aerobic.heterotroph.Other,Aerobic.heterotroph.In.situ),
                HET = ifelse(HET %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", HET),
                HET = ifelse(HET %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", HET),
                HET = ifelse(HET %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", HET),               
                HET = ifelse(HET == "NT NT", "NT", HET),
                NO2.R = paste(Nitrite.reduction.Other,Nitrite.reduction.In.situ),
                NO2.R = ifelse(NO2.R %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", NO2.R),
                NO2.R = ifelse(NO2.R %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", NO2.R),
                NO2.R = ifelse(NO2.R %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", NO2.R),               
                NO2.R = ifelse(NO2.R == "NT NT", "NT", NO2.R),   
                NO3.R = paste(Nitrate.reduction.Other,Nitrate.reduction.In.situ),
                NO3.R = ifelse(NO3.R %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", NO3.R),
                NO3.R = ifelse(NO3.R %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", NO3.R),
                NO3.R = ifelse(NO3.R %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", NO3.R),               
                NO3.R = ifelse(NO3.R == "NT NT", "NT", NO3.R),
                AN.Other = paste(Other.anaerobic.Other,Other.anaerobic.In.situ),
                AN.Other = ifelse(AN.Other %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", AN.Other),
                AN.Other = ifelse(AN.Other %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", AN.Other),
                AN.Other = ifelse(AN.Other %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", AN.Other),               
                AN.Other = ifelse(AN.Other == "NT NT", "NT", AN.Other),   
                FER = paste(Fermentation.Other,Fermentation.In.situ),
                FER = ifelse(FER %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", FER),
                FER = ifelse(FER %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", FER),
                FER = ifelse(FER %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", FER),               
                FER = ifelse(FER == "NT NT", "NT", FER)             
  )
  return(MiF)
}
