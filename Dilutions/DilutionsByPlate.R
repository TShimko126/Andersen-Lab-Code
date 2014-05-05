##############
#You need to edit these first two lines of code to reflect your directory setup and file name
##############
#Set working directory to wherever Dilutions2.Rmd and DrugDoseReponses.csv are located
setwd("~/Dropbox/HTA/Dilutions/2014 mapping/")

#Edit this line to have your file name
concs = read.csv("20140423_RIAILmapping1_forMD.csv")
##############

#Dilution calculation
dilute = function(startConc, finConc, volume){
  (finConc/startConc)*volume
}

makePlates = function(drug, numWells, numPlates, startConc, finConc){
  stock1 = startConc*1000
  drug = paste0(as.character(drug)," - ", finConc, " uM")
  tenfolds = c(stock1/10000, stock1/1000, stock1/100, stock1/10, stock1)
  volume = (50*numWells)*numPlates
  adjVolume = 1.12*volume
  ulDrug = dilute(stock1, finConc, adjVolume)
  if (ulDrug<.25){
    tenfolds = c(stock1/10000, stock1/1000, stock1/100, stock1/10, stock1)
    for (i in seq(1,length(tenfolds))){
      newUl = dilute(tenfolds[i], finConc, adjVolume)
      if (newUl>.25 & newUl<(.1*adjVolume)){
        newStock = tenfolds[i]
      }
    }
  }
  if (exists("newStock")){
    stockDilute = stock1/newStock
    stockDilution = paste0("Use ", stock1/1000, " mM stock.")
    ulDrug = dilute(newStock, finConc, adjVolume)
    stock = newStock/1000
  } else{
    stock = stock1/1000
    stockDilution = "Already made."
  }
  ulLysateMix = adjVolume - ulDrug
  lysateAmt = ulLysateMix*.1
  smedAmt = ulLysateMix*.9
  lysateMixRec = paste0("Mix ", lysateAmt, " ul lysate with ", smedAmt, " ul S Medium. Total Volume: ", lysateAmt+smedAmt, " ul")
  addDrug = paste0("Add ", ulDrug, " ul of the ", stock, " mM stock solution of ", drug, ".")
  saveRDS(stockDilution,file="stockDilution.rds")
  saveRDS(lysateMixRec,file="lysateMixRec.rds")
  saveRDS(addDrug,file="addDrug.rds")
  saveRDS(drug,file="drug.rds")
  require(knitr) # required for knitting from rmd to md
  require(markdown) # required for md to html 
  knit('DilutionsByPlate.rmd', 'DilutionsByPlate.md') # creates md file
  date = Sys.Date()
  date = as.character(format(date, format = "%m%d%Y"))
  dirname = paste0(date,"-DiluteByPlate")
  dir.create(paste0("~/Dropbox/HTA/Dilutions/2014 mapping/", dirname), showWarnings = FALSE)
  markdownToHTML('DilutionsByPlate.md', paste0("~","/Dropbox/HTA/Dilutions/2014 mapping/",dirname,"/",gsub(" ", "-", drug), "_", finConc, ".html")) # creates html file
  madefiles = c("stockDilution.rds","lysateMixRec.rds","addDrug.rds","drug.rds")
  file.remove(madefiles)
}

for (i in seq(1,nrow(concs))){
  compound = as.character(concs[i,1])
  nwells = as.numeric(as.character(concs[i,2]))
  nplates = as.numeric(as.character(concs[i,3]))
  stconc = as.numeric(as.character(concs[i,4]))
  ficonc = as.numeric(as.character(concs[i,5]))
  makePlates(compound,nwells,nplates,stconc,ficonc)
}
















