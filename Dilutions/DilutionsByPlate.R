##############
#You need to edit these first two lines of code to reflect your directory setup and file name
##############
#Set working directory to wherever Dilutions2.Rmd and DrugDoseReponses.csv are located
setwd("~/Dropbox/HTA/Dilutions/2014 mapping/")


#Edit this line to have your file name
concs = read.csv("20140423_RIAILmapping1_forMD.csv")
##############


a = c(1,2,3,4,5)

concs[,6] = rep("DMSO", nrow(concs))
concs[,7] = sample(a,23,replace = TRUE)

#Dilution calculation
dilute = function(startConc, finConc, volume){
    (finConc/startConc)*volume
}

makePlates = function(drug, numWells, numPlates, startConc, finConc, diluent, concentration){
    stock1 = startConc*1000
    drug1= paste0(as.character(drug)," - ", finConc, " µM")
    tenfolds = c(stock1/10000, stock1/1000, stock1/100, stock1/10, stock1)
    volume = (50*numWells)*numPlates
    adjVolume = 1.12*volume
    ulDrug = dilute(stock1, finConc, adjVolume)
    if (ulDrug<.25){
        tenfolds = c(stock1/10000, stock1/1000, stock1/100, stock1/10, stock1)
        for (i in seq(1,length(tenfolds))){
            newUl = dilute(tenfolds[i], finConc, adjVolume)
            if (newUl>.25 & newUl<((concentration/100)*adjVolume)){
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
        stockDilution = paste0("Use ", stock, " mM stock.")
    }
    ulDrug = round(ulDrug, 2)
    ulDiluent = round(((concentration/100)*adjVolume) - ulDrug, 2)
    ulLysateMix = adjVolume - (ulDrug + ulDiluent)
    lysateAmt = round(ulLysateMix*.1, 0)
    smedAmt = round(ulLysateMix*.9, 0)
    lysateMixRec = paste0(lysateAmt+smedAmt, " µl")
    addDrug = paste0("Add ", round(ulDrug, 2), " ul of ", drug, " and ", ulDiluent, " µL ", diluent, ".")
    dilConc = paste0(concentration, "%")
    saveRDS(stockDilution,file="stockDilution.rds")
    saveRDS(lysateMixRec,file="lysateMixRec.rds")
    saveRDS(addDrug,file="addDrug.rds")
    saveRDS(drug1,file="drug.rds")
    saveRDS(dilConc, "dilConc.rds")
    require(knitr) # required for knitting from rmd to md
    require(markdown) # required for md to html 
    knit('DilutionsByPlate.rmd', 'DilutionsByPlate.md') # creates md file
    date = Sys.Date()
    date = as.character(format(date, format = "%m%d%Y"))
    dirname = paste0(date,"-DiluteByPlate")
    dir.create(paste0("~/Dropbox/HTA/Dilutions/2014 mapping/", dirname), showWarnings = FALSE)
    markdownToHTML('DilutionsByPlate.md', paste0("~","/Dropbox/HTA/Dilutions/2014 mapping/",dirname,"/",gsub(" ", "-", drug), "_", finConc, ".html")) # creates html file
    madefiles = c("stockDilution.rds","lysateMixRec.rds","addDrug.rds","drug.rds", "dilConc.rds", 'DilutionsByPlate.md')
    file.remove(madefiles)
}

for (i in seq(1,nrow(concs))){
    compound = as.character(concs[i,1])
    nwells = as.numeric(as.character(concs[i,2]))
    nplates = as.numeric(as.character(concs[i,3]))
    stconc = as.numeric(as.character(concs[i,4]))
    ficonc = as.numeric(as.character(concs[i,5]))
    diluent = concs[i,6]
    concentration = concs[i,7]
    makePlates(compound,nwells,nplates,stconc,ficonc,diluent,concentration)
}
















