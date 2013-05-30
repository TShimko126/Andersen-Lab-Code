miligrams = function(molweight, totalvol, milimol){
  x = milimol/1000
  x = x*molweight*totalvol
  x
}

dilute = function(startConc, finConc, vol){
  (finConc/startConc)*vol
}

makeDilute = function(dilutions, stock){
  volume = 1200
  ul = c()
  for (i in seq(1,length(dilutions))){
    dil = dilute(stock, dilutions[i], volume)
    ul[i] = dil
  }
  lysate = .1*volume
  lysate2 = rep(lysate, each = length(dilutions))
  solven = rep(solvent, each = length(dilutions))
  solv = c()
  for (i in seq(1,length(dilutions))){
    solv[i] = volume - ul[i] - lysate
  }
  concentrations = c()
  for (i in seq(1,length(dilutions))){
    concentrations[i] = dilutions[i]
  }
  df = cbind(solvent,concentrations,solv,ul,lysate2)
}

compound = readline("Which compound?: ")
solvent = readline("Which diluent?: ")
makeSol = readline("Do you need to make the stock solution? (yes/no): ")
if (makeSol == "yes"){
  makeStock = as.numeric(as.character(readline("What concentration does the stock solution need to be? (in mM): ")))
  stock = makeStock*1000
  molMass = as.numeric(as.character(readline("What is the molecular weight of the compound?: ")))
  finalVol = as.numeric(as.character(readline("What is the final volume you would like to make? (in mL): ")))
  mgs = miligrams(molMass, finalVol, makeStock)
  stockmake = paste("To make a", makeStock, "mM solution, dissolve", mgs, "mg in", finalVol, "mL of", solvent)
}
if (makeSol == "no"){
  stockmake = "Already made"
  stock = as.numeric(readline("What is stock concentration? (in mM): "))
  stock = stock*1000
}

dilutes = as.character(readline("What are the final concentrations you would like to make? (in µM): "))
dilutes = gsub(" ", "", dilutes)
dilutes = as.numeric(unlist(strsplit(dilutes, ",")))
date = Sys.Date()
date = as.character(format(date, format = "%m%d%Y"))

recipe = as.data.frame(makeDilute(dilutes,stock))
diluentAmount = c()
for (i in seq(1,length(recipe[,1]))){
  diluentAmount[i] = max(as.numeric(as.character(recipe[,4])))-as.numeric(as.character(recipe$ul[i]))
}
recipe = cbind(recipe, diluentAmount)
recipe$solv = rep(min(as.numeric(as.character(recipe$solv)))+120, length(recipe[,1]))

colnames(recipe) = c("Diluent", "   Conc. (µM)","   S med./Lysate Mix (µL)", paste("   ", compound, "at", stock/1000, "mM"), "   100 mg/mL Bact. Lysate (µL)", paste("   ", recipe$solvent[1],"(µL)"))

saveRDS(stockmake,file="stockmake.rds")
saveRDS(recipe,file="recipe.rds")
saveRDS(compound,file="compound.rds")

require(knitr) # required for knitting from rmd to md
require(markdown) # required for md to html 
knit('Dilutions2.rmd', 'Dilutions2.md') # creates md file
markdownToHTML('Dilutions2.md', paste0(compound, "_", date, ".html")) # creates html file
