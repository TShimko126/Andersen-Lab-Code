library(R2HTML)

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
  stockmake = paste("To make a", makeStock, "mM solution, dilute", mgs, "mg in", finalVol, "mL of", solvent, ".")
}
if (makeSol == "no"){
  stock = as.numeric(readline("What is stock concentration? (in mM): "))
  stock = stock*1000
}

dilutes = as.character(readline("What are the final concentrations you would like to make? (in µM): "))
dilutes = gsub(" ", "", dilutes)
dilutes = as.numeric(unlist(strsplit(dilutes, ",")))
date = Sys.Date()
date = as.character(format(date, format = "%m%d%Y"))

recipe = as.data.frame(makeDilute(dilutes,stock))
colnames(recipe) = c("Diluent Type", "Final Concentration (µM)","S medium (µL)", paste(compound, "at", stock/1000, "mM"), "100 mg/mL Bacterial Lysate (µL)")

print(stockmake)
print(recipe)

file = paste0(compound, "_", date)
HTML.data.frame(recipe,file=)

filename = paste0("/Users/tylershimko/Desktop/Dilutions/", compound, "_", date, ".txt")
file.create(filename)
capture.output(print(recipe, print.gap=3), file=filename)
#write.table(recipe, file = filename, sep = "\t")