#Set working directory to wherever Dilutions2.Rmd and DrugDoseReponses.csv are located
#setwd("/Users/tylershimko/Desktop/Dose_Response_Dilution_Calculator")


#Function to determine miligrams necessary to make stock solution
miligrams = function(molweight, totalvol, milimol){
  x = milimol/1000
  x = x*molweight*totalvol
  x
}

#Dilution calculation
dilute = function(startConc, finConc, volume){
  (finConc/startConc)*volume
}

#Function to make set of necessary dilutions 
makeDilute = function(dilutions, stock, volume){
  ul = c()
  #Caluculate the dilution for each of the requested dilutions
  for (i in seq(1,length(dilutions))){
    dil = dilute(stock, dilutions[i], volume)
    ul[i] = dil
  }
  #The total amount of lysate solution should be 10% of the total volume
  lysate = .1*volume
  #Repeat the amount of lysate, as it will be necessary tas the final column in 
  lysate2 = rep(lysate, each = length(dilutions))
  #
  solv = c()
  for (i in seq(1,length(dilutions))){
    solv[i] = volume - ul[i] - lysate
  }
  #Get a vector of the final concentrations
  concentrations = c()
  for (i in seq(1,length(dilutions))){
    concentrations[i] = dilutions[i]
  }
  #Make a data frame of the type of solvent used (DMSO, water, etc.), the final concentrations, the amount of stock drug/compound solution to added in uL, and the total amount of lysate added for that dilution step
  df = cbind(solvent,concentrations,solv,ul,lysate2)
}

#Get the compound names, molecular weights, solvents, and dilution curves from a csv file entitled "DrugDoseResponses.csv"
curves = read.csv("DrugDoseResponses.csv")

#Clean up the curves data so that only rows with complete, valid entries are used
curves2 = curves[curves[,5] != "" & curves[,3] != "NA" & curves[,4] != "",]
curves2 = curves2[is.na(curves2[,1]) == FALSE,]

volume = as.numeric(as.character(readline("What volume would you like to use? (in uL): ")))

#Initialize 
for (i in seq(1,nrow(curves2))){
  #Pull compound and solvent info from csv
  compound = as.character(curves2[i,1])
  solvent = as.character(curves2[i,3])
  ####Next section is an artifact of an older version of this code that calculated the amount of powder necessary to make a stock solution
  makeSol = "no"
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
    stock = as.numeric(as.character(curves2[i,4]))
    stock1 = stock*1000
  }
  ####End artifact
  
  #Read in the dilution curves and save as a numeric vector
  dilutes = as.character(as.character(curves2[i,5]))
  dilutes = gsub(" ", "", dilutes)
  dilutes = as.numeric(unlist(strsplit(dilutes, ",")))
  
  #Get the date
  date = Sys.Date()
  date = as.character(format(date, format = "%m%d%Y"))
  
  #This section checks a series of ten-fold dilutions to see which stock concentration is most effecient at creating a dilution curve that is "pipettable" (not pipetting less than .5 ul), and not above 1% diluent so as to control for diluent effects
  tenfolds = c(stock1/10000, stock1/1000, stock1/100, stock1/10, stock1)
  count = c()
  for (j in seq(1,length(tenfolds))) {
    #Create a recipe data frame for each of the 10 fold dilutions listed above
    recipe = as.data.frame(makeDilute(dilutes,tenfolds[j],volume))
    #Count the total dilution amounts and check how many fall in our magic range of 0.5<x<12, 12 is 1% of the 1200 uL volume we are making up
    count[j] = sum(as.numeric(as.character(recipe[,4])) >=.5 & as.numeric(as.character(recipe[,4])) <= 12)
  }
  #Check if all of the counts are zero, if this is true, you cannot make of the dilutions with the stock concentration you already have
  if (all(count == 0L) == TRUE){
    #Sets brake variable equal to TRUE, this is used in the markdown file to tell us that our stock is not concentrated enough for this curve
    brake = TRUE
  } else {
    brake = FALSE
  }
  #Determine which dilution has the highest number of dilutions that can be made
  index = which(count == max(count))
  stock2 = tenfolds[index]
  #If multiple stocks can make the same number of dilutions, take the lower concentration stock, this is arbitrary, as anyone could be used. However, using the lowest working stock controls pretty well for only needing to use a higher stock later on.
  stock2 = stock2[1]
  
  #If the stock you need is more dilute than the one you have, this will tell you how to dilute the one you have
  if (stock2 != stock1){
    stockmake = paste0("You must make a ", stock1/stock2, "-fold dilution of the ", stock, " mM stock that is already made.")
  }
  
  #Make the recipe
  recipe = as.data.frame(makeDilute(dilutes,stock2, volume))
  
  goodrows = c()
  tooHigh = c()
  #Check which rows have under 1% diluent and which have over
  for (k in seq(1,nrow(recipe))){
    if (as.numeric(as.character(recipe[k,4])) <= 12 & as.numeric(as.character(recipe[k,4])) >= .2 | as.numeric(as.character(recipe[k,4])) == 0){
      goodrows = append(goodrows, as.numeric(k))
    }
    if (as.numeric(as.character(recipe[k,4])) > 12){
      tooHigh = append(tooHigh, as.numeric(k))
    }
  }
  
  #Separate the rows
  high = recipe[tooHigh,]
  recipe = recipe[goodrows,]
  
  #Get the concentrations that don't fall under 1%
  higherConcs = as.numeric(as.character(high[,2]))
  
  #Check if you have concentrations that don't fall under 1%
  if (length(higherConcs) > 0){
    higherStocks = c()
    currentStock = stock2
    #For every step between your working stock and concentrated stock, ceoncetrate 10-fold
    while (currentStock <= stock1){
      higherStocks = append(higherStocks, currentStock)
      currentStock = currentStock*10
    }
    #For each of the higher concentrated stocks, see if all of the more concentrated solutions fall in the magic range
    workStocks = c()
    for (l in seq(1,length(higherStocks))){
      highRecipe = as.data.frame(makeDilute(higherConcs,higherStocks[l],volume))
      highRecipe[,2:5] = sapply(highRecipe[,2:5],as.character)
      highRecipe[,2:5] = sapply(highRecipe[,2:5],as.numeric)
      if (all(highRecipe[,4] >= .2) & all(highRecipe[,4] <= 12)){
        workStocks = append(workStocks, higherStocks[l])
      }
    }
    
    #Of the working stocks (ones that fit the above criteria), take the most dilute
    finalHigherStock = min(workStocks)
    
    #Make the recipe for the higher concentrations
    highRecipe = as.data.frame(makeDilute(higherConcs, finalHigherStock, volume))
    for (row in seq(1,nrow(highRecipe))){
      for (j in seq(2,length(highRecipe))){
        highRecipe[row,j] = as.numeric(as.character(highRecipe[row,j]))
      }
    }
  }
  
  
  #Data frame finalization for recipe
  diluentAmount = c()
  for (i in seq(1,length(recipe[,1]))){
    diluentAmount[i] = max(as.numeric(as.character(recipe[,4])))-as.numeric(as.character(recipe$ul[i]))
  }
  recipe = cbind(recipe, diluentAmount)
  recipe$solv = rep(min(as.numeric(as.character(recipe$solv)))+120, length(recipe[,1]))
  recipe[,4] = as.numeric(as.character(recipe[,4]))
  colnames(recipe) = c("Dil.", "Conc. (µM)","S med./Lys Mix (µL)", paste(compound, "at", stock2/1000, "mM"), "100 mg/mL Bact. Lys. (µL)", paste(recipe$solvent[1],"(µL)"))
  
  #Finalize the data frame for the recipe requiring the higher concentration stock
  if (exists("highRecipe")){
    diluentAmount = c()
    for (i in seq(1,nrow(highRecipe))){
      diluentAmount[i] = max(as.numeric(as.character(recipe[,4])))-as.numeric(as.character(highRecipe[i,4]))
    }
    highRecipe = cbind(highRecipe, diluentAmount)
    highRecipe[,3] = rep(min(as.numeric(as.character(recipe[,3]))), length(highRecipe[,1]))
    highRecipe[,4] = as.numeric(as.character(highRecipe[,4]))
    colnames(highRecipe) = c("Dil.", "Conc. (µM)","S med./Lys Mix (µL)", paste(compound, "at", finalHigherStock/1000, "mM"), "100 mg/mL Bact. Lys. (µL)", paste(recipe[1,1],"(µL)"))
  }
  
  #Create temporary files to hand off to the markdown file
  saveRDS(stockmake,file="stockmake.rds")
  saveRDS(recipe,file="recipe.rds")
  saveRDS(compound,file="compound.rds")
  saveRDS(brake,file="brake.rds")
  if (exists("highRecipe")){
    saveRDS(highRecipe,file="highRecipe.rds")
  }
  
  #Knit to HTML and save file
  require(knitr) # required for knitting from rmd to md
  require(markdown) # required for md to html 
  knit('Dilutions2.rmd', 'Dilutions2.md') # creates md file
  markdownToHTML('Dilutions2.md', paste0("~","/Desktop/Dilutions/",gsub(" ", "-", compound), "_", date, ".html")) # creates html file
  madefiles = c("compound.rds", "Dilutions2.md", "recipe.rds", "stockmake.rds", "brake.rds")
  if (exists("highRecipe")){
    madefiles = append(madefiles, "highRecipe.rds")
  }
  #Remove the temporary files
  file.remove(madefiles)
  
  #Remove conditionally made variables that may interfere with with the functioning of the next iteration of the loop
  f = c("makeSol", "count", "stock1", "stock2", "recipe", "higherConcs")
  if (exists("highRecipe")){
    f = append(f, "highRecipe")
  }
  rm(list = f)
}