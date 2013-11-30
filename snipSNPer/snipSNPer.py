import urllib2, re, string, csv, os
from rpy2.robjects import r as r
from EnzymeGetter import getter as getter
from Bio.Restriction import *
from Bio.Seq import *



#Get position
pos = raw_input("Enter your position/region of interest (Ex: V:9966404 or V:9966404..9976226): ")

#Change : to .. for splitting in next line
pos = re.sub(",", "", pos)
pos = re.sub(":","..",pos)
pos = pos.split("..")
chrom = pos[0]


#Write position data to csv
file1 = open("positionlist.csv", "wb")
writer1 = csv.writer(file1)
lenpos = len(pos)
for i in range(0,lenpos):
	writer1.writerow([pos[i]])
file1.close()

#Determine if entered position data is a position or a region
if len(pos) <= 2:
	position = True
	reg = False
else:
	position = False
	reg = True
	
#Get position/region data
if position:
	locat = int(pos[1])
if reg:
	upstream = int(pos[1])
	downstream = int(pos[2])

#Get enzyme list
choice = raw_input("How would you like to choose restriction enzymes? (file/list/standards): ")

if choice == "file":
	import tkFileDialog
	filename = tkFileDialog.askopenfilename()
	file = open(filename, "r")
	enz = file.readline()
	enz = re.sub(" ","",enz)
	enz = enz.split(",")
	file.close()

if choice == "list":    
	enz = raw_input("Which enzymes would you like to use? (Use proper capitalization, separate with commas): ")
	enz = re.sub(" ","",enz)
	enz = enz.split(",")


cutsites = []

#Define function to get cut sites from Rebase (maintained by NEB)
def cutFinder(enzList):
	for i in range(len(enzList)):
		dbAddress = "http://rebase.neb.com/rebase/enz/"+enzList[i]+".html"
		head = {'User-Agent': 'Mozilla/5.0'}
		page = urllib2.Request(dbAddress, headers = head)
		page2 = urllib2.urlopen(page)
		text = page2.read()
		rsite = re.search(">[ATGCRYWSKMBDHVN\^]+\^?[ATGCRYWSKMBDHVN\^]+ ?<?", text)
		rsite = rsite.group()
		rsite = re.sub("\^", "", rsite)
		rsite = re.sub("[<>]+", "", rsite)
		print enzList[i] + " - " + rsite
		cutsites.append([enzList[i], rsite])

#Call function to get cut sites
print ""
print "Getting the enzyme cut sites..."

if choice == "file" or choice == "list":
	cutFinder(enz)

if choice == "standards":
	cutsites = getter()


print ""
print "Checking for and removing isoschizomers..."

finalSites = []
seen = []
for i in range(len(cutsites)):
	if cutsites[i][1] in seen:
		continue
	else:
		seen.append(cutsites[i][1])
		finalSites.append(cutsites[i])

cutsites = finalSites



	

#Import SNP file and make first line (strain names) a list
print ""
print "Importing necessary files..."

inSNP = open("SNPsetfixed.txt","r")
SNPline = inSNP.readline()
SNPlist = SNPline.split()

#Import genome
print ""
print "Importing the genome..."

CHR1 = open("ChrInospace.txt", "r")
CHR1list = CHR1.read()
CHR2 = open("ChrIInospace.txt", "r")
CHR2list = CHR2.read()
CHR3 = open("ChrIIInospace.txt", "r")
CHR3list = CHR3.read()
CHR4 = open("ChrIVnospace.txt", "r")
CHR4list = CHR4.read()
CHR5 = open("ChrVnospace.txt", "r")
CHR5list = CHR5.read()
CHRX = open("ChrXnospace.txt", "r")
CHRXlist = CHRX.read()

CHR1list2 = re.sub(r'\s', '', CHR1list)
CHR2list2 = re.sub(r'\s', '', CHR2list)
CHR3list2 = re.sub(r'\s', '', CHR3list)
CHR4list2 = re.sub(r'\s', '', CHR4list)
CHR5list2 = re.sub(r'\s', '', CHR5list)
CHRXlist2 = re.sub(r'\s', '', CHRXlist)



#Get two strains for comparison
print ""
print "Please enter the strains you would like to use."
strain = raw_input("Input strain 1: ")
strain2 = raw_input("Input strain 2: ")
reference = SNPlist.index(strain)
query = SNPlist.index(strain2)

snipSNPs = []
snipList = []
SNPset = open("SNPsetfixed.txt","r")
SNPline2 = SNPset.readline()
SNPlist2 = SNPline2.split()

#Define function to get all snipSNP sites on the chromosome
def getSites():
	snipSNPs = []
	snipList = []
	SNPset = open("SNPsetfixed.txt","r")
	SNPline2 = SNPset.readline()
	SNPlist2 = SNPline2.split()
	for entry in range(len(cutsites)):
		rsite = ""
		cutter = cutsites[entry][1]
		for i in cutter:
			if i == "A":
				rsite += "A"
			if i == "T":
				rsite += "T"
			if i == "G":
				rsite += "G"
			if i == "C":
				rsite += "C"
			if i == "R":
				rsite += "[AG]"
			if i == "Y":
				rsite += "[CT]"
			if i == "S":
				rsite += "[GC]"
			if i == "W":
				rsite += "[AT]"
			if i == "K":
				rsite += "[GT]"
			if i == "M":
				rsite += "[AC]"
			if i == "B":
				rsite += "[CGT]"
			if i == "D":
				rsite += "[AGT]"
			if i == "H":
				rsite += "[ACT]"
			if i == "V":
				rsite += "[ACG]"
			if i == "N":
				rsite += "[ATGC]"
		length = len(rsite)
		length2 = length-1
		
	#Go through each SNP and look for snip-SNPs
		for line in SNPset:
			SNPlist2 = line.split()

	#Is there a SNP?
			if SNPlist2[query+1] != SNPlist2[reference+1]:
				location = SNPlist2[0]
				location2 = location.split('_')
				basepair = int(location2[1])

	#What's around it in the genome?
				if location2[0] == chrom:
					if location2[0] == "I":      
						context = CHR1list2[basepair-length:basepair-1]+SNPlist2[reference+1]+CHR1list2[basepair+1:basepair+length2]
						newcontext = CHR1list2[basepair-length:basepair-1]+SNPlist2[query+1]+CHR1list2[basepair+1:basepair+length2]
					if location2[0] == "II":      
						context = CHR2list2[basepair-length:basepair-1]+SNPlist2[reference+1]+CHR2list2[basepair+1:basepair+length2]
						newcontext = CHR2list2[basepair-length:basepair-1]+SNPlist2[query+1]+CHR2list2[basepair+1:basepair+length2]
					if location2[0] == "III":      
						context = CHR3list2[basepair-length:basepair-1]+SNPlist2[reference+1]+CHR3list2[basepair+1:basepair+length2]
						newcontext = CHR3list2[basepair-length:basepair-1]+SNPlist2[query+1]+CHR3list2[basepair+1:basepair+length2]
					if location2[0] == "IV":      
						context = CHR4list2[basepair-length:basepair-1]+SNPlist2[reference+1]+CHR4list2[basepair+1:basepair+length2]
						newcontext = CHR4list2[basepair-length:basepair-1]+SNPlist2[query+1]+CHR4list2[basepair+1:basepair+length2]
					if location2[0] == "V":      
						context = CHR5list2[basepair-length:basepair-1]+SNPlist2[reference+1]+CHR5list2[basepair+1:basepair+length2]
						newcontext = CHR5list2[basepair-length:basepair-1]+SNPlist2[query+1]+CHR5list2[basepair+1:basepair+length2]
					if location2[0] == "X":      
						context = CHRXlist2[basepair-length:basepair-1]+SNPlist2[reference+1]+CHRXlist2[basepair+1:basepair+length2]
						newcontext = CHRXlist2[basepair-length:basepair-1]+SNPlist2[query+1]+CHRXlist2[basepair+1:basepair+length2]

					#Is it in a restriction site?
					cut1 = 0
					cut2 = 0
					present1 = bool(re.search(rsite, context))
					if present1:
						cut1 = 1
					present2 = bool(re.search(rsite, newcontext))
					if present2:
						cut2 = 1
					if cut1 != cut2:
						#If it's a snip-SNP, add enzyme and position data to list
						row = [cutsites[entry][0], chrom, basepair, rsite, cut1, cut2, SNPlist2]
						snipSNPs.append(row)

			#Reopen file (python prevents looping through an open file multiple times)
			SNPset = open("SNPsetfixed.txt","r")
			SNPline2 = SNPset.readline()
			SNPlist2 = SNPline2.split()
	return snipSNPs

#Call getSites and return the array of all snip-SNPs with requested enzymes  
print ""
print "Looking for snipSNP sites between", strain, "and", strain2 + "..."
snppr = getSites()

#Initialize lists needed for finding the closest snip-SNPs
distance = []
distanceup = []
distancedown = []

#Find the closest snip-SNPs if position is entered
if position == True:
	print ""
	print "Looking for sites closest to the entered position..."
	file = open("cutsiteslist.csv", "wb")
	writer = csv.writer(file)
	writer.writerow(["Enzyme", "Position"])
	upst = []
	upst2 = []
	downst = []
	downst2 = []
	for element in snppr:
		dist = int(element[2]) - int(locat)
		distance.append(dist)
	for i in distance: 
		if i < 0:
			upst.append(i)
			upst2 = upst
			upst2.sort(reverse = True)
		if i > 0:
			downst.append(i)
			downst2 = downst
			downst2.sort()
		
	localSites = []
	if len(upst2) >= 5:
		for i in range(5):
			upper = distance.index(upst2[4-i])
			writer.writerow([str(snppr[upper][0]), str(snppr[upper][2])])
			localSites.append([snppr[upper][0], snppr[upper][2], snppr[upper][3], snppr[upper][4], snppr[upper][6]])
	else:
		for i in range(0, len(upst2)):
			upper = distance.index(upst2[upst2.length()-i])
			writer.writerow([str(snppr[upper][0]), str(snppr[upper][2])])
			localSites.append([snppr[upper][0], snppr[upper][2], snppr[upper][3], snppr[upper][4], snppr[upper][6]])
	if len(downst2) >= 5:
		for i in range(5):
			downer = distance.index(downst2[i])
			writer.writerow([str(snppr[downer][0]), str(snppr[downer][2])])
			localSites.append([snppr[downer][0], snppr[downer][2], snppr[downer][3], snppr[downer][4], snppr[downer][6]])
	else:
		for i in range(0, len(downst2)):
			downer = distance.index(downst2[i])
			writer.writerow([str(snppr[downer][0]), str(snppr[downer][2])])
			localSites.append([snppr[downer][0], snppr[downer][2], snppr[downer][3], snppr[downer][4], snppr[downer][6]])
	
	file.close()

#Find the closest snip-SNPs if a region is entered
if reg == True:
	print ""
	print "Looking for sites closest to the entered position..."
	file = open("cutsiteslist.csv", "wb")
	writer = csv.writer(file)
	writer.writerow(["Enzyme", "Position"])
	upst = []
	upst2 = []
	downst = []
	downst2 = []
	regSize = int(downstream) - int(upstream)
	for element in snppr:
		distup = int(element[2]) - int(upstream)
		distdown = int(element[2]) - int(downstream)
		distanceup.append(distup)
		distancedown.append(distdown)
	for i in distanceup: 
		if i < 0:
			upst.append(i)
			upst2 = upst
			upst2.sort(reverse = True)
	for i in distancedown: 
		if i > 0:
			downst.append(i)
			downst2 = downst
			downst2.sort()

	localSites = []
	if len(upst2) >= 5:
		for i in range(5):
			upper = distanceup.index(upst2[4-i])
			writer.writerow([str(snppr[upper][0]), str(snppr[upper][2])])
			localSites.append([snppr[upper][0], snppr[upper][2], snppr[upper][3], snppr[upper][4], snppr[upper][6]])
	else:
		for i in range(0, len(upst2)):
			upper = distanceup.index(upst2[upst2.length()-i])
			writer.writerow([str(snppr[upper][0]), str(snppr[upper][2])])
			localSites.append([snppr[upper][0], snppr[upper][2], snppr[upper][3], snppr[upper][4], snppr[upper][6]])
	if len(downst2) >= 5:
		for i in range(5):
			downer = distancedown.index(downst2[i])
			writer.writerow([str(snppr[downer][0]), str(snppr[downer][2])])
			localSites.append([snppr[downer][0], snppr[downer][2], snppr[downer][3], snppr[downer][4], snppr[downer][6]])
	else:
		for i in range(0, len(downst2)):
			downer = distancedown.index(downst2[i])
			writer.writerow([str(snppr[downer][0]), str(snppr[downer][2])])
			localSites.append([snppr[downer][0], snppr[downer][2], snppr[downer][3], snppr[downer][4], snppr[downer][6]])
	for i in distanceup:
		if i > 0 and i < regSize:
			inReg = distanceup.index(i)
			writer.writerow([str(snppr[inReg][0]), str(snppr[inReg][2])])
			localSites.append([snppr[inReg][0], snppr[inReg][2], snppr[inReg][3], snppr[inReg][4], snppr[inReg][6]])
	file.close()
	
#Map the snip-SNP sites in relation to the region or position of interest
cwd = os.getcwd()
r.assign("cwd", cwd)
save = raw_input("Would you like to save the chromosome map? (yes/no): ")
if save == "yes" or save == "y":
	filename = []
	filename.append(raw_input("What would you like to name the file (do not include filetype extension): "))
	print filename
	filename.append(".png")
	print filename
	filename = cwd + "/Output/" + filename
	print filename
else:
	filename = "standinfilename.png"
print ""
print "Plotting the closest sites..."
cwd = os.getcwd()
r.assign("cwd", cwd)
r.assign("filename", filename)
r("source('RPlottingScript.R')")

#Open the plot to front of screen
cd = "cd " + cwd
os.system(cd)
preview = "open -a preview " + filename
os.system(preview)

#Determine which sites will be used
print ""
use = raw_input("Which site numbers would you like to use? (separate with commas):")
filename = raw_input("What would you like to name the output file (do not include filetype extension): ")
filename = filename + ".csv"
cd = os.getcwd()
filename = cd + "/Output/" + filename
finalCSV = open(filename, "wb")
outputWriter = csv.writer(finalCSV)
outputWriter.writerow(["Enzyme", "Chromosome", "Position", "Forward Primer", "Reverse Primer" ,"Cut Strain", "Cut Band Sizes", "Uncut Strain", "Uncut Band Size(s)"])
use = re.sub(" ", "", use)
use = use.split(",")
for i in range(len(use)):
	use[i] = int(use[i])
	use[i] = use[i]-1
for site in use:
	enzymeName = localSites[site][0]
	basepair = localSites[site][1]
	recognitionSite = localSites[site][2]
	siteLocation = localSites[site][4]
	if localSites[site][3] == 0:
		cutWorm = strain
		uncutWorm = strain2
	else:
		cutWorm = strain2
		uncutWorm = strain

	if chrom == "I":      
		context = CHR1list2[basepair-700:basepair-1]+siteLocation[reference+1]+CHR1list2[basepair+1:basepair+700]
		newcontext = CHR1list2[basepair-700:basepair-1]+siteLocation[query+1]+CHR1list2[basepair+1:basepair+700]
	if chrom == "II":
		context = CHR2list2[basepair-700:basepair-1]+siteLocation[reference+1]+CHR2list2[basepair+1:basepair+700]
		newcontext = CHR2list2[basepair-700:basepair-1]+siteLocation[query+1]+CHR2list2[basepair+1:basepair+700]
	if chrom == "III":
		context = CHR3list2[basepair-700:basepair-1]+siteLocation[reference+1]+CHR3list2[basepair+1:basepair+700]
		newcontext = CHR3list2[basepair-700:basepair-1]+siteLocation[query+1]+CHR3list2[basepair+1:basepair+700]
	if chrom == "IV":
		context = CHR4list2[basepair-700:basepair-1]+siteLocation[reference+1]+CHR4list2[basepair+1:basepair+700]
		newcontext = CHR4list2[basepair-700:basepair-1]+siteLocation[query+1]+CHR4list2[basepair+1:basepair+700]
	if chrom == "V":
		context = CHR5list2[basepair-700:basepair-1]+siteLocation[reference+1]+CHR5list2[basepair+1:basepair+700]
		newcontext = CHR5list2[basepair-700:basepair-1]+siteLocation[query+1]+CHR5list2[basepair+1:basepair+700]
	if chrom == "X":      
		context = CHRXlist2[basepair-700:basepair-1]+siteLocation[reference+1]+CHRXlist2[basepair+1:basepair+700]
		newcontext = CHRXlist2[basepair-700:basepair-1]+siteLocation[query+1]+CHRXlist2[basepair+1:basepair+700]

	fill = {"region" : context, "outfile" : "out.pr3", "productSize" : "500-1200", "start" : "600", "length" : "200"}
	makeFile = '''SEQUENCE_ID=example
SEQUENCE_TEMPLATE=%(region)s
SEQUENCE_TARGET=%(start)s,%(length)s
PRIMER_TASK=pick_detection_primers
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=18
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=21
PRIMER_MAX_NS_ACCEPTED=1
PRIMER_PRODUCT_SIZE_RANGE=%(productSize)s
P3_FILE_FLAG=1
PRIMER_EXPLAIN_FLAG=1
='''

	makeFile = makeFile % fill
	fil = open("input.txt", "wb")
	fil.write(makeFile)
	fil.close()	
	command = "primer3_core -output=out.txt input.txt"
	os.system(command)

	lines = [line.strip() for line in open("out.txt", "r")]
	lines = lines[0:-1]
	output = {}
	for line in lines:
		newLine = line.split("=")
		key = newLine[0]
		value = newLine[1]
		output[key] = value

	refSeq = context
	queSeq = newcontext

	for pairNumber in range(int(output["PRIMER_PAIR_NUM_RETURNED"])):
		primerLeftKey = "PRIMER_LEFT_" + str(pairNumber)
		primerLeft = output[primerLeftKey]
		primerLeft = primerLeft.split(",")
		primerLeftSeqKey = primerLeftKey + "_SEQUENCE"
		primerLeftSeq = output[primerLeftSeqKey]
		start = int(primerLeft[0])

		primerRightKey = "PRIMER_RIGHT_" + str(pairNumber)
		primerRight = output[primerRightKey]
		primerRight = primerRight.split(",")
		primerRightSeqKey = primerRightKey + "_SEQUENCE"
		primerRightSeq = output[primerRightSeqKey]
		stop = int(int(primerRight[0]) + int(primerRight[1]))
	
		
		refSeq = Seq(context[start:stop])
		queSeq = Seq(newcontext[start:stop])
	
		rb = RestrictionBatch([enzymeName])
		refAnalysis = Analysis(rb, refSeq, linear=True)
		queAnalysis = Analysis(rb, queSeq, linear=True)
		refCutsites = refAnalysis.full()
		queCutsites = queAnalysis.full()

		print refCutsites
		print queCutsites
	
		refBandSizes = []
		refBands = refCutsites.keys()[0].catalyze(refSeq)

		refCutsites.keys()[0].catalyze(refSeq)

		for i in range(len(refBands)):
			refBandSizes.append(len(refBands[i]))
	
		refBandSizes.sort()
	


		queBandSizes = []
		queBands = queCutsites.keys()[0].catalyze(queSeq)

		queCutsites.keys()[0].catalyze(queSeq)

		for i in range(len(queBands)):
			queBandSizes.append(len(queBands[i]))
	
		queBandSizes.sort()

		
	
		if(len(queBandSizes) < len(refBandSizes)):
			cutBands = refBandSizes
			uncutBands = queBandSizes
	
		else:
			cutBands = queBandSizes
			uncutBands = refBandSizes


		outputWriter.writerow([enzymeName, chrom, basepair, primerLeftSeq, primerRightSeq, cutWorm, cutBands, uncutWorm, uncutBands])
finalCSV.close()

try:
	os.remove("standinfilename.png")
except:
	pass

os.remove("positionlist.csv")
os.remove("cutsiteslist.csv")
os.remove("input.txt")
os.remove("out.txt")


print ""
print "Final file output to " + filename
	
	
	
	





