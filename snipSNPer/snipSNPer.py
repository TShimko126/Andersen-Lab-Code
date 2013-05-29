from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
import urllib2, re, string, csv

#Get position
pos = raw_input("Enter your position/region of interest (Ex: V:9966404 or V:9966404..9976226): ")
#Change : to .. for splitting in next line
pos = re.sub(",", "", pos)
pos = re.sub(":","..",pos)
pos = pos.split("..")
chrom = pos[0]
#Determine if entered position data is a position or a region
if len(pos) <= 2:
    position = True
    reg = False
else:
    position = False
    reg = True
#Get position/region data
if position == True:
    locat = int(pos[1])
if reg == True:
    upstream = int(pos[1])
    downstream = int(pos[2])

#Get enzyme list
choice = raw_input("How would you like to choose restriction enzymes? (file/list): ")
if choice == "file":
    import tkFileDialog
    filename = tkFileDialog.askopenfilename()
    enz = open(filename, "r")
    enz = enz.readline()
    enz = re.sub(" ","",enz)
    enz = enz.split(",")
if choice == "list":    
    enz = raw_input("Which enzymes would you like to use? (Use proper capitalization, separate with commas): ")
    enz = re.sub(" ","",enz)
    enz = enz.split(",")

print " "
print "Loading restriction site information from Rebase..."
print " "

cutsites = {}

#Define function to get cut sites from Rebase (maintained by NEB)
def cutFinder(enzList):
    for i in enzList:
        dbAddress = "http://rebase.neb.com/rebase/enz/"+i+".html"
        head = {'User-Agent': 'Mozilla/5.0'}
        page = urllib2.Request(dbAddress, headers = head)
        page2 = urllib2.urlopen(page)
        text = page2.read()
        rsite = re.search(">[ATGCRYWSKMBDHVN\^]+\^?[ATGCRYWSKMBDHVN\^]+ ?<?", text)
        rsite = rsite.group()
        rsite = re.sub("\^", "", rsite)
        rsite = re.sub("[<>]+", "", rsite)
        print i + " - " + rsite
        cutsites[i] = rsite

#Call function to get cut sites
cutFinder(enz)

#Import SNP file and make first line (strain names) a list
inSNP = open("SNPsetfixed.txt","r")
SNPline = inSNP.readline()
SNPlist = SNPline.split()

#Import genome
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

print ""

#Get two strains for comparison
strain = raw_input("Input strain 1: ")
strain2 = raw_input("Input strain 2: ")
query = SNPlist.index(strain)
reference = SNPlist.index(strain2)

#Define function to get all snip-snip sites on the chromosome
def getSites():
    snipSNPs = []
    snipList = []
    SNPset = open("SNPsetfixed.txt","r")
    SNPline2 = SNPset.readline()
    SNPlist2 = SNPline2.split()
    for site in cutsites:
        rsite = ""
        cutter = cutsites[site]
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
                        context = CHR1list2[basepair-length:basepair-1]+SNPlist2[reference+1]+CHR1list2[basepair:basepair+length2]
                        newcontext = CHR1list2[basepair-length:basepair-1]+SNPlist2[query+1]+CHR1list2[basepair:basepair+length2]
                        region = CHR1list2[basepair-500:basepair+500]
                    if location2[0] == "II":      
                        context = CHR2list2[basepair-length:basepair-1]+SNPlist2[reference+1]+CHR2list2[basepair:basepair+length2]
                        newcontext = CHR2list2[basepair-length:basepair-1]+SNPlist2[query+1]+CHR2list2[basepair:basepair+length2]
                        region = CHR2list2[basepair-500:basepair+500]
                    if location2[0] == "III":      
                        context = CHR3list2[basepair-length:basepair-1]+SNPlist2[reference+1]+CHR3list2[basepair:basepair+length2]
                        newcontext = CHR3list2[basepair-length:basepair-1]+SNPlist2[query+1]+CHR3list2[basepair:basepair+length2]
                        region = CHR3list2[basepair-500:basepair+500]
                    if location2[0] == "IV":      
                        context = CHR4list2[basepair-length:basepair-1]+SNPlist2[reference+1]+CHR4list2[basepair:basepair+length2]
                        newcontext = CHR4list2[basepair-length:basepair-1]+SNPlist2[query+1]+CHR4list2[basepair:basepair+length2]
                        region = CHR4list2[basepair-500:basepair+500]
                    if location2[0] == "V":      
                        context = CHR5list2[basepair-length:basepair-1]+SNPlist2[reference+1]+CHR5list2[basepair:basepair+length2]
                        newcontext = CHR5list2[basepair-length:basepair-1]+SNPlist2[query+1]+CHR5list2[basepair:basepair+length2]
                        region = CHR5list2[basepair-500:basepair+500]
                    if location2[0] == "X":      
                        context = CHRXlist2[basepair-length:basepair-1]+SNPlist2[reference+1]+CHRXlist2[basepair:basepair+length2]
                        newcontext = CHRXlist2[basepair-length:basepair-1]+SNPlist2[query+1]+CHRXlist2[basepair:basepair+length2]
                        region = CHRXlist2[basepair-500:basepair+500]

                    #Is it in a restriction site?
                    cut1 = 0
                    cut2 = 0
                    present1 = bool(re.search(rsite, context))
                    if present1 == True:
                        cut1 = 1
                    present2 = bool(re.search(rsite, newcontext))
                    if present2 == True:
                        cut2 = 1
                    if cut1 != cut2:
                        #If it's a snip-SNP, add enzyme and position data to list
                        row = [site, chrom, basepair]
                        snipSNPs.append(row)
            #Reopen file (python prevents looping through an open file multiple times)
            SNPset = open("SNPsetfixed.txt","r")
            SNPline2 = SNPset.readline()
            SNPlist2 = SNPline2.split()
    return snipSNPs

#Call getSites and return the array of all snip-SNPs with requested enzymes    
snppr = getSites()

#Initialize lists needed for finding the closest snip-SNPs
distance = []
distanceup = []
distancedown = []

#Find the closest snip-SNPs if position is entered
if position == True:
    upst = []
    downst = []
    for element in snppr:
        dist = int(element[2]) - int(locat)
        distance.append(dist)
    for i in distance: 
        if i < 0:
            upst.append(i)
        if i > 0:
            downst.append(i)
    upbound = distance.index(max(upst))
    downbound = distance.index(min(downst))
    
    print ""
    print "The closest upstream snip-SNP is located at position " + str(snppr[upbound][1]) + ":" + str(snppr[upbound][2]) + " and is cut by " + str(snppr[upbound][0]) +"."
    print "The closest downstream snip-SNP is located at position " + str(snppr[downbound][1]) + ":" + str(snppr[downbound][2]) + " and is cut by " + str(snppr[downbound][0]) +"."        
    print ""

#Find the closest snip-SNPs if a region is entered
if reg == True:
    writer = csv.writer(open("cutsiteslist.csv", "wb"))
    writer.writerow(["Enzyme", "Position"])
    upst = []
    downst = []
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

    for i in range(5):
        upper = distanceup.index(upst2[i])
        downer = distancedown.index(downst2[i])
        writer.writerow([str(snppr[upper][0]), str(snppr[upper][2])])
        writer.writerow([str(snppr[downer][0]), str(snppr[downer][2])])
        
    
    
    #upbound = distanceup.index(max(upst))
    #downbound = distancedown.index(min(downst))
    
    
    
    #print ""
    #print "The closest upstream snip-SNP is located at position " + str(snppr[upbound][1]) + ":" + str(snppr[upbound][2]) + " and is cut by " + str(snppr[upbound][0]) +"."
    #print "The closest downstream snip-SNP is located at position " + str(snppr[downbound][1]) + ":" + str(snppr[downbound][2]) + " and is cut by " + str(snppr[downbound][0]) +"."
    #print ""
    
    








