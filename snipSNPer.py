import urllib2, re, string

#Get position
pos = raw_input("Enter your position/region of interest (Ex: V:9966404 or V:9966404..9976226): ")
pos = re.sub(":","..",pos)
pos = pos.split("..")
chrom = pos[0]
if len(pos) <= 2:
    position = True
    reg = False
else:
    position = False
    reg = True
if position == True:
    locat = pos[1]
if reg == True:
    upstream = pos[1]
    downstream = pos[2]

#Get enzyme list
enz = raw_input("Which enzymes would you like to use? (Use proper capitalization, separate with commas): ")
enz = re.sub(" ","",enz)
enz = enz.split(",")

print enz

cutsites = {}

#Define function to get cut sites from Rebase (maintained by NEB)
def cutFinder(enzList):
    for i in enzList:
        dbAddress = "http://rebase.neb.com/rebase/enz/"+i+".html"
        head = {'User-Agent': 'Mozilla/5.0'}
        page = urllib2.Request(dbAddress, headers = head)
        page2 = urllib2.urlopen(page)
        text = page2.read()
        rsite = re.search("[ATGCRYN]+\^[ATGCRYN]+", text)
        rsite = rsite.group()
        rsite = re.sub("\^", "", rsite)
        cutsites[i] = rsite

cutFinder(enz)

print cutsites

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

#Get two strains for comparison
strain = raw_input("Input strain 1: ")
strain2 = raw_input("Input strain 2: ")
query = SNPlist.index(strain)
reference = SNPlist.index(strain2)

#get ready to write the output file
#import csv
#with open(strain+"_"+strain2+".csv", 'w') as csvfile:
    #spamwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    #spamwriter.writerow(["location"] + [strain2] + [strain])
      

def getSites():
    snipSNPs = []
    snipList = []
    SNPset = open("SNPsetfixed.txt","r")
    SNPline2 = SNPset.readline()
    SNPlist2 = SNPline2.split()
    for site in cutsites:
        rsite = cutsites[site]
        length = len(rsite)
        length2 = length-1
        
    #go through each SNP and look for snip-SNPs
        for line in SNPset:
            SNPlist2 = line.split()

    #is there a SNP?
            if SNPlist2[query+1] != SNPlist2[reference+1]:
                location = SNPlist2[0]
                location2 = location.split('_')
                basepair = int(location2[1])

    #what's around it in the genome?
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

    #is it in a restriction site?
                    cut1 = 0
                    cut2 = 0
                    present1 = bool(re.search(rsite, context))
                    if present1 == True:
                        cut1 = 1
                    present2 = bool(re.search(rsite, newcontext))
                    if present2 == True:
                        cut2 = 1
                    if cut1 != cut2:
                        #####Output positions here, find closest position
                        row = [site, chrom, basepair]
                        snipSNPs.append(row)
            SNPset = open("SNPsetfixed.txt","r")
            SNPline2 = SNPset.readline()
            SNPlist2 = SNPline2.split()
    return snipSNPs


    
snppr = getSites()
print snppr