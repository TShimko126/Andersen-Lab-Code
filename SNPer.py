#Get restriction site from Rebase (maintained by NEB)
import urllib2, re
def cutFinder():
    enzyme = raw_input("Which enzyme would you like to use? (Use proper capitalization, no spaces): ")
    dbAddress = "http://rebase.neb.com/rebase/enz/"+enzyme+".html"
    head = {'User-Agent': 'Mozilla/5.0'}
    page = urllib2.Request(dbAddress, headers = head)
    page2 = urllib2.urlopen(page)
    text = page2.read()
    rsite = re.search("[ATGCRYN]+\^[ATGCRYN]+", text)
    rsite = rsite.group()
    rsite = re.sub("\^", "", rsite)
    return rsite

rsite = cutFinder()

print rsite

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
import re
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

#Get restriction site
import string
rsite2 = rsite.translate(string.maketrans('TAGC', 'ATCG'))
rsite3 = rsite2[::-1]

print rsite3

#Get ready to write the output file
import csv
with open(strain+"_"+strain2+"_"+enzyme+".csv", 'w') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow(["location"] + [strain2] + [strain])
    
#Go through each SNP and look for snip-SNPs
for line in inSNP:
    SNPlist2 = line.split()

#Is there a SNP?
    if SNPlist2[query+1] != SNPlist2[reference+1]:
        location = SNPlist2[0]
        location2 = location.split('_')
        basepair = int(location2[1])

#What's around it in the genome?
        if location2[0] == "I":      
            context = CHR1list2[basepair-6:basepair-1]+SNPlist2[reference+1]+CHR1list2[basepair:basepair+5]
            newcontext = CHR1list2[basepair-6:basepair-1]+SNPlist2[query+1]+CHR1list2[basepair:basepair+5]
            region = CHR1list2[basepair-500:basepair+500]
        if location2[0] == "II":      
            context = CHR2list2[basepair-6:basepair-1]+SNPlist2[reference+1]+CHR2list2[basepair:basepair+5]
            newcontext = CHR2list2[basepair-6:basepair-1]+SNPlist2[query+1]+CHR2list2[basepair:basepair+5]
            region = CHR2list2[basepair-500:basepair+500]
        if location2[0] == "III":      
            context = CHR3list2[basepair-6:basepair-1]+SNPlist2[reference+1]+CHR3list2[basepair:basepair+5]
            newcontext = CHR3list2[basepair-6:basepair-1]+SNPlist2[query+1]+CHR3list2[basepair:basepair+5]
            region = CHR3list2[basepair-500:basepair+500]
        if location2[0] == "IV":      
            context = CHR4list2[basepair-6:basepair-1]+SNPlist2[reference+1]+CHR4list2[basepair:basepair+5]
            newcontext = CHR4list2[basepair-6:basepair-1]+SNPlist2[query+1]+CHR4list2[basepair:basepair+5]
            region = CHR4list2[basepair-500:basepair+500]
        if location2[0] == "V":      
            context = CHR5list2[basepair-6:basepair-1]+SNPlist2[reference+1]+CHR5list2[basepair:basepair+5]
            newcontext = CHR5list2[basepair-6:basepair-1]+SNPlist2[query+1]+CHR5list2[basepair:basepair+5]
            region = CHR5list2[basepair-500:basepair+500]
        if location2[0] == "X":      
            context = CHRXlist2[basepair-6:basepair-1]+SNPlist2[reference+1]+CHRXlist2[basepair:basepair+5]
            newcontext = CHRXlist2[basepair-6:basepair-1]+SNPlist2[query+1]+CHRXlist2[basepair:basepair+5]
            region = CHRXlist2[basepair-500:basepair+500]
            
#Is it in a restriction site?
        cut1 = 0
        cut2 = 0
        for i in range(6):
            if context[i:i+6] == rsite or context[i:i+6] == rsite3:
                cut1 = 1
            if newcontext[i:i+6] == rsite or newcontext[i:i+6] == rsite3:
                 cut2 = 1
        if cut1 != cut2:
            with open(strain+"_"+strain2+"_"+rsite+".csv", 'a') as csvfile:
                spamwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                spamwriter.writerow([location] + [context] + [newcontext])