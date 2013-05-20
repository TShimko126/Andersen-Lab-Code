#Get restriction site from Rebase (maintained by NEB)
import urllib2, re, matplotlib, string

chrom = raw_input("Enter your chromosome of interest (use numbers, not roman numerals): ")
chrom = chrom.upper()
print chrom

pos = raw_input("Enter your position/region of interest (separate regional endpoints with ..): ")
pos = pos.split("..")
posLen = len(pos)
print pos

enz = raw_input("Which enzymes would you like to use? (Use proper capitalization, separate with commas): ")
enz = re.sub(" ","",enz)
enz = enz.split(",")

print enz

cutsites = {}

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

def nearestCut(enzyme, position):
    