import os
import csv
from Bio.Seq import *
from Bio.Restriction import *

os.chdir("/Users/tylershimko/Andersen-Lab-Code/snipSNPer/")
outputWriter = csv.writer(open("out.csv", "wb"))
lines = [line.strip() for line in open("out.txt", "r")]
lines = lines[0:-1]
output = {}
for line in lines:
	newLine = line.split("=")
	key = newLine[0]
	value = newLine[1]
	output[key] = value

originalSeq = output["SEQUENCE_TEMPLATE"]

seq = Seq(originalSeq)

primerLeft = output["PRIMER_LEFT_0"]
primerLeft = primerLeft.split(",")
start = int(primerLeft[0])

primerRight = output["PRIMER_RIGHT_0"]
primerRight = primerRight.split(",")
stop = int(primerRight[0]) + int(primerRight[1])

rb = RestrictionBatch(["EcoRI"])
Ana = Analysis(rb, seq, linear=True)
cutsiteNumber = Ana.full()
cutsiteNumber = {EcoRI : [5, 50, 357, 593, 898, 1983]}


refBandSizes = []
for i in range(0,len(cutsiteNumber[EcoRI])):
	print i
	if i == 0:
		refBandSizes.append(cutsiteNumber[EcoRI][i])
	else:
		refBandSizes.append(cutsiteNumber[EcoRI][i] - cutsiteNumber[EcoRI][i-1])

refBandSizes.sort()
print refBandSizes

#ampSeq = Seq(originalSeq[start:stop])

#outputWriter.writerow([enzymeName, chrom, basepair, output["PRIMER_LEFT_0_SEQUENCE"], output["PRIMER_RIGHT_0_SEQUENCE"], cutWorm, cutBands, uncutWorm, uncutBands])