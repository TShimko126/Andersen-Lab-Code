import os
from rpy2.robjects import r as r
cwd = os.getcwd()
r.assign("cwd", cwd)
save = raw_input("Would you like to save the chromosome map? (yes/no): ")
if save == "yes" or save == "y":
	filename = []
	filename.append(raw_input("What would you like to name the file (do not include filetype extension): "))
	filename.append(".png")
	filename = "".join(filename)
else:
	filename = "standinfilename.png"
cwd = os.getcwd()
r.assign("cwd", cwd)
r.assign("filename", filename)
r("source('RPlottingScript.R')")
