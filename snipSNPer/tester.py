import re


print ""
use = raw_input("Which site numbers would you like to use? (separate with commas):")
use = re.sub(" ", "", use)
use = use.split(",")
for i in range(len(use)):
	use[i] = int(use[i])
print use
print type(use[1])
for site in use:
	print localSites[use]
