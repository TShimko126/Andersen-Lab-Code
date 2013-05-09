enzyme = raw_input("Which enzyme would you like to use? (Use proper capitalization, no spaces): ")

#Get restriction site from Rebase (maintained by NEB)
import urllib2, re
dbAddress = "http://rebase.neb.com/rebase/enz/"+enzyme+".html"
head = {'User-Agent': 'Mozilla/5.0'}
page = urllib2.Request(dbAddress, headers = head)
page2 = urllib2.urlopen(page)
text = page2.read()
cutSite = re.search("[ATGCRYN]+\^[ATGCRYN]+", text)
cutSite = cutSite.group()
cutSite = re.sub("\^", "", cutSite)

print cutSite


