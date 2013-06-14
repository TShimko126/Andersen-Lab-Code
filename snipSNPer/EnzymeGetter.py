def getter():
	from BeautifulSoup import BeautifulSoup
	import urllib2, re

	address = "https://www.neb.com/tools-and-resources/selection-charts/time-saver-qualified-restriction-enzymes"
	site = urllib2.urlopen(address)
	soup = BeautifulSoup(site)

	tables = soup.findAll("table")
	table = tables[1]
	# print table

	rows = table.findAll("tr")
	rows = rows[1:]

	enzlist = []
	for row in rows:
		cols = row.findAll("td")
		enzcut = []
		for ind in range(2):
			match = re.search("<a", str(cols[ind].contents))
			if match:
				links = cols[ind].findAll("a")
				contents = str(links[0].contents)
				contents = re.sub("\[u\'", "", contents)
				contents = re.sub("-HF&reg;", "", contents)
				contents = re.sub(" RE-Mix&reg;", "", contents)
				contents = re.sub("&nbsp;", "", contents)
				contents = re.sub("\', <sup>&alpha;</sup>, u\'", "", contents)
				contents = re.sub("-HF&trade;", "", contents)
				contents = re.sub("\']", "", contents)
				enzcut.append(contents)
			else:
				contents = str(cols[ind].contents)
				contents = re.sub("\[u\'", "", contents)
				contents = re.sub("/", "", contents)
				contents = re.sub("\']", "", contents)
				enzcut.append(contents)
		enzlist.append(enzcut)

	indices = []
	for i in range(len(enzlist)):
		match1 = re.search("\(-?[0-9]+(/?-?[0-9]+)?\)", enzlist[i][1])
		match2 = re.search("[RYWSKMBDHVN]", enzlist[i][1])
		if len(enzlist[i][1]) > 4:
			if match1:
				continue
			else:	
				if match2:
					continue
				else:
					ind = enzlist.index(enzlist[i])
					indices.append(ind)
	
	indices = list(set(indices))



	finalEnzList = []
	for i in range(len(indices)):
		finalEnzList.append(enzlist[indices[i]])
	return finalEnzList
