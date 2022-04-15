import urllib.request

page = urllib.request.urlopen('https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&retmode=raw&withparts=on&basic_feat=on&id=NC_020356')

data = page.read()
data = data.decode("utf-8") 
data = data.split("\n")
print(data)