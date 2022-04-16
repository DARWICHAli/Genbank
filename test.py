import urllib.request
from Bio import Entrez
import random
import string

# page = urllib.request.urlopen('https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&retmode=raw&withparts=on&basic_feat=on&id=NC_020356')

# data = page.read()
# data = data.decode("utf-8") 
# data = data.split("\n")
# print(data)



Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20)) + '@random.com'
handle = Entrez.efetch(db="nucleotide", id="NC_015060", rettype="gb", retmode="text")
print(type(handle.read()))