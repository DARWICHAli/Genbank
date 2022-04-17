import urllib.request
from Bio import Entrez
import random
import string
import re

#filename_expr = re.compile(r"join().txt")

page = urllib.request.urlopen('https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&retmode=raw&withparts=on&basic_feat=on&id=NC_020356')

data = page.read()
data = data.decode("utf-8") 
data = data.split("\n")

s = 0
e = 0
founds = False
for i in data:
    if not founds:
        s+=1
    e+=1
    if i.startswith("FEATURES"):
        founds = True
    if i.startswith("ORIGIN"):
        e-=2
        break

print(data[s])
print(data[e])