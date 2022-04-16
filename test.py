from Bio import Entrez
import random
import string
import pickle
import os
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation
from Bio.Seq import Seq
from Bio import SeqFeature
import os


Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20)) + '@random.com'


handle = Entrez.efetch(db="nucleotide", id="NC_015063", rettype="gb", retmode="text")
print(handle.read())
handle_read = SeqIO.read(handle, "gb")


organism = handle_read.annotations['organism']
features = handle_read.features
regions = ["CDS", "centromere", "intron", "mobile_element", "ncRNA", "rRNA", "telomere", "tRNA", "3'UTR", "5'UTR"]



for f in features:
    if f.type in regions:
        
        seq_result =  f.type + ' ' + organism + ' ' + str(handle_read.id)
        temp_seq = handle_read.seq
        filename = "{}_{}_NC_015063.txt".format(f.type, organism.replace(" ","_"))
        result = open(filename, "a")

        if f.location.strand == -1: # complement detection
            seq_result += ' complement('

            if f.location_operator == "join": # complement(join) detection
                seq_result += ' join({}..{}'.format(f.location.parts[0].start, f.location.parts[0].end)
                for l in f.location.parts[1:]:
                    seq_result += ', {}..{}'.format(l.start, l.end)
                seq_result += ')'
                temp_seq = f.location.extract(handle_read.seq)

            else:
                seq_result += str(f.location.start) + '..' + str(f.location.end)
                fc = SeqFeature.FeatureLocation(f.location.start, f.location.end)
                temp_seq = fc.extract(temp_seq)
            seq_result += ')'

        elif f.location.strand == 1:
            seq_result += ' {}..{}'.format(f.location.start, f.location.end)
            fc = SeqFeature.FeatureLocation(f.location.start, f.location.end)
            temp_seq = fc.extract(temp_seq)

        elif f.location.strand == 0:
            path = os.getcwd()
            os.remove(path + filename)
        
        result.writelines(seq_result+'\n')
        result.writelines(temp_seq+'\n')
        result.close()
