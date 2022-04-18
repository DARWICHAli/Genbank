from Bio import Entrez
import random
import string
import os
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

import os
import glob
import re # compiling the pattern for alphanumeric 

class Parser:


    @classmethod
    def parse_NC(cls, NC, path, region_choice):
        
        Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20)) + '@random.com'

        handle = Entrez.efetch(db="nucleotide", id=NC, rettype="gbwithparts", retmode="text")
        handle_read = SeqIO.read(handle, "gb")

        organism = handle_read.annotations['organism']
        features = handle_read.features
        #regions = ["CDS", "centromere", "intron", "mobile_element", "ncRNA", "rRNA", "telomere", "tRNA", "3'UTR", "5'UTR"]
        visited_regions = [False for i in range(len(region_choice))]

        bornes_expr = re.compile(r"[1-9][0-9]*")

        for f in features:
            # si cette région a été selectionnée par l'utilisateut
            if f.type in region_choice:

                if cls.error_check(f.location.parts, bornes_expr, handle_read):
                    return False

                header =  f.type + ' ' + organism + ' ' + str(handle_read.id)
                final_seq = ""

                filename = "{}/{}_{}_{}.txt".format(path ,f.type, organism.replace(" ","_"), handle_read.id)

                index = region_choice.index(f.type)
                if(visited_regions[index] == False):
                    try: os.remove(filename)
                    except: pass
                
                result = open(filename, "a")
                visited_regions[index] = True
                
                # complement detection
                if f.location.strand == -1: 
                    
                    if f.location_operator == "join": # complement(join) detection
                        final_seq = cls.join(header, handle_read.seq, f.location)

                    elif f.location_operator is None: # complément sans join
                        header += ' complement('+ str(f.location.start+1) + '..' + str(f.location.end)
                        fc = FeatureLocation(f.location.start, f.location.end+1)
                        header += ')\n'
                        final_seq = header + fc.extract(handle_read.seq)

                    # todo: supprimer les fichiers créés
                    else:
                        print("error")
                        break
                
                # pas de complémént
                elif f.location.strand == 1:

                    if f.location_operator == "join": # complement(join) detection
                        final_seq = cls.join(header, handle_read.seq, f.location)

                    elif f.location_operator is None:
                        header += ' {}..{}'.format(f.location.start+1, f.location.end)
                        fc = FeatureLocation(f.location.start, f.location.end+1)
                        final_seq = header + '\n' + fc.extract(handle_read.seq)

                    # todo: supprimer les fichiers créés
                    else:
                        print("error")
                        break
                
                # mauvais strand
                elif f.location.strand == 0:
                    os.remove(filename)

                # résultat final pour cette région
                result.writelines(final_seq+'\n')
                result.close()
        return True

    @classmethod
    def join(cls, header, sequence, location):
        final_seq = ""
        isComplement = location.strand
        header += ' complement(' if(isComplement == -1) else ""
        index = 0 if(isComplement == 1) else -1
        
        a = FeatureLocation(location.parts[index].start, location.parts[index].end+1, strand = isComplement)

        header += ' ' if(isComplement == 1) else ""
        header += 'join({}..{}'.format(a.start+1, a.end-1)
        seq_join_interm = [str(a.extract(sequence))]

        for l in (location.parts[::isComplement])[1:]:
            b = FeatureLocation(l.start, l.end+1, strand = isComplement)
            header += ',{}..{}'.format(b.start+1, b.end-1)
            seq_join_interm.append(str(b.extract(sequence)))

        header += '))' if(isComplement == -1) else ')' 

        seq_join_interm = seq_join_interm[::isComplement]

        final_seq += header + '\n' + "".join(seq_join_interm)

        for (i,s) in enumerate(seq_join_interm):
            final_seq += '\n' + header + ' Exon ' + str(i+1) + '\n' + str(s)

        return final_seq


    @classmethod
    def error_check(cls, location_parts, bornes_expr, handle_read):
        for l in location_parts:
            if l.start > l.end or re.fullmatch(bornes_expr, str(l.start)) is None or re.fullmatch(bornes_expr, str(l.end)) is None:
                filespath = glob.glob("*_{}.txt".format(handle_read.id))
                for file in filespath:
                    try: os.remove(file)
                    except: print("cant remove file")
                return True
        return False
