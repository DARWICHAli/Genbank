import random
import string
import os
from Bio import SeqIO, Entrez
from Bio.SeqFeature import FeatureLocation
from datetime import datetime, timedelta, date
import time
import re # compiling the pattern for alphanumeric 
import glob
import warnings


class ParserClass:
    @classmethod
    def manage_errors(cls, f, len_seq, signal):
        if f.location.start < 0:
            signal.emit("Error : sequence must not start with 0 or less. {}".format(f.location))
            print("Error : sequence must not start with 0 or less. {}".format(f.location))
            return True
        elif f.location.start > len_seq or f.location.end > len_seq:
            signal.emit("Error : sequence borders must be less or equal than total sequence size. {}".format(f.location))
            print("Error : sequence borders must be less or equal than total sequence size. {}".format(f.location))
            return True
        #elif str(f.location.start)[0] == '<' or str(f.location.end)[0] == '>'\
        #        or str(f.location.start)[0] == '>' or str(f.location.end)[0] == '<':
        try:
            int(str(f.location.start))
            int(str(f.location.end))
        except:
            signal.emit("Error : only numerical arguments are accepted for sequence limits. {}".format(f.location))
            print("Error : only numerical arguments are accepted for sequence limits. {}".format(f.location))
            return True
        return False

    @classmethod
    def parse_NC(cls, NC, path, region_choice, signal):
        Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20)) + '@random.com'
        handle = Entrez.efetch(db="nucleotide", id=NC, rettype="gbwithparts", retmode="text")
        # with open("input.gb", "w") as f:
        #    f.write(handle.read())
        #    f.close()

        with warnings.catch_warnings(record=True) as w:
            handle_read = SeqIO.read(handle, "gb")
            #handle_read = SeqIO.read("test_input.gb", "gb")
            if w:
                for warning in w:
                    signal.emit(warning.message)
                    print(warning.message)


        organism = handle_read.annotations['organism']
        features = handle_read.features

        # delete 'intron' from regions beceause it's not actually a region but just a parsing option
        regions = ["CDS", "centromere", "mobile_element", "ncRNA", "rRNA", "telomere", "tRNA", "3'UTR", "5'UTR"]
        try: 
            file_regions = handle_read.annotations['structured_comment']['Genome-Annotation-Data']['Features Annotated'].split("; ")
        except:
            print("except region")
            file_regions = regions
            
        visited_regions = [False, False, False, False, False, False, False, False, False, False]

        options = region_choice
        selected_regions = []
        intron_is_selected = False
        cds_is_selected = False
        for option in options:
            if option in file_regions:
                if(option == 'CDS'):
                    cds_is_selected = True
                selected_regions.append(option)
            elif option == "intron":
                if "CDS" not in options:
                    selected_regions.append('CDS')
                intron_is_selected = True
            else:
                signal.emit('Pass : {} does not contain selected option \'{}\''.format(NC, option))
                print('Pass : {} does not contain selected option \'{}\''.format(NC, option))

        count_complements = 0
        nb_introns = 0
        for f in features:
            if f.type in selected_regions:
                if f.location:
                    if(cls.manage_errors(f, len(handle_read.seq), signal)):
                        continue

                    intron_seq = ""
                    cds_seq = ""
                    final_seq = ""
                    header = f.type + ' ' + organism + ' ' + str(handle_read.id)

                    ogranism_str = organism.replace(' ', '_').replace('/', '_')
                    ogranism_str = ogranism_str.replace('[', '').replace(']', '').replace(':', '_').replace('\'', '')
                    if f.type == "CDS":
                        if intron_is_selected:
                            intron_filename = path + "/intron_{}_{}.txt".format(ogranism_str, handle_read.id)
                        if cds_is_selected:
                            cds_filename = path + "/CDS_{}_{}.txt".format(ogranism_str, handle_read.id)
                    else:
                        filename = path + "/{}_{}_{}.txt".format(f.type, ogranism_str, handle_read.id)
                    index = regions.index(f.type)
                    if (visited_regions[index] == False):
                        try:
                            if f.type == "CDS":
                                if intron_is_selected:
                                    os.remove(intron_filename)
                                if cds_is_selected:
                                    os.remove(cds_filename)
                            else:
                                os.remove(filename)
                        except:
                            pass

                    if f.type == "CDS":
                        if intron_is_selected:
                            intron_file = open(intron_filename, "a")
                        if cds_is_selected:
                            cds_file = open(cds_filename, "a")
                    else:
                        result = open(filename, "a")

                    visited_regions[index] = True

                    if f.location.strand == -1:  # complement detection
                        count_complements += 1
                        if f.location_operator == "join":  # complement(join) detection
                            if f.type == "CDS":
                                if intron_is_selected:
                                    try:
                                        [intron_seq, nb] = cls.join(header, handle_read.seq, f.location, signal, True)
                                        nb_introns += nb
                                    # gestion erreur
                                    except:
                                        continue

                                if cds_is_selected:
                                    cds_seq = cls.join(header, handle_read.seq, f.location, signal)[0]
                                    # gestion erreur
                                    if not cds_seq:
                                        continue
                            else:
                                final_seq = cls.join(header, handle_read.seq, f.location, signal)[0]
                                # gestion erreur
                                if not final_seq:
                                    continue

                        else:
                            header += ' complement(' + str(f.location.start + 1) + '..' + str(f.location.end)
                            fc = FeatureLocation(f.location.start, f.location.end, -1)
                            header += ')\n'

                            if f.type == "CDS":
                                if cds_is_selected:
                                    cds_seq = header + fc.extract(handle_read.seq)
                            else:
                                final_seq = header + fc.extract(handle_read.seq)

                    elif f.location.strand == 1:
                        count_complements -= 1
                        if f.location_operator == "join":  # join detection
                            if f.type == "CDS":
                                if intron_is_selected:
                                    try:
                                        [intron_seq, nb] = cls.join(header, handle_read.seq, f.location, signal, True)
                                        nb_introns += nb
                                    # gestion erreur
                                        if(not intron_seq): continue
                                    except: continue
                                if cds_is_selected:
                                    cds_seq = cls.join(header, handle_read.seq, f.location, signal)[0]
                                    # gestion erreur
                                    if not cds_seq: continue
                            else:
                                final_seq = cls.join(header, handle_read.seq, f.location, signal)[0]
                                # gestion erreur
                                if not final_seq: continue

                        # final_seq = header + '\n' + f.extract(handle_read.seq)
                        else:
                            header += ' {}..{}'.format(f.location.start + 1, f.location.end)
                            fc = FeatureLocation(f.location.start, f.location.end, 1)
                            if f.type == "CDS":
                                #if intron_is_selected:
                                #    intron_seq = header + '\n' + fc.extract(handle_read.seq)
                                if cds_is_selected:
                                    cds_seq = header + '\n' + fc.extract(handle_read.seq)
                            else:
                                final_seq = header + '\n' + fc.extract(handle_read.seq)

                    elif f.location.strand == 0:
                        signal.emit('Error : noisy strand')
                        print('Error : noisy strand')
                        continue

                    else:
                        signal.emit('Error : cannot join complementory and normal strands')
                        print('Error : cannot join complementory and normal strands')

                    if f.type == "CDS":
                        if intron_is_selected:
                            if intron_seq:
                                intron_file.writelines(intron_seq + '\n')
                            intron_file.close()
                        if cds_is_selected:
                            if cds_seq:
                                cds_file.writelines(cds_seq + '\n')
                            cds_file.close()
                    elif final_seq:
                        result.writelines(final_seq + '\n')
                        result.close()
        if(nb_introns == 0 and intron_is_selected):
            intron_file.close()
            os.remove(intron_file.name)
        #print("number of introns found: {}".format(nb_introns))
        return True

    @classmethod
    def join(cls, header, sequence, location, signal, intron=False):
        final_seq = ""
        isComplement = location.strand
        header += ' complement(' if (isComplement == -1) else " "
        index = 0 if (isComplement == 1) else -1

        first_feature = location.parts[index]
        header += 'join({}..{}'.format(first_feature.start + 1, first_feature.end)
        last_end = first_feature.end
        if intron:
            seq_join_intern = []
        if not intron:
            seq_join_intern = [str(first_feature.extract(sequence))]

        nb_introns = 0
        for l in (location.parts[::isComplement])[1:]:
            if (l.start + 1) < last_end:
                signal.emit("Error : invalid join sequence order")
                print("Error : invalid join sequence order")
                return (None,0)
            header += ',{}..{}'.format(l.start + 1, l.end)
            if intron:
                f = FeatureLocation(last_end, l.start+1, strand=isComplement)
                seq_join_intern.append(str(f.extract(sequence)))
                nb_introns += 1
            else:
                seq_join_intern.append(str(l.extract(sequence)))
            last_end = l.end
            
        header += '))' if (isComplement == -1) else ')'

        seq_join_intern = seq_join_intern[::isComplement]

        final_seq += header + '\n' + "".join(seq_join_intern)

        header = header.replace("CDS", "intron")
        for (i, s) in enumerate(seq_join_intern):
            key_word = ' Intron ' if intron else ' Exon '
            final_seq += '\n' + header + key_word + str(i + 1) + '\n' + str(s)

        return (final_seq, nb_introns)

    @classmethod
    def verify_modification_date(path, handle_date):
        file_date = os.path.getmtime(path)
        file_date = time.strftime('%Y%m%d', time.localtime(file_date))
        handle_date = datetime.strptime(handle_date, "%d-%b-%Y")
        handle_date =  str(handle_date).split(' ')[0].replace('-','')
        # fichier déjà parsé et non modifié:
        if int(file_date) >= int(handle_date):
            return True
        return False