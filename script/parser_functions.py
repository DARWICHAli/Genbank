import time
import numpy as np
import os.path
from ftp_downloader import *
import random
import string
import os
from Bio import SeqIO, Entrez
from Bio.SeqFeature import FeatureLocation, CompoundLocation, ExactPosition, SeqFeature
from datetime import datetime
import time
import warnings
import re

green = [0,255,0,255]
white = [255,255,255,255]
purple = [255,0,255,255]

stop = False
count = 0
bdd_path = "last_opened_paths.txt"
class ParserFunctions:
    verbose = True
    # This is a generator function!
    def __init__(self, regions_choice, mutex, mutex_fetch, queue, stop_event):
        self.region_choice = regions_choice
        self.mutex = mutex
        self.mutex_fetch = mutex_fetch
        self.queue = queue
        self.stop_event = stop_event

    def write_bdd(self, path, file, mutex):
        mutex.acquire()
        with open(path, "a") as bdd:
            bdd.writelines(file+'\n')
        mutex.release()

    def copy_bdd(self, path, filenames, mutex):

        mutex.acquire()
        with open(path, "r") as bdd:
            lines = bdd.readlines()
            lines = lines.copy()
        new_bdd = open(path, "w")
        for l in lines:
            if l.strip('\n') not in filenames:
                new_bdd.writelines(l)
        new_bdd.close()
        mutex.release()
    def emit_log(self, log_signal, msg, color):
        if self.verbose:
            self.mutex.acquire()
            log_signal.put({'msg':msg, 'color':color, 'type':1}, block=True)
            self.mutex.release()



    def manage_errors(self, f, len_seq, log_signal):
        if f.location.start < 0:
            self.emit_log(log_signal, "Pass: sequence must not start with 0 or less. {}".format(f.location), purple)
            return True
        elif f.location.start > len_seq or f.location.end > len_seq:
            self.emit_log(log_signal, "Pass: sequence borders must be less or equal than total sequence size. {}".format(f.location), purple)
            return True

        try:
            int(str(f.location.start))
            int(str(f.location.end))
        except:
            self.emit_log(log_signal, "Pass: only numerical arguments are accepted for sequence limits. {}".format(f.location), purple)
            return True
        return False

    def parse_NC(self, row_df):#, region_choice, log_signal, progress_signal, organism_df, mutex, mutex_fetch, mutex_count, mutex_stop):
        global count

        mutex = self.mutex
        mutex_fetch = self.mutex_fetch
        i = 0
        global stop
        index_df, organism_name, path, NC, file_features = row_df
        #print(f"### Organism: {organism_name}")
        cp_region_choice = list(self.region_choice)
        #for NC in NC_LIST:
        if(self.stop_event.is_set()):
            self.emit_log(self.queue, f'Pass {NC}: Mutex Stop', purple)
            #print(f'######################## \n \n Pass {NC}: Mutex Stop')
            return
        i+=1
        msg = f"Parsing {str(NC)} in: {str(path)}"
        self.emit_log(self.queue, msg, white)
        Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20)) + '@random.com'
        Entrez.api_key = 'd190df79af669bbea636278753c3236afa08'

        # Premier efetch pour récupérer juste la date.
        # le read prends prend plus de temps que la partie parsing, donc il faut vérifier la date avant
        t = time.time()
        # gestion timeout
        nb_tries = 10
        Entrez.max_tries = nb_tries
        try:
            #print(f"Try Number : {i}, NC: {NC}")
            #mutex_fetch.acquire()
            handle = Entrez.efetch(db="nucleotide", id=NC, rettype="gb", retmode="text", datetype='mdat')
            #mutex_fetch.release()

        except Exception as e:
            self.emit_log(self.queue, f'Pass {NC}: Efetch error {e}', purple)
            return

        print(f" ## Efetch 1: {time.time() - t}")

        handle_date = handle.readline().strip().split(' ')
        handle_date = handle_date[len(handle_date)-1]
        NC_modified_init = self.verify_modification_date(path, handle_date)
        NC_modified = NC_modified_init
        new_region_choice = []

        for option in cp_region_choice:
            if(len(file_features) and option not in file_features):
                self.emit_log(self.queue, "Skip: option " + option + " in " + NC + " is not found.", purple)
                continue
            NC_modified = NC_modified_init

            filename = "{}_{}_{}".format(option, organism_name, NC)
            pattern = re.compile(filename + ".\d+.txt")
            file_exists = False

            for filepath in os.listdir(path):
                if pattern.match(filepath):
                    file_exists = True
                    break

            NC_modified = NC_modified or not file_exists

            if(NC_modified):
                new_region_choice.append(option)
            else:
                self.emit_log(self.queue, "Skip: option " + option + " in " + NC + " is up to date.", green)

        if(len(new_region_choice) == 0):
            self.emit_log(self.queue, f"Pass {NC}: No regions to parse", purple)
            self.queue.put({'type':2})
            return

        cp_region_choice = new_region_choice

        Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20)) + '@random.com'
        Entrez.api_key = '24f75a54c3b96e27fcae236d031f3cfb4909'

        t = time.time()
        Entrez.max_tries = nb_tries
        # deuxieme efectch car on perd les informations du premiers
        try:
            #mutex_fetch.acquire()
            handle = Entrez.efetch(db="nucleotide", id=NC, rettype="gbwithparts", retmode="text", datetype='mdat')
            #mutex_fetch.release()
        except Exception as e:
            self.emit_log(self.queue, f'Pass {NC}: Efetch error {e}', purple)
            return

        print(f" ## Efetch 2: {time.time() - t}")

        t = time.time()
        nb_tries = 2
        with warnings.catch_warnings(record=True) as w:
            for i in range(nb_tries):
                try:
                    #handle_read = SeqIO.read("./saved.gb", "gb")
                    handle_read = SeqIO.read(handle, "gb")
                    break
                except Exception as e:
                    t_start = time.time()
                    while(True):
                        if self.stop_event.is_set():
                            self.emit_log(self.queue, f'Pass {NC}: Mutex Stop', purple)
                            return
                        if(time.time()-t_start > nb_tries):
                            print(f'Pass {NC}: Error Efetch Read {e}\nTry Again ... {i}/{nb_tries}')
                            break
            #SeqIO.write(handle_read, "saved.gb", "gb")

            if (i == (nb_tries - 1)):
                self.emit_log(self.queue, f'{time.ctime(time.time())}: Pass {NC}: Error Handle Read.', purple)
                return
            if w:
                for warning in w:
                    self.emit_log(self.queue, warning.message, purple)

        print(f" ## Handle read: {time.time() - t}")

        organism = handle_read.annotations['organism']
        features = handle_read.features

        # delete 'intron' from regions beceause it's not actually a region but just a parsing option
        regions = ["CDS", "centromere", "mobile_element", "ncRNA", "rRNA", "telomere", "tRNA", "3'UTR", "5'UTR"]
        try:
            file_regions = handle_read.annotations['structured_comment']['Genome-Annotation-Data']['Features Annotated'].split("; ")
        except:
            file_regions = regions

        selected_regions = []
        intron_is_selected = False
        cds_is_selected = False
        intron_filename = ""
        cds_filename = ""
        filename = ""

        # we update the dataframe
        mutex.acquire()
        self.queue.put({'type':5, 'col':'features', 'row':index_df, 'value':list(set([f for f in file_regions]))})
        mutex.release()
        ### organism_df['features'][index_df] = [f for f in file_regions]
        # if 0 introns found, this will be deleted from the dataframe, so that we don't verify it again in
        # in the next executions
        ###if("CDS" in file_regions):
        ###    organism_df['features'][index_df].append("intron")
        if("CDS" in file_regions):
            mutex.acquire()
            self.queue.put({'type':4, 'col':'features', 'row':index_df})
            mutex.release()

        for option in cp_region_choice:
            if option in file_regions:
                if(option == 'CDS'):
                    cds_is_selected = True
                selected_regions.append(option)
            elif option == "intron":
                if "CDS" not in cp_region_choice:
                    selected_regions.append('CDS')
                intron_is_selected = True
            else:
                self.emit_log(self.queue,f'Pass : {NC} does not contain selected option \'{option}\'', purple)

        if(not len(selected_regions)):
            self.emit_log(self.queue,'Pass : {} does not contain any of the selected options \'{}\''.format(NC, option), purple)
            #print(f'######################## \n \n Pass : {NC} does not contain selected option \'{option}\'')
            return

        visited_regions = [False for i in selected_regions]

        count_complements = 0
        nb_introns = 0
        str_features = list(map(lambda f: str(f.location) + ',' + str(f.type), features))
        _ , unique_indexes = np.unique(np.array(str_features), return_index=True)
        unique_features = np.array(features)[unique_indexes]


        for f in unique_features:
            if(self.stop_event.is_set()):
                self.emit_log(self.queue, f'Pass {NC}: Mutex Stop', purple)
                return
            if f.type in selected_regions:
                index = selected_regions.index(f.type)

                if f.location:
                    if(self.manage_errors(f, len(handle_read.seq), self.queue)):
                        continue

                    intron_seq = ""
                    cds_seq = ""
                    final_seq = ""
                    header = f.type + ' ' + organism + ' ' + str(handle_read.id)+':'

                    if f.type == "CDS":
                        if intron_is_selected:
                            intron_filename = path + "intron_{}_{}.txt".format(organism_name, handle_read.id)
                        if cds_is_selected:
                            cds_filename = path + "CDS_{}_{}.txt".format(organism_name, handle_read.id)
                    else:
                        filename = path + "{}_{}_{}.txt".format(f.type, organism_name, handle_read.id)

                    index = selected_regions.index(f.type)
                    if (visited_regions[index] == False):
                        try:
                            if f.type == "CDS":
                                if intron_is_selected:
                                    if os.path.isfile(intron_filename):
                                        os.remove(intron_filename)
                                if cds_is_selected:
                                    if os.path.isfile(cds_filename):
                                        os.remove(cds_filename)
                            else:
                                if os.path.isfile(filename):
                                    os.remove(filename)
                        except Exception as e:
                            print(f"Error Remove Filename: {e}")
                            pass


                    if f.type == "CDS":
                        if intron_is_selected:
                            if(visited_regions[index] == False):
                                self.write_bdd(bdd_path, intron_filename, mutex)
                            intron_file = open(intron_filename, "a")
                        if cds_is_selected:
                            if(visited_regions[index] == False):
                                self.write_bdd(bdd_path, cds_filename, mutex)
                            cds_file = open(cds_filename, "a")
                    else:
                        if(visited_regions[index] == False):
                            self.write_bdd(bdd_path, filename, mutex)
                        result = open(filename, "a")

                    visited_regions[index] = True

                    if f.location.strand == -1:  # complement detection
                        count_complements += 1
                        if f.location_operator == "join":  # complement(join) detection
                            if f.type == "CDS":
                                if intron_is_selected:
                                    try:
                                        (intron_seq, nb) = self.join(header, handle_read.seq, f.location, self.queue, True)
                                        nb_introns += nb
                                    # gestion erreur
                                    except Exception as e:
                                        print(f"Error Intron strand -1 {e}")
                                        continue

                                if cds_is_selected:
                                    cds_seq = self.join(header, handle_read.seq, f.location, self.queue)[0]
                                    # gestion erreur
                                    if not cds_seq:
                                        continue
                            else:
                                final_seq = self.join(header, handle_read.seq, f.location, self.queue)[0]
                                # gestion erreur
                                if not final_seq:
                                    continue

                        else:
                            header += f' complement({str(f.location.start + 1)}..{str(f.location.end)}'
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
                                        (intron_seq, nb) = self.join(header, handle_read.seq, f.location, self.queue, True)
                                        nb_introns += nb
                                    # gestion erreur
                                        if(not intron_seq): continue
                                    except Exception as e:
                                        print(f"Intron Error strand 1: {e}")
                                        continue
                                if cds_is_selected:
                                    cds_seq = self.join(header, handle_read.seq, f.location, self.queue)[0]
                                    # gestion erreur
                                    if not cds_seq: continue
                            else:
                                final_seq = self.join(header, handle_read.seq, f.location, self.queue)[0]
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
                        self.emit_log(self.queue,'Pass ' + NC + ': noisy strand', purple)
                        continue

                    else:
                        self.emit_log(self.queue,'Pass ' + NC + ': cannot join complementory and normal strands', purple)

                    file_names = []

                    if f.type == "CDS":
                        if intron_is_selected:
                            if intron_seq:
                                intron_file.writelines(intron_seq + '\n')
                            intron_file.close()
                            file_names.append(intron_filename)

                        if cds_is_selected:
                            if cds_seq:
                                cds_file.writelines(cds_seq + '\n')
                            cds_file.close()
                            file_names.append(cds_filename)

                    elif final_seq:
                        result.writelines(final_seq + '\n')
                        result.close()
                        file_names.append(filename)

                    self.copy_bdd(bdd_path, file_names, mutex)

        if(nb_introns == 0 and intron_is_selected):
            # we don't want to parse intron from this NC ever again!!!! (unless you redownload the reports)
            mutex.acquire()
            self.queue.put({'type':3, 'col':'features', 'row':index_df})
            mutex.release()
            ###if("intron" in organism_df['features'][index_df]):
            ###    organism_df['features'][index_df].remove("intron")
            try:
                intron_file.close()
            except Exception as e:
                #print(e)
                pass

            try:
                os.remove(intron_filename)
            except Exception as e:
                #print(f"Error Remove Intron File : {e}")
                pass
            self.emit_log(self.queue, f"No Intron found in {NC}.", purple)

        self.emit_log(self.queue, f"Parsing of {NC} done successfully.", green)
        mutex.acquire()
        self.queue.put({'type':2})
        mutex.release()
        return True




    def join(self, header, sequence, location, log_signal, intron=False):
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
            if (l.start + 1) <= last_end:
                self.emit_log(log_signal, "Pass: invalid join sequence order", purple)
                return (None,0)
            header += ',{}..{}'.format(l.start + 1, l.end)
            if intron:
                f = FeatureLocation(last_end, l.start, strand=isComplement)
                seq_join_intern.append(str(f.extract(sequence)))
                nb_introns += 1
            else:
                seq_join_intern.append(str(l.extract(sequence)))
            last_end = l.end

        header += '))' if (isComplement == -1) else ')'

        seq_join_intern = seq_join_intern[::isComplement]

        sep = ''
        if (len(seq_join_intern) > 1):
            if(not intron):
                final_seq += header + '\n' + "".join(seq_join_intern)
            sep = '\n'

        if intron:
            header = header.replace("CDS", "intron")
        for (i, s) in enumerate(seq_join_intern):
            key_word = ' Intron ' if intron else ' Exon '
            if not intron :
                final_seq += sep + header + key_word + str(i + 1) + '\n' + str(s)
            else :
                final_seq += header + key_word + str(i + 1) + '\n' + str(s)
                final_seq += sep if i < len(seq_join_intern)-1 else ""

        return (final_seq, nb_introns)



    def verify_modification_date(self, path, handle_date):
        if(not os.path.isdir(path)):
            return True
        file_date = os.path.getmtime(path)
        file_date = time.strftime('%Y%m%d', time.localtime(file_date))
        handle_date = datetime.strptime(handle_date, "%d-%b-%Y")
        handle_date =  str(handle_date).split(' ')[0].replace('-','')
        # fichier déjà parsé et non modifié:
        if int(file_date) >= int(handle_date):
            return False
        return True
