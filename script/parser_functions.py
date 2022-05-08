from PyQt5 import QtCore
import time
import pickle
from nbformat import write
import pandas as pd
import os.path
from ftp_downloader import *
import random
import string
import os
from Bio import SeqIO, Entrez, AlignIO
from Bio.SeqFeature import FeatureLocation
from datetime import datetime
import time
import warnings
from Bio.SeqIO import AbiIO
from Bio.SeqIO import AceIO
from Bio.SeqIO import FastaIO
from Bio.SeqIO import GckIO
from Bio.SeqIO import IgIO  # IntelliGenetics or MASE format
from Bio.SeqIO import InsdcIO  # EMBL and GenBank
from Bio.SeqIO import NibIO
from Bio.SeqIO import PdbIO
from Bio.SeqIO import PhdIO
from Bio.SeqIO import PirIO
from Bio.SeqIO import QualityIO  # FastQ and qual files
from Bio.SeqIO import SeqXmlIO
from Bio.SeqIO import SffIO
from Bio.SeqIO import SnapGeneIO
from Bio.SeqIO import SwissIO
from Bio.SeqIO import TabIO
from Bio.SeqIO import TwoBitIO
from Bio.SeqIO import UniprotIO
from Bio.SeqIO import XdnaIO

from threading import Thread, Lock

green = [0,255,0,255]
white = [255,255,255,255]
purple = [255,0,255,255] 
red = [255,0,0,255] 

stop = False
count = 0
bdd_path = "last_opened_paths.txt"
_FormatToIterator = {
    "abi": AbiIO.AbiIterator,
    "abi-trim": AbiIO._AbiTrimIterator,
    "ace": AceIO.AceIterator,
    "fasta": FastaIO.FastaIterator,
    "fasta-2line": FastaIO.FastaTwoLineIterator,
    "ig": IgIO.IgIterator,
    "embl": InsdcIO.EmblIterator,
    "embl-cds": InsdcIO.EmblCdsFeatureIterator,
    "gb": InsdcIO.GenBankIterator,
    "gck": GckIO.GckIterator,
    "genbank": InsdcIO.GenBankIterator,
    "genbank-cds": InsdcIO.GenBankCdsFeatureIterator,
    "imgt": InsdcIO.ImgtIterator,
    "nib": NibIO.NibIterator,
    "cif-seqres": PdbIO.CifSeqresIterator,
    "cif-atom": PdbIO.CifAtomIterator,
    "pdb-atom": PdbIO.PdbAtomIterator,
    "pdb-seqres": PdbIO.PdbSeqresIterator,
    "phd": PhdIO.PhdIterator,
    "pir": PirIO.PirIterator,
    "fastq": QualityIO.FastqPhredIterator,
    "fastq-sanger": QualityIO.FastqPhredIterator,
    "fastq-solexa": QualityIO.FastqSolexaIterator,
    "fastq-illumina": QualityIO.FastqIlluminaIterator,
    "qual": QualityIO.QualPhredIterator,
    "seqxml": SeqXmlIO.SeqXmlIterator,
    "sff": SffIO.SffIterator,
    "snapgene": SnapGeneIO.SnapGeneIterator,
    "sff-trim": SffIO._SffTrimIterator,  # Not sure about this in the long run
    "swiss": SwissIO.SwissIterator,
    "tab": TabIO.TabIterator,
    "twobit": TwoBitIO.TwoBitIterator,
    "uniprot-xml": UniprotIO.UniprotIterator,
    "xdna": XdnaIO.XdnaIterator,
}

class ParserFunctions:
    def parse(cls, handle, format, mutex, alphabet=None):
        global stop
        # Try and give helpful error messages:
        if not isinstance(format, str):
            raise TypeError("Need a string for the file format (lower case)")
        if not format:
            raise ValueError("Format required (lower case string)")
        if not format.islower():
            raise ValueError("Format string '%s' should be lower case" % format)
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")

        iterator_generator = _FormatToIterator.get(format)
        if iterator_generator:
            return iterator_generator(handle)
        if format in AlignIO._FormatToIterator:
            # Use Bio.AlignIO to read in the alignments
            for alignment in AlignIO.parse(handle, format):
                for r in alignment:
                    if cls.check_stopping(mutex):
                        return None
                    return r
        raise ValueError("Unknown format '%s'" % format)


    def read(cls, handle, format, mutex, alphabet=None):
        iterator = cls.parse(handle, format, mutex, alphabet)
        if iterator is None:
            return None
        try:
            record = next(iterator)
        except StopIteration:
            raise ValueError("No records found in handle") from None
        try:
            next(iterator)
            raise ValueError("More than one record found in handle")
        except StopIteration:
            pass
        return record

    def reinit(self):
        global stop
        global count
        stop = 0
        count = 0

    def set_stop(self, mutex):
        global stop
        mutex.acquire()
        stop = True
        mutex.release()
        
    def get_count(self, mutex):
        global count
        mutex.acquire()
        true_count = count
        mutex.release()
        return true_count

    def write_bdd(self, path, file, mutex):
        mutex.acquire()
        bdd = open(path, "a")
        bdd.writelines(file+'\n')
        bdd.close()
        mutex.release()

    def copy_bdd(self, path, filenames, mutex):

        mutex.acquire()
        bdd = open(path, "r")
        lines = bdd.readlines()
        lines = lines.copy()
        bdd.close()
        new_bdd = open(path, "w")
        for l in lines:
            if l.strip('\n') not in filenames:
                new_bdd.writelines(l)
        new_bdd.close()
        mutex.release()

    def check_stopping(self, mutex):
        mutex.acquire()
        stopping = stop
        mutex.release()
        return stopping


    def manage_errors(self, f, len_seq, log_signal):
        if f.location.start < 0:
            log_signal.emit("Pass: sequence must not start with 0 or less. {}".format(f.location), purple)
            ##print("Pass: sequence must not start with 0 or less. {}".format(f.location))
            return True
        elif f.location.start > len_seq or f.location.end > len_seq:
            log_signal.emit("Pass: sequence borders must be less or equal than total sequence size. {}".format(f.location), purple)
            ##print("Pass: sequence borders must be less or equal than total sequence size. {}".format(f.location))
            return True
        #elif str(f.location.start)[0] == '<' or str(f.location.end)[0] == '>'\
        #        or str(f.location.start)[0] == '>' or str(f.location.end)[0] == '<':
        try:
            int(str(f.location.start))
            int(str(f.location.end))
        except:
            log_signal.emit("Pass: only numerical arguments are accepted for sequence limits. {}".format(f.location), purple)
            ##print("Pass: only numerical arguments are accepted for sequence limits. {}".format(f.location))
            return True
        return False

    def parse_NC(self, row_df, region_choice, log_signal, progress_signal, organism_df, mutex, mutex_fetch, mutex_count, mutex_stop):
        global count
        mutex_count.acquire()
        count = count + 1
        mutex_count.release()
        i = 0
        global stop
        index_df, organism_name, path, NC_LIST, file_features = row_df
        cp_region_choice = [i for i in region_choice]
        for NC in NC_LIST:
            if(self.check_stopping(mutex_stop)):
                mutex_count.acquire()
                count = count - 1
                mutex_count.release()
                return
            i+=1
            msg = "Parsing " + str(NC) + ' in: ' + str(path)
            log_signal.emit(msg, white)
            ##print(msg)
            Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20)) + '@random.com'
            Entrez.api_key = 'd190df79af669bbea636278753c3236afa08'

            # Premier efetch pour récupérer juste la date.
            # le read prends prend plus de temps que la partie parsing, donc il faut vérifier la date avant


            try:
                handle = Entrez.efetch(db="nucleotide", id=NC, rettype="gb", retmode="text", datetype='mdat')
            except:
                print ("Error fetching")
                mutex_fetch.acquire()
                timer = 0.0
                while(timer != 5.0):
                    if self.check_stopping(mutex_stop):
                        return
                    timer += 0.1
                handle = Entrez.efetch(db="nucleotide", id=NC, rettype="gb", retmode="text", datetype='mdat')
                mutex_fetch.release()

            handle_date = handle.readline().strip().split(' ')
            handle_date = handle_date[len(handle_date)-1]
            NC_modified_init = self.verify_modification_date(path, handle_date)
            NC_modified = NC_modified_init
            new_region_choice = []


            for option in cp_region_choice:
                if(len(file_features) and option not in file_features):
                    continue
                NC_modified = NC_modified_init
                filename = path + "{}_{}_{}.{}.txt".format(option, organism_name, NC, i )
                NC_modified = NC_modified or not os.path.isfile(filename)
                if(NC_modified):
                    new_region_choice.append(option)
                else:
                    log_signal.emit("Skip: option " + option + " in " + NC + " is up to date.", purple)

            if(not len(new_region_choice)):
                continue

            cp_region_choice = new_region_choice

            Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20)) + '@random.com'
            Entrez.api_key = '24f75a54c3b96e27fcae236d031f3cfb4909'
            # deuxieme efectch car on perd les informations du premiers
            try:
                handle = Entrez.efetch(db="nucleotide", id=NC, rettype="gbwithparts", retmode="text", datetype='mdat')
            except:
                mutex_fetch.acquire()
                timer = 0.0
                while(timer != 5.0):
                    if self.check_stopping(mutex_stop):
                        return
                    timer += 0.1
                handle = Entrez.efetch(db="nucleotide", id=NC, rettype="gbwithparts", retmode="text", datetype='mdat')
                mutex_fetch.release()

            with warnings.catch_warnings(record=True) as w:
                handle_read = self.read(handle, "gb", mutex_stop)
                #handle_read = SeqIO.read("test_input.gb", "gb")
                if handle_read is None:
                    mutex_count.acquire()
                    count = count - 1
                    mutex_count.release()
                    return
                if w:
                    for warning in w:
                        log_signal.emit(warning.message, purple)
                        #print(warning.message)
                        

            organism = handle_read.annotations['organism']
            features = handle_read.features

            # delete 'intron' from regions beceause it's not actually a region but just a parsing option
            regions = ["CDS", "centromere", "mobile_element", "ncRNA", "rRNA", "telomere", "tRNA", "3'UTR", "5'UTR"]
            try:
                file_regions = handle_read.annotations['structured_comment']['Genome-Annotation-Data']['Features Annotated'].split("; ")
            except:
                #print("except region")
                file_regions = regions
            
            selected_regions = []
            intron_is_selected = False
            cds_is_selected = False
            intron_filename = ""
            cds_filename = ""
            filename = ""

            # we update the dataframe
            organism_df['features'][index_df] = [f for f in file_regions]
            # if 0 introns found, this will be deleted from the dataframe, so that we don't verify it again in
            # in the next executions
            if("CDS" in file_regions):
                organism_df['features'][index_df].append("intron")

            #print(region_choice)
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
                    log_signal.emit('Pass : {} does not contain selected option \'{}\''.format(NC, option), purple)
                    #print('Pass : {} does not contain selected option \'{}\''.format(NC, option))
            
            visited_regions = [False for i in selected_regions]

            count_complements = 0
            nb_introns = 0
            for f in features:
                if(self.check_stopping(mutex_stop)):
                    mutex_count.acquire()
                    count = count - 1
                    mutex_count.release()
                    return
                if f.type in selected_regions:
                    #log_signal.emit(NC + ": Parsing " + f.type)
                    ##print(NC+ ": Parsing " + f.type)
                    index = selected_regions.index(f.type)

                    if f.location:
                        if(self.manage_errors(f, len(handle_read.seq), log_signal)):
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
                                        os.remove(intron_filename)
                                    if cds_is_selected:
                                        os.remove(cds_filename)
                                else:
                                    os.remove(filename)
                            except:
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
                                            [intron_seq, nb] = self.join(header, handle_read.seq, f.location, log_signal, True)
                                            nb_introns += nb
                                        # gestion erreur
                                        except:
                                            continue

                                    if cds_is_selected:
                                        cds_seq = self.join(header, handle_read.seq, f.location, log_signal)[0]
                                        # gestion erreur
                                        if not cds_seq:
                                            continue
                                else:
                                    final_seq = self.join(header, handle_read.seq, f.location, log_signal)[0]
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
                                            [intron_seq, nb] = self.join(header, handle_read.seq, f.location, log_signal, True)
                                            nb_introns += nb
                                        # gestion erreur
                                            if(not intron_seq): continue
                                        except: continue
                                    if cds_is_selected:
                                        cds_seq = self.join(header, handle_read.seq, f.location, log_signal)[0]
                                        # gestion erreur
                                        if not cds_seq: continue
                                else:
                                    final_seq = self.join(header, handle_read.seq, f.location, log_signal)[0]
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
                            log_signal.emit('Pass ' + NC + ': noisy strand', purple)
                            #print('Error : noisy strand')
                            continue

                        else:
                            log_signal.emit('Pass ' + NC + ': cannot join complementory and normal strands', purple)
                            #print('Error : cannot join complementory and normal strands')

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
                if("intron" in organism_df['features'][index_df]):
                    organism_df['features'][index_df].remove("intron")
                try:
                    intron_file.close()
                    os.remove(intron_filename)
                except:
                    try:
                        os.remove(intron_filename)
                    except:
                        pass
                
            progress_signal.emit(0)
            log_signal.emit("Parsing of " + NC + " done successfully.", green)
            #print("number of introns found: {}".format(nb_introns))
        mutex_count.acquire()
        count = count - 1
        mutex_count.release()
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
                log_signal.emit("Pass: invalid join sequence order", purple)
                #print("Error : invalid join sequence order")
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

        if intron:
            header = header.replace("CDS", "intron")
        for (i, s) in enumerate(seq_join_intern):
            key_word = ' Intron ' if intron else ' Exon '
            final_seq += '\n' + header + key_word + str(i + 1) + '\n' + str(s)

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