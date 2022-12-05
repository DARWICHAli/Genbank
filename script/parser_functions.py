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
from Bio.Align import MultipleSeqAlignment
from Bio.File import as_handle
from threading import Thread, Lock, get_ident


green = [0,255,0,255]
white = [255,255,255,255]
purple = [255,0,255,255]
red = [255,0,0,255]
yellow = [0, 255, 255]

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
    verbose = True
    # This is a generator function!
    def _SeqIO_to_alignment_iterator(cls, handle, format, mutex, seq_count=None):

        if format not in SeqIO._FormatToIterator:
            raise ValueError("Unknown format '%s'" % format)

        if seq_count:
            # Use the count to split the records into batches.
            seq_record_iterator = cls.parse(handle, format, mutex)
            if(seq_record_iterator is None):
                return None
            records = []
            for record in seq_record_iterator:
                if(cls.check_stopping()):
                    return None
                records.append(record)
                if len(records) == seq_count:
                    yield MultipleSeqAlignment(records)
                    records = []
            if records:
                raise ValueError("Check seq_count argument, not enough sequences?")
        else:
            # Must assume that there is a single alignment using all
            # the SeqRecord objects:
            records = list(cls.parse(handle, format, mutex))
            if records:
                yield MultipleSeqAlignment(records)
            else:
                return None

    def parse_2(cls, handle, format, mutex, seq_count=None):

        # Try and give helpful error messages:
        if not isinstance(format, str):
            raise TypeError("Need a string for the file format (lower case)")
        if not format:
            raise ValueError("Format required (lower case string)")
        if format != format.lower():
            raise ValueError("Format string '%s' should be lower case" % format)
        if seq_count is not None and not isinstance(seq_count, int):
            raise TypeError("Need integer for seq_count (sequences per alignment)")

        with as_handle(handle) as fp:
            # Map the file format to a sequence iterator:
            if format in _FormatToIterator:
                iterator_generator = _FormatToIterator[format]
                i = iterator_generator(fp, seq_count)

            elif format in SeqIO._FormatToIterator:
                # Exploit the existing SeqIO parser to the dirty work!
                i = cls._SeqIO_to_alignment_iterator(fp, format, mutex, seq_count=seq_count)
                if i is None:
                    return None
            else:
                raise ValueError("Unknown format '%s'" % format)

            yield from i

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
            v = cls.parse_2(handle, format, mutex)
            if v is None:
                return None
            for alignment in v:
                if cls.check_stopping(mutex):
                    return None
                for r in alignment:
                    return r
            return None
        raise ValueError("Unknown format '%s'" % format)


    def read(cls, handle, format, mutex, alphabet=None):
        iterator = cls.parse(handle, format, mutex, alphabet)
        if iterator is None:
            return None
        try:
            record = next(iterator)
            # print(f"## ANN: {record.annotations['organism']}")
            # print(f"## REC: {record.features}")
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

    def emit_log(self, log_signal, msg, color):
        if self.verbose:
            log_signal.emit(msg, color)


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

    def parse_NC(self, row_df, region_choice, log_signal, progress_signal, organism_df, mutex, mutex_fetch, mutex_count, mutex_stop):
        
        folder_color = white

        global count
        print(str(get_ident()))
        mutex_count.acquire()
        count = count + 1
        mutex_count.release()
        i = 0
        global stop
        index_df, organism_name, path, NC, file_features = row_df

        cp_region_choice = [i for i in region_choice]
        #for NC in NC_LIST:
        if(self.check_stopping(mutex_stop)):
            mutex_count.acquire()
            count = count - 1
            mutex_count.release()
            self.emit_log(log_signal, f'Pass {NC}: Mutex Stop', purple)
            return
        i+=1
        msg = "Parsing " + str(NC) + ' in: ' + str(path)
        log_signal.emit(msg, white)
        Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20)) + '@random.com'
        Entrez.api_key = 'd190df79af669bbea636278753c3236afa08'

        # Premier efetch pour récupérer juste la date.
        # le read prends prend plus de temps que la partie parsing, donc il faut vérifier la date avant
        # gestion timeout
        try:
            handle = Entrez.efetch(db="nucleotide", id=NC, rettype="gb", retmode="text", datetype='mdat')
        except:
            mutex_fetch.acquire()
            timer = 0.0
            while(timer != 5.0):
                if self.check_stopping(mutex_stop):
                    self.emit_log(log_signal, f'Pass {NC}: Mutex Stop', purple)
                    return
                timer += 0.1
            handle = Entrez.efetch(db="nucleotide", id=NC, rettype="gb", retmode="text", datetype='mdat')
            mutex_fetch.release()


        # On vérifie si le NC a été modidifé dans la base de données Genbank
        handle_date = handle.readline().strip().split(' ')
        handle_date = handle_date[len(handle_date)-1]
        # renvoie True si le NC a été modifié
        NC_modified_init = self.verify_modification_date(path, handle_date)
        NC_modified = NC_modified_init
        # nouvel Array qui va contenir les régions qu'on veut parser
        new_region_choice = []

        parse = False
        # Si le NC n'est pas modifié, on vérifie si l'une des options à parser n'existe pas dans les options déjà parsées
        # if(not NC_modified):
        #     for i in cp_region_choice:
        #         if(i not in file_features):
        #             parse = True
        #             break
        

        # A chaque tour de boucle, on vérifie si cette option a déja été parsée, donc que le fichier existe dans le dossier
        # Si le fichier n'existe pas ou que le NC a été modifié, on ajoute cette option aux options à parser
        # Sinon, on ne l'ajoute pas
        for option in cp_region_choice:
            # Si la colonne file_features a été mise à jour lors de l'execution précedente,
            # on vérifie que ce NC contient bien cette option, sinon ce n'est pas la peine de la chercher.
            if(len(file_features) and option not in file_features):
                continue
            # On vérifie si le NC a été modifie ou si la région a déjà été parsée
            filename = path + f"{option}_{organism_name}_{NC}.{i}.txt"
            NC_modified = NC_modified_init or not os.path.isfile(filename)
            if(NC_modified):
                new_region_choice.append(option)
            else:
                self.emit_log(log_signal, f"Skip: option {option} in {NC} is up to date.", purple)

        if(not len(new_region_choice)):
            self.emit_log(log_signal, f"Pass {NC}: All regions are up to date", green)
            ## EMIT SIGNAL 
            return

        # update region choice
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
                    mutex_count.acquire()
                    count = count - 1
                    mutex_count.release()
                    self.emit_log(log_signal, f'Pass {NC}: Mutex Stop', purple)
                    return
                timer += 0.1
            handle = Entrez.efetch(db="nucleotide", id=NC, rettype="gbwithparts", retmode="text", datetype='mdat')
            mutex_fetch.release()

        with warnings.catch_warnings(record=True) as w:
            if(not self.check_stopping(mutex_stop)):
                handle_read = self.read(handle, "gb", mutex_stop)
            else:
                handle_read = None
            #handle_read = SeqIO.read("test_input.gb", "gb")
            if handle_read is None:
                mutex_count.acquire()
                count = count - 1
                mutex_count.release()
                self.emit_log(log_signal, f'Pass {NC}: Mutex Stop', purple)
                return
            if w:
                for warning in w:
                    log_signal.emit(warning.message, purple)

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

        # we update the dataframe with the new extracted features for this NC
        organism_df['features'][index_df] = [f for f in file_regions]

        # if 0 introns found, this will be deleted from the dataframe, so that we don't verify it again in the next executions
        if("CDS" in file_regions):
            organism_df['features'][index_df].append("intron")

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
                self.emit_log(log_signal,f'Pass option {option} in {NC}: does not contain selected option', purple)

        if(not len(selected_regions)):
            self.emit_log(log_signal,f'Pass {NC}: does not contain any of the selected options', purple)
            return

        visited_regions = [False for i in selected_regions]

        count_complements = 0
        nb_introns = 0

        for f in features:
            if(self.check_stopping(mutex_stop)):
                mutex_count.acquire()
                count = count - 1
                mutex_count.release()
                self.emit_log(log_signal, 'Mutex Stop', purple)
                return
            if f.type in selected_regions:
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
                        self.emit_log(log_signal,'Pass ' + NC + ': noisy strand', purple)
                        continue

                    else:
                        self.emit_log(log_signal,'Pass ' + NC + ': cannot join complementory and normal strands', purple)

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

        
        self.emit_log(log_signal, f"Parsing of {NC} done successfully.", green)
        progress_signal.emit(0)

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
                self.emit_log(log_signal, "Pass: invalid join sequence order", purple)
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

        sep = ''
        if (len(seq_join_intern) > 1):
            final_seq += header + '\n' + "".join(seq_join_intern)
            sep = '\n'

        if intron:
            header = header.replace("CDS", "intron")
        for (i, s) in enumerate(seq_join_intern):
            key_word = ' Intron ' if intron else ' Exon '
            final_seq += sep + header + key_word + str(i + 1) + '\n' + str(s)

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
