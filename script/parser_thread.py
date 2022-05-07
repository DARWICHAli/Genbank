from PyQt5 import QtCore
import time
import shutil
import pandas as pd
import os.path
from ftp_downloader import *
from parser_class import ParserClass
import shutil
import os
import random
import string
import os
from Bio import SeqIO, Entrez
from Bio.SeqFeature import FeatureLocation
from datetime import datetime
import time
import warnings

save_pickle = False
DEBUG = False
VERBOSE = False


class ParserThread(QtCore.QThread):
	
	progress_signal = QtCore.pyqtSignal(int)
	log_signal = QtCore.pyqtSignal(str)
	dataframe_result = QtCore.pyqtSignal(pd.core.frame.DataFrame)
	end_signal = QtCore.pyqtSignal(str)
	

	def __init__(self, parent, index, organism_df, path_choice, regions_choice):
		super(ParserThread, self).__init__(parent)
		self.index = index
		self.parent = parent
		self.regions_choice = regions_choice
		self.path_choice =  path_choice
		self.organism_df = organism_df
		self.isRunning = False
		
		self.nb_NC = 0
		self.nb_parsed = 0
		self.current_file = 0
		self.current_path = ""

################################################################################
################################################################################


	def run(self):
		self.isRunning = True
		
		# resetting time
		start_time = time.time()

		print("Parsing Started...")
		parsing_choice = "->".join(self.path_choice.split('/')[2:])
		self.log_signal.emit("Parsing of " + parsing_choice + " started...")

		total = len([path for path in self.organism_df['path'] if path.startswith(self.path_choice)])
		print(type(total))
		self.nb_NC = total
		self.parent.mainwindow.progressBar.setFormat("0/"+str(total)+" NC")
		# About 157421 files to parse in total, we test with the first 10
		for (index, path, NC_LIST) in self.organism_df.itertuples():
			for NC in NC_LIST:
				if(path.startswith(self.path_choice)):
					self.current_path = path
					msg = "Parsing " + str(NC) + ' in: ' + str(path)
					self.log_signal.emit(msg)
					print(msg)
					self.parse_NC(NC, path, self.regions_choice, self.log_signal)
					self.progress_signal.emit(self.nb_NC)

		msg = "Parsing finished in: "
		self.log_signal.emit(msg)
		print(msg)
		self.log_signal.emit(str(time.time() - start_time))
		print(str(time.time() - start_time))
		self.end_signal.emit("----------- End -----------")


################################################################################
################################################################################


	def stop(self):
		# try:
		# 	self.current_file.close()
		# 	os.remove(self.current_path)
		# except:
		# 	try:
		# 		os.remove(self.current_path)
		# 	except:
		# 		pass

		self.isRunning = False
		print("Stopping thread...",self.index)
		msg = "Stopping thread..." + str(self.index)
		self.log_signal.emit(msg)
		self.terminate()


########################################################################################################################
##############################################  Parsing Functions ######################################################


	def manage_errors(self, f, len_seq, signal):
		if f.location.start < 0:
			signal.emit("Pass: sequence must not start with 0 or less. {}".format(f.location))
			print("Pass: sequence must not start with 0 or less. {}".format(f.location))
			return True
		elif f.location.start > len_seq or f.location.end > len_seq:
			signal.emit("Pass: sequence borders must be less or equal than total sequence size. {}".format(f.location))
			print("Pass: sequence borders must be less or equal than total sequence size. {}".format(f.location))
			return True
		#elif str(f.location.start)[0] == '<' or str(f.location.end)[0] == '>'\
		#        or str(f.location.start)[0] == '>' or str(f.location.end)[0] == '<':
		try:
			int(str(f.location.start))
			int(str(f.location.end))
		except:
			signal.emit("Pass: only numerical arguments are accepted for sequence limits. {}".format(f.location))
			print("Pass: only numerical arguments are accepted for sequence limits. {}".format(f.location))
			return True
		return False

	def parse_NC(self, NC, path, region_choice, signal):
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
					
		handle_date = handle_read.annotations['date']
		NC_modified_init = self.verify_modification_date(path, handle_date)
		NC_modified = NC_modified_init

		organism = handle_read.annotations['organism']
		features = handle_read.features

		# delete 'intron' from regions beceause it's not actually a region but just a parsing option
		regions = ["CDS", "centromere", "mobile_element", "ncRNA", "rRNA", "telomere", "tRNA", "3'UTR", "5'UTR"]
		try:
			file_regions = handle_read.annotations['structured_comment']['Genome-Annotation-Data']['Features Annotated'].split("; ")
		except:
			print("except region")
			file_regions = regions
		
		selected_regions = []
		intron_is_selected = False
		cds_is_selected = False
		intron_filename = ""
		cds_filename = ""
		filename = ""
		
		
		for option in region_choice:
			if option in file_regions:
				if(option == 'CDS'):
					cds_is_selected = True
				selected_regions.append(option)
			elif option == "intron":
				if "CDS" not in region_choice:
					selected_regions.append('CDS')
				intron_is_selected = True
			else:
				signal.emit('Pass : {} does not contain selected option \'{}\''.format(NC, option))
				print('Pass : {} does not contain selected option \'{}\''.format(NC, option))
		
		for option in file_regions:
			if option not in selected_regions:
				file_regions.remove(option)
		
		visited_regions = [False for i in selected_regions]
		createdNow = [False for i in selected_regions]

		count_complements = 0
		nb_introns = 0
		for f in features:
			if f.type in file_regions:
				index = selected_regions.index(f.type)
				NC_modified = NC_modified_init
				if f.location:
					if(self.manage_errors(f, len(handle_read.seq), signal)):
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
							if(not os.path.isfile(intron_filename)):
								createdNow[index] = True
							NC_modified = NC_modified or not os.path.isfile(intron_filename) or createdNow[index]
						if cds_is_selected:
							cds_filename = path + "/CDS_{}_{}.txt".format(ogranism_str, handle_read.id)
							if(not os.path.isfile(cds_filename)):
								createdNow[index] = True
							NC_modified = NC_modified or not os.path.isfile(cds_filename) or createdNow[index]
					else:
						filename = path + "/{}_{}_{}.txt".format(f.type, ogranism_str, handle_read.id)
						if(not os.path.isfile(filename)):
							createdNow[index] = True
						NC_modified = NC_modified or not os.path.isfile(filename) or createdNow[index]
					
					if(not NC_modified):
						self.log_signal.emit("Pass: " + NC + "_" + f.type + " is up to date.")
						file_regions.remove(f.type)						
						continue

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
										[intron_seq, nb] = self.join(header, handle_read.seq, f.location, signal, True)
										nb_introns += nb
									# gestion erreur
									except:
										continue

								if cds_is_selected:
									cds_seq = self.join(header, handle_read.seq, f.location, signal)[0]
									# gestion erreur
									if not cds_seq:
										continue
							else:
								final_seq = self.join(header, handle_read.seq, f.location, signal)[0]
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
										[intron_seq, nb] = self.join(header, handle_read.seq, f.location, signal, True)
										nb_introns += nb
									# gestion erreur
										if(not intron_seq): continue
									except: continue
								if cds_is_selected:
									cds_seq = self.join(header, handle_read.seq, f.location, signal)[0]
									# gestion erreur
									if not cds_seq: continue
							else:
								final_seq = self.join(header, handle_read.seq, f.location, signal)[0]
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
			try:
				intron_file.close()
				os.remove(intron_file.name)
			except:
				try:
					os.remove(intron_file.name)
				except:
					pass

		#print("number of introns found: {}".format(nb_introns))
		return True




	def join(self, header, sequence, location, signal, intron=False):
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



		