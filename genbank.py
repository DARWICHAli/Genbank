#!/usr/bin/env python3
#from Bio import Entrez
import pandas as pd
import os
import shutil
#from Bio import SeqIO
import pickle
import random
import string
from ftplib import FTP
import time
import os.path
import datetime
import asyncio

from pip import main

save_pickle = False
DEBUG = False
VERBOSE = False

class Genbank:


	@classmethod
	def get_kingdom_choice(cls, mainwindow):
		selected_kingdoms = []
		if(mainwindow.checkBox_prokaryota.isChecked()):
			selected_kingdoms = selected_kingdoms + ["prokaryota"]
		if(mainwindow.checkBox_archaea.isChecked()):
			selected_kingdoms = selected_kingdoms + ["archaea"]
		if(mainwindow.checkBox_bacteria.isChecked()):
			selected_kingdoms = selected_kingdoms + ["bacteria"]
		if(mainwindow.checkBox_eukaryota.isChecked()):
			selected_kingdoms = selected_kingdoms + ["eukaryota"]
		if(mainwindow.inputKingdom.toPlainText() != ""):
			selected_kingdoms = selected_kingdoms + [mainwindow.inputKingdom.toPlainText()]
		return selected_kingdoms


################################################################################
################################################################################


	@classmethod
	def get_region_choice(cls, mainwindow):
		selected_regions = []
		if(mainwindow.checkBox_rrna.isChecked()):
			selected_regions = selected_regions + ["rrna"]
		if(mainwindow.checkBox_cds.isChecked()):
			selected_regions = selected_regions + ["cds"]
		if(mainwindow.checkBox_trna.isChecked()):
			selected_regions = selected_regions + ["trna"]
		if(mainwindow.checkBox_centromere.isChecked()):
			selected_regions = selected_regions + ["centromere"]
		if(mainwindow.checkBox_telomere.isChecked()):
			selected_regions = selected_regions + ["telomere"]
		if(mainwindow.checkBox_3utr.isChecked()):
			selected_regions = selected_regions + ["3utr"]
		if(mainwindow.checkBox_5utr.isChecked()):
			selected_regions = selected_regions + ["5utr"]
		if(mainwindow.checkBox_mobile_element.isChecked()):
			selected_regions = selected_regions + ["mobile element"]
		if(mainwindow.checkBox_mobile_ncrna.isChecked()):
			selected_regions = selected_regions + ["ncrna"]
		if(mainwindow.checkBox_mobile_intron.isChecked()):
			selected_regions = selected_regions + ["intron"]
		if(mainwindow.inputRegion.toPlainText() != ""):
			selected_regions = selected_regions + [mainwindow.inputRegion.toPlainText()]
		return selected_regions


################################################################################
################################################################################


	@classmethod
	def start(cls, mainwindow):
		print("start")
		cls.log(mainwindow, "start")
		selected_kingdoms = cls.get_kingdom_choice(mainwindow)
		selected_regions = cls.get_region_choice(mainwindow)
		print(selected_kingdoms)
		print(selected_regions)
		cls.log(mainwindow, "Selected Kingdoms:  \n" + str(selected_kingdoms))
		cls.log(mainwindow, "Selected Regions:  \n" + str(selected_regions))
		start_time = time.time()
		directory = os.getcwd()
		cls.reset_tree(directory, mainwindow)
		# cls.load_df_from_pickle(mainwindow)
		print("Arborescence DONE !")
		print("--- %s seconds ---" % (time.time() - start_time))
		cls.log(mainwindow, "Arborescence done in " + str((time.time() - start_time)) + " seconds ")

################################################################################
################################################################################


	@classmethod
	def pause(cls):
		print("pause")

################################################################################
################################################################################


	@classmethod
	def reset(cls):
		print("reset")


################################################################################
################################################################################


	@classmethod
	def reset_tree(cls, root, mainwindow):
		# Reset the tree stored locally and download new one from ftp server
		if os.path.exists('./GENOME_REPORTS/overview.txt') and os.path.isfile("./pickle/organism_df"):
			time = os.path.getmtime("./GENOME_REPORTS/overview.txt") #get the time of last modification
			pred_time = datetime.datetime.fromtimestamp(time)
			now = datetime.datetime.now()
			cls.log(mainwindow, "Time: " + str(now))

		# if pred_time.month == now.month and pred_time.day == now.day and pred_time.year == now.year : # No update in this case
				# root.destroy() # Exit
				# return

		try: shutil.rmtree('./GENOME_REPORTS') # remove all files/subdirectories
		except: pass
		os.mkdir("./GENOME_REPORTS")
		os.chdir('./GENOME_REPORTS')
		os.mkdir('IDS')

		with FTP('ftp.ncbi.nlm.nih.gov') as ftp:
			ftp.login()  # Connexion to the FTP server
			cls.log(mainwindow, "Connecting to FTP server")
			# Retrieve Genome overview and IDS
			cls.log(mainwindow, "Retrieving Genome overview and IDS")
			ftp.cwd('genomes/GENOME_REPORTS')
			ftp.retrbinary('RETR overview.txt', open("overview.txt",'wb').write)
			ftp.cwd('IDS')
			for filename in ftp.nlst():
				ftp.retrbinary('RETR '+ filename, open("IDS/" + filename, 'wb').write)
		if VERBOSE:
			print("Genome overview and IDS downloaded !")
		cls.log(mainwindow, "Genome overview and IDS downloaded !")


		# Parsing of "overview.txt"
		organism_names = []
		organism_paths = []

		os.chdir('../')
		with open('./GENOME_REPORTS/overview.txt') as f:
			first_row = True
			count_rows = 1
			for row in f:
				if DEBUG:
					print(count_rows, " / 59674")
				count_rows += 1
				if first_row:
					first_row=False
					continue
				parsed_row = row.split('\t')

				cls.log(mainwindow, "extracting tree components")
				cls.update_progress_bar(mainwindow, 5)
				# Extraction of the tree components
				try :
					organism = parsed_row[0].replace(' ','_').replace('/','_')
					kingdom = parsed_row[1].replace(' ','_').replace('/','_')
					group = parsed_row[2].replace(' ','_').replace('/','_')
					subgroup = parsed_row[3].replace(' ','_').replace('/','_')
					path = './Results/' + kingdom +'/' + group +'/' + subgroup +'/' + organism
					organism_names.append(parsed_row[0])
					organism_paths.append('./Results/' + kingdom +'/' + group +'/' + subgroup +'/')
				except IndexError : pass

		# Parsing of the IDS
		ids_files = os.listdir('./GENOME_REPORTS/IDS/')
		if VERBOSE:
			print('overview done !')
		cls.log(mainwindow, "Overview done !")
		cls.update_progress_bar(mainwindow, 15)
		organism_names_ids = []
		organism_paths_ids = []
		organism_NC_ids = []
		i = 0

		for ids in ids_files:
			i += 1
			#progressbar['value'] = 0
			#root.update_idletasks()
			#cls.log(mainwindow, "Downloading IDs ")
			with open('./GENOME_REPORTS/IDS/' + ids) as f:
				n_line = sum(1 for _ in f)
			with open('./GENOME_REPORTS/IDS/' + ids) as f:
				index = -1
				if VERBOSE:
					print("ids")
				for row in f:
					parsed_row = row.replace('\n', '').split('\t')

					if (parsed_row[1][0:2] != 'NC'): # We need only the NC
						continue
					
					try:
						index = organism_names.index(parsed_row[5]) # Retrieve Kingdom
					except ValueError:
						parsed_name = parsed_row[5].split(' ')[::-1]
						try_name = parsed_row[5]
						for word in parsed_name :
							try_name = try_name.replace(' '+word, '')
							try:
								index = organism_names.index(try_name)
								break
							except : pass

				#progressbar['value'] += 1/n_line*100
				#root.update_idletasks()
				#cls.update_progress_bar(mainwindow, 1/n_line*100)

				if index != -1:
					try:
						organism_NC_ids[organism_names_ids.index(organism_names[index])].append(parsed_row[1])
					except ValueError:
						organism_names_ids.append(organism_names[index])
						organism_paths_ids.append(organism_paths[index])
						organism_NC_ids.append([parsed_row[1]])
						name = organism_names[index].replace(" ", "_")
						name = name.replace("[", "_")
						name = name.replace("]", "_")
						name = name.replace(":", "_")
						path = organism_paths[index] + name + "/"
						if not os.path.exists(path):
							os.makedirs(path)

		# Store the organisms in pandas DataFrame
		organism_df = pd.DataFrame({
						"name":organism_names_ids,
						"path":organism_paths_ids,
						"NC":organism_NC_ids})

			# Create a pickle file to save the dataframe in local
		if not os.path.exists("./pickle"):
			os.makedirs("./pickle")
		with open("./pickle/organism_df", 'wb') as f:
			pickle.dump(organism_df, f)

			#root.destroy() # EXIT

################################################################################
################################################################################


	@classmethod
	def load_df_from_pickle(cls, mainwindow):
		# Loading pickle dataframe
		try:
			with open("./pickle/organism_df", 'rb') as f:
				organism_df = pickle.load(f)
		except IOError:
			print("Pickle file not accessible")
			directory = os.getcwd()
			cls.log(mainwindow, "Pickle File not accessible, resetting tree now.")
			return cls.reset_tree(directory, mainwindow)

		for i in range(len(organism_df)):
			name = organism_df["name"][i].replace(" ", "_")
			name = name.replace("[", "_")
			name = name.replace("]", "_")
			name = name.replace(":", "_")
			path = organism_df["path"][i] + name + "/"
			if not os.path.exists(path):
				os.makedirs(path)
		return organism_df

################################################################################
################################################################################

	# Logger
	@classmethod
	def log(cls, mainwindow, str):
		mainwindow.logOutput.insertPlainText(str+'\n')
		sb = mainwindow.logOutput.verticalScrollBar()
		sb.setValue(sb.maximum())
		
	@classmethod
	def update_progress_bar(cls, mainwindow, value):
		mainwindow.progressBar.setProperty("value", mainwindow.progressBar.value() + value)