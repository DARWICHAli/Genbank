from PyQt5 import QtCore
import datetime
import time
import os
import pandas as pd
import os.path
import pickle
import shutil
from ftplib import FTP
from ftp_downloader import *
from parser import Parser
import asyncio

save_pickle = False
DEBUG = False
VERBOSE = False
kingdoms_choice = []
regions_choice = []

class ThreadClass(QtCore.QThread):
	
	progress_signal = QtCore.pyqtSignal(float)
	any_signal = QtCore.pyqtSignal(str)
	dataframe_result = QtCore.pyqtSignal(pd.core.frame.DataFrame)
	time_signal = QtCore.pyqtSignal(float)
	end_signal = QtCore.pyqtSignal(str)
	

	def __init__(self, parent=None, index = 0):
		super(ThreadClass, self).__init__(parent)
		self.index = index
		self.parent = parent
		self.region_choice = []
		parent.region_signal.connect(self.get_region_choice)
		self.organism_df = 0
		self.isRunning = False

################################################################################
################################################################################


	def run(self):
		self.isRunning = True
		nb_parsed=0

		start_time = time.time()

		self.load_tree()
		
		msg = "Download overview and IDS time : "
		self.any_signal.emit(msg)
		self.time_signal.emit(time.time() - start_time)

		start_time = time.time()

		msg = "Started Parsing files..."
		print(msg)
		self.any_signal.emit(msg)

		# About 157421 files to parse in total, we test with the first 10
		for (index, names, path, NC_LIST) in self.organism_df.itertuples():
			for NC in NC_LIST:
				
				if(nb_parsed==10): break

				msg = "Parsing " + str(NC) + '...\n In: ' + str(path)
				self.any_signal.emit(msg)
				print(msg)
				if(Parser.parse_NC(NC, path, self.region_choice) == False):
					msg = "Erreur Parsing " + str(NC) + ". Fichier supprimé."
				else:
					msg = "Parsing de " + str(NC) + " réussis."
				self.any_signal.emit(msg)
				print(msg)
				nb_parsed+=1

		msg = "Parsing finished in: "
		self.any_signal.emit(msg)
		print(msg)
		self.time_signal.emit(time.time() - start_time)
		print(str(time.time() - start_time))
		self.end_signal.emit("End")

################################################################################
################################################################################

	def get_region_choice(self, region_choice):
		self.region_choice = []
		index = self.sender().index
		if(index == 0):
			self.region_choice = region_choice
			print(self.region_choice)

	
################################################################################
################################################################################

	def download_files(self):
		# Parsing of "overview.txt"
		organism_names = []
		organism_paths = []

		#os.chdir('../script')
		with open('../GENOME_REPORTS/overview.txt') as f:
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

				# Extraction of the tree components
				try :
					organism = parsed_row[0].replace(' ','_').replace('/','_')
					kingdom = parsed_row[1].replace(' ','_').replace('/','_')
					group = parsed_row[2].replace(' ','_').replace('/','_')
					subgroup = parsed_row[3].replace(' ','_').replace('/','_')
					path = '../Results/' + kingdom +'/' + group +'/' + subgroup +'/' + organism + '/' 
					organism_names.append(parsed_row[0])
					organism_paths.append(path)
				except IndexError : pass

		# Parsing of the IDS
		ids_files = os.listdir('../GENOME_REPORTS/IDS/')
		if VERBOSE:
			print('overview done !')
		msg = "Overview Done. Retrieving IDS..."
		print(msg)
		self.any_signal.emit(msg) 
		organism_names_ids = []
		organism_paths_ids = []
		organism_NC_ids = []
		i = 0
		for ids in ids_files:
			i += 1

			msg = "Retrieving " + str(ids) + "..."
			self.progress_signal.emit(1)
			self.any_signal.emit(msg)
			with open('../GENOME_REPORTS/IDS/' + ids) as f:
				n_line = sum(1 for _ in f)
			with open('../GENOME_REPORTS/IDS/' + ids) as f:
				if VERBOSE:
					print("ids")
				for row in f:
					parsed_row = row.replace('\n', '').split('\t')

					if (parsed_row[1][0:2] != 'NC'): # We need only the NC
						continue
					try:
						index = organism_names.index(parsed_row[5])
					except ValueError:
						parsed_name = parsed_row[5].split(' ')[::-1]
						try_name = parsed_row[5]
						for word in parsed_name :
							try_name = try_name.replace(' '+word, '')
							try:
								index = organism_names.index(try_name)
								break
							except : pass

					try:
						organism_NC_ids[organism_names_ids.index(organism_names[index])].append(parsed_row[1])
					except ValueError:
						organism_names_ids.append(organism_names[index])
						organism_paths_ids.append(organism_paths[index])
						organism_NC_ids.append([parsed_row[1]])
						# name = organism_names[index].replace(" ", "_")
						# name = name.replace("[", "_")
						# name = name.replace("]", "_")
						# name = name.replace(":", "_")
						path = organism_paths[index] 
						
						if not os.path.exists(path):
							os.makedirs(path)

		# Store the organisms in pandas DataFrame
		organism_df = pd.DataFrame({
					"name":organism_names_ids,
					"path":organism_paths_ids,
					"NC":organism_NC_ids})

		# Create a pickle file to save the dataframe in local
		if not os.path.exists("../pickle"):
			os.makedirs("../pickle")
		with open("../pickle/organism_df", 'wb') as f:
			pickle.dump(organism_df, f)
		
		msg = "Loading data from pickle..."
		print(msg)
		self.any_signal.emit(msg) 

	
################################################################################
################################################################################


	def load_df_from_pickle(self):

		msg = "Loading dataframe from pickle..."
		print(msg)
		self.any_signal.emit(msg)

		try:
			with open("../pickle/organism_df", 'rb') as f:
				organism_df = pickle.load(f)
		except IOError:
			msg = "pickle file not accessible"
			print(msg)
			self.any_signal.emit(msg) 
			return

		for i in range(len(organism_df)):
			name = organism_df["name"][i].replace(" ", "_")
			name = name.replace("[", "_")
			name = name.replace("]", "_")
			name = name.replace(":", "_")
			path = organism_df["path"][i] + "/"
			if not os.path.exists(path):
				os.makedirs(path)

		msg = "Finished loading hierarchy"
		print(msg)
		self.any_signal.emit(msg)

		self.dataframe_result.emit(organism_df)

		return organism_df



################################################################################
################################################################################


	def stop(self):
		self.parent.region_signal.disconnect(self.get_region_choice)
		self.isRunning = False
		print("Stopping thread...",self.index)
		msg = "Stopping thread..." + str(self.index)
		self.any_signal.emit(msg)
		self.terminate()
		
	

########################################################################################################################
########################################################################################################################


	def load_tree(self):
		files = ["Bacteria.ids", "Eukaryota.ids", "Archaea.ids", "Viruses.ids"]

		if os.path.isdir("../pickle") and os.path.isfile("../pickle/organism_df"):
			# Initialization
			ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
			ftp.login()
			ftp.cwd('genomes/GENOME_REPORTS/IDS')

			last_ftp_change = 0
			last_local_change = 1e50

			# ftp files last modification timestamp
			for f in files:
				remote_datetime = ftp.voidcmd("MDTM " + f)[4:].strip()
				remote_timestamp = time.mktime(time.strptime(remote_datetime, '%Y%m%d%H%M%S'))
				if int(remote_timestamp) > int(last_ftp_change):
					last_ftp_change = remote_timestamp
			#print(last_ftp_change)

			# local files timestamp
			last_local_change = os.path.getmtime("../pickle/organism_df")
			#print(last_local_change)

			# Download the newest files and create the tree
			if int(last_ftp_change) > int(last_local_change):
				self.download_ftp()
				self.download_files()
				#print("loaded from ftp")
			else:
				self.organism_df = self.load_df_from_pickle()
				#print("loaded from pickle")
		else:
			self.download_ftp()
			self.download_files()
			#print("loaded from ftp (file or directory doesn't exist)")



########################################################################################################################
########################################################################################################################


	def download_ftp(self):

		try: shutil.rmtree('../GENOME_REPORTS') # remove all files/subdirectories
		except: pass
		os.mkdir("../GENOME_REPORTS")
		os.chdir('../GENOME_REPORTS')

		directory = os.getcwd()

		if os.path.exists(os.path.join(directory, "IDS")):
			shutil.remove(os.path.join(directory, "IDS"))
		os.mkdir(os.path.join(directory, "IDS"))

		files = [("", "overview.txt"), ("IDS", "Bacteria.ids"),
				("IDS", "Eukaryota.ids"), ("IDS", "Archaea.ids"),
				("IDS", "Viruses.ids")]

		# Adding destination
		files = [(directory,) + f for f in files]

		msg = "Downloading files..."
		print(msg)
		self.any_signal.emit(msg)

		with Pool(5) as p:
			p.map(download_ftp_file, files)
