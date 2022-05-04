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
from parser_class import ParserClass
import asyncio

save_pickle = False
DEBUG = False
VERBOSE = False


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
		self.regions_choice = []
		self.kingdoms_choice = []
		parent.region_signal.connect(self.get_region_choice)
		parent.kingdom_signal.connect(self.get_kingdom_choice)
		self.isRunning = False
		self.organism_df = 0
		self.nb_NC = 0

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

		self.nb_NC = sum( len(n) for n in self.organism_df["NC"] )
		print(self.nb_NC)
		# About 157421 files to parse in total, we test with the first 10
		self.nb_NC=100 
		for (index, path, NC_LIST) in self.organism_df.itertuples():
			# We only parse kingdoms selected by user
			kingdom = path.split('/')[2]
			if(kingdom not in self.kingdoms_choice): continue

			for NC in NC_LIST:	
				if(nb_parsed==50): break
				msg = "Parsing " + str(NC) + '...\n In: ' + str(path)
				self.any_signal.emit(msg)
				print(msg)
				if(ParserClass.parse_NC(NC, path, self.regions_choice) == False):
					msg = "Erreur Parsing " + str(NC) + ". Fichier supprimé."
				else:
					msg = "Parsing de " + str(NC) + " réussis."
				self.any_signal.emit(msg)
				print(msg)
				nb_parsed+=1
				self.progress_signal.emit(100./self.nb_NC)

		msg = "Parsing finished in: "
		self.any_signal.emit(msg)
		print(msg)
		self.time_signal.emit(time.time() - start_time)
		print(str(time.time() - start_time))
		self.end_signal.emit("End")

################################################################################
################################################################################

	def get_region_choice(self, regions_choice):
		index = self.sender().index
		if(index == 0):
			self.regions_choice = regions_choice


	def get_kingdom_choice(self, kingdoms_choice):
		self.get_kingdom_choice = []
		index = self.sender().index
		if(index == 0):
			self.kingdoms_choice = kingdoms_choice

################################################################################
################################################################################

	def download_files(self):
		# Parsing of "overview.txt"
		organism_names = []
		organism_paths = []

		# Retrieving the information from overview.txt to build the tree
		with open('../GENOME_REPORTS/overview.txt') as f:
			# we drop the header
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
					kingdom = parsed_row[1].replace(' ','_').replace('/','_')
					# we only want the kingdoms that were selected by the user
					if(kingdom not in self.kingdoms_choice): 
						continue
					organism = parsed_row[0].replace(' ','_').replace('/','_')
					organism = organism.replace('[','').replace(']','').replace(':','_').replace('\'','')
					group = parsed_row[2].replace(' ','_').replace('/','_')
					subgroup = parsed_row[3].replace(' ','_').replace('/','_')

					# organism_names is used only for comparing with the organisms names we find later in the IDS files
					# to get the right index in dataframe
					organism_names.append(parsed_row[0])

					organism_paths.append('../Results/' + kingdom +'/' + group +'/' + subgroup +'/' + organism)

				except IndexError : pass

		# Parsing of the IDS
		ids_files = os.listdir('../GENOME_REPORTS/IDS/')
		if VERBOSE:
			print('overview done !')
		msg = "Overview Done. Parsing IDS to store in pickle dataframe..."
		print(msg)
		self.any_signal.emit(msg) 

		#organism_names_ids = [] # we will store the organisms names here
		organism_paths_dataframe = []	# we will store the organisms paths here
		organism_NC_dataframe = []	# we will store the NC to parse here
		i = 0
		found = 0
		not_found = 0
		# looping through the kingdoms ids files (viruses.ids, archaea.ids, etc ..)
		for ids in ids_files:
			i += 1

			msg = "Parsing " + str(ids) + "..."
			self.any_signal.emit(msg)
			with open('../GENOME_REPORTS/IDS/' + ids) as f:
				#n_line = sum(1 for _ in f)
				if VERBOSE:
					print("ids")
				for row in f:
					parsed_row = row.replace('\n', '').split('\t')

					if (parsed_row[1][0:2] != 'NC'): # We need only the NC
						continue
					try:
						# we need the index in the previous array that we built from overview
						# so that we can put the NC and path in the same index in our dataframe
						index = organism_names.index(parsed_row[5])
						found += 1
					except ValueError:
						# Organism does not exist in overview, we move to next NC
						found = False
						#continue
						parsed_name = parsed_row[5].split(' ')[::-1]
						try_name = parsed_row[5]
						for word in parsed_name :
							try_name = try_name.replace(' '+word, '')
							try:
								index = organism_names.index(try_name)
								found = True
								break
							except: 
								found = False
						if not found: continue

					try:
						# checking if we have already added this organism to the array
						# If so, we simply add another NC
						organism_NC_dataframe[organism_paths_dataframe.index(organism_paths[index])].append(parsed_row[1])
						self.nb_NC += 1
					except ValueError:
						# We encouter this organism for the first time, we add it to the array with the NC
						#organism_names_ids.append(organism_names[index])
						try:
							organism_paths_dataframe.append(organism_paths[index])
							organism_NC_dataframe.append([parsed_row[1]])
							self.nb_NC += 1
						except: pass
					
		print("found: " + str(found))
		print("not found " + str(not_found))
		# Store the organisms in pandas DataFrame
		organism_df = pd.DataFrame({
					#"name":organism_names_ids,
					"path":organism_paths_dataframe,
					"NC":organism_NC_dataframe})

		# Create a pickle file to save the dataframe in local
		if not os.path.exists("../pickle"):
			os.makedirs("../pickle")
		with open("../pickle/organism_df", 'wb') as f:
			pickle.dump(organism_df, f)
	

	
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

		msg = "Creating hierarchy..."
		print(msg)
		self.any_signal.emit(msg)
		for i in range(len(organism_df)):
			# name = organism_df["name"][i].replace(" ", "_")
			# name = name.replace("[", "_")
			# name = name.replace("]", "_")
			# name = name.replace(":", "_")
			# path = organism_df["path"][i] + name + "/"
			if not os.path.exists(organism_df["path"][i]):
				os.makedirs(organism_df["path"][i])

		msg = "Finished loading hierarchy"
		print(msg)
		self.any_signal.emit(msg)

		self.dataframe_result.emit(organism_df)

		self.organism_df = organism_df



################################################################################
################################################################################


	def stop(self):
		self.parent.region_signal.disconnect(self.get_region_choice)
		self.parent.kingdom_signal.disconnect()
		self.isRunning = False
		print("Stopping thread...",self.index)
		msg = "Stopping thread..." + str(self.index)
		self.any_signal.emit(msg)
		self.terminate()
		
	

########################################################################################################################
########################################################################################################################


	def load_tree(self):

		# TO DO: Add a to download array containing only the kingdoms to update instead of updating everything again
		for f in self.kingdoms_choice:
			if not os.path.isfile("../GENOME_REPORTS/IDS/" + f + '.ids'):
				self.download_ftp()
				self.download_files()
				self.load_df_from_pickle()
				return

		if os.path.isdir("../pickle") and os.path.isfile("../pickle/organism_df"):
			# Initialization
			ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
			ftp.login()
			ftp.cwd('genomes/GENOME_REPORTS/IDS')

			last_ftp_change = 0
			last_local_change = 1e50

			# ftp files last modification timestamp
			for f in self.kingdoms_choice:
				remote_datetime = ftp.voidcmd("MDTM " + f +'.ids')[4:].strip()
				remote_timestamp = time.mktime(time.strptime(remote_datetime, '%Y%m%d%H%M%S'))
				if int(remote_timestamp) > int(last_ftp_change):
					last_ftp_change = remote_timestamp

			# local files timestamp
			last_local_change = os.path.getmtime("../pickle/organism_df")

			# Download the newest files and create the tree
			if int(last_ftp_change) > int(last_local_change):
				self.download_ftp()
				self.download_files()
			
		else:
			self.download_ftp()
			self.download_files()
		
		self.load_df_from_pickle()
		#print("loaded from ftp (file or directory doesn't exist)")



########################################################################################################################
########################################################################################################################


	def download_ftp(self):

		try: shutil.rmtree('../GENOME_REPORTS') # remove all files/subdirectories
		except: 
			print("exception in download_ftp(): rmtree  ../GENOME_REPORTS failed.")

		try:
			os.mkdir("../GENOME_REPORTS")
		except: pass
		
		os.chdir('../GENOME_REPORTS')

		directory = os.getcwd()

		# Pourquoi ce test ? IDS n'a pas été supprimé en supprimant GENOME_REPORTS ?
		if os.path.exists(os.path.join(directory, "IDS")):
			shutil.remove(os.path.join(directory, "IDS"))
		os.mkdir(os.path.join(directory, "IDS"))

		files = [("", "overview.txt")]

		for k in self.kingdoms_choice:
			files.append(("IDS", k + '.ids'))

		# Adding destination
		files = [(directory,) + f for f in files]

		msg = "Downloading files..."
		print(msg)
		self.any_signal.emit(msg)
		print(files)
		with Pool(5) as p:
			p.map(download_ftp_file, files)
