from PyQt5 import QtCore
import time
import shutil
import pandas as pd
import os.path
from multiprocessing import Pool
import pickle
from ftp_downloader import *
from multiprocessing import Pool
import shutil
import os
import ftplib
import socket
from ftplib import error_temp
from time import sleep 
from parser_thread import bdd_path

save_pickle = False
DEBUG = False
VERBOSE = False
green = [0,255,0,255]
white = [255,255,255,255]
purple = [255,0,255,255] 
red = [255,0,0,255] 

class DownloaderThread(QtCore.QThread):
	
	progress_signal = QtCore.pyqtSignal(int)
	log_signal = QtCore.pyqtSignal(str, list)
	dataframe_result = QtCore.pyqtSignal(pd.core.frame.DataFrame)
	end_signal = QtCore.pyqtSignal(str)

	def __init__(self, parent=None, index = 0):
		super(DownloaderThread, self).__init__(parent)
		self.index = index
		self.parent = parent
		self.kingdoms_choice = ['Eukaryota','Archaea','Viruses','Bacteria']
		self.isRunning = False
		self.organism_df = 0


################################################################################
################################################################################


	def run(self):
		self.isRunning = True

		start_time = time.time()		
		
		self.remove_ill_terminated()
		self.load_tree()

		self.log_signal.emit("Téléchargement de l'arborescence en " + str(round(time.time() - start_time,3)) + " s.", green)
		self.end_signal.emit("Terminé avec succès.")

################################################################################
################################################################################

	def remove_ill_terminated(self):
		try:
			bdd = open(bdd_path,"r")
			lines = bdd.readlines()
			bdd.close()
			for l in lines:
				try:
					os.remove(l.strip('\n'))
				except: pass
			try: os.remove(bdd_path)
			except: pass
		except: pass

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

					organism = parsed_row[0].replace(' ','_').replace('/','_').replace('[','').replace(']','').replace(':','_').replace('\'','')
					group = parsed_row[2].replace(' ','_').replace('/','_')
					subgroup = parsed_row[3].replace(' ','_').replace('/','_')

					# organism_names is used only for comparing with the organisms names we find later in the IDS files
					# to get the right index in dataframe
					organism_names.append(parsed_row[0])

					organism_paths.append('../Results/' + kingdom +'/' + group +'/' + subgroup +'/' + organism +'/')

				except IndexError : pass

		# Parsing of the IDS
		ids_files = os.listdir('../GENOME_REPORTS/IDS/')
		if VERBOSE:
			print('overview done !')
		msg = "Overview terminée.\nRécupération des IDs à stocker dans le dataframe..."
		print(msg)
		self.log_signal.emit(msg,green) 

		#organism_names_ids = [] # we will store the organisms names here
		organism_paths_dataframe = []	# we will store the organisms paths here
		organism_NC_dataframe = []	# we will store the NC to parse here
		organism_names_dataframe = []
		i = 0
		found = 0

		# looping through the kingdoms ids files (viruses.ids, archaea.ids, etc ..)
		for ids in ids_files:
			i += 1

			msg = "Récupération des IDs des " + str(ids) + "..."
			self.log_signal.emit(msg, white)
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
					except:
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
					except ValueError:
						# We encouter this organism for the first time, we add it to the array with the NC
						#organism_names_ids.append(organism_names[index])
						try:
							organism_paths_dataframe.append(organism_paths[index])
							organism_NC_dataframe.append([parsed_row[1]])
							organism_names_dataframe.append(organism_names[index].replace(' ','_').replace('/','_').replace('[','').replace(']','').replace(':','_').replace('\'',''))
						except: pass
					

		# Store the organisms in pandas DataFrame
		organism_df = pd.DataFrame({
					"name":organism_names_dataframe,
					"path":organism_paths_dataframe,
					"NC":organism_NC_dataframe,
					"features":[ [] for i in organism_paths_dataframe]
					})

		# Create a pickle file to save the dataframe in local
		if not os.path.exists("../pickle"):
			os.makedirs("../pickle")
		with open("../pickle/organism_df", 'wb') as f:
			pickle.dump(organism_df, f)
	

	
################################################################################
################################################################################


	def load_df_from_pickle(self):

		msg = "Récupération des données depuis le dataframe..."
		print(msg)
		self.log_signal.emit(msg, white)

		try:
			with open("../pickle/organism_df", 'rb') as f:
				organism_df = pickle.load(f)
		except IOError:
			msg = "pickle file not accessible"
			print(msg)
			self.log_signal.emit(msg, red) 
			return

		msg = "Construction de l'arborescence..."
		print(msg)
		self.log_signal.emit(msg, white)
		for i in range(len(organism_df)):
			if not os.path.exists(organism_df["path"][i]):
				os.makedirs(organism_df["path"][i])

		msg = "Construction de l'arborescence terminée avec succès."
		print(msg)
		self.log_signal.emit(msg, green)

		self.dataframe_result.emit(organism_df)

		self.organism_df = organism_df



################################################################################
################################################################################


	def stop(self):

		self.isRunning = False
		print("Stopping thread...",self.index)
		msg = "Thread " + str(self.index) + " finished."
		self.log_signal.emit(msg, white)
		self.terminate()
		
	

########################################################################################################################
########################################################################################################################


	def load_tree(self):

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

		files = [("", "overview.txt"),("IDS","Eukaryota.ids"),("IDS", "Archaea.ids"), ("IDS", "Bacteria.ids"), ("IDS", "Viruses.ids")]
		files = [(directory,) + f for f in files]

		msg = "Téléchargement des fichiers..."
		print(msg)
		self.log_signal.emit(msg, white)
		with Pool(5) as p:
			p.map(download_ftp_file, files)


########################################################################################################################
########################################################################################################################

		