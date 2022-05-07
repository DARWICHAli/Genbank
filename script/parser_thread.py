from signal import signal
from PyQt5 import QtCore
import time
import pickle
import pandas as pd
from ftp_downloader import *
import time
from parser_functions import ParserFunctions, bdd_path
import threading
from functools import partial
from itertools import repeat
from multiprocessing.pool import ThreadPool as Pool

save_pickle = False
DEBUG = False
VERBOSE = False

from threading import Thread, Lock


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
		self.mutex = Lock()
		self.mutex_fetch = Lock()
		self.mutex_count = Lock()
		self.mutex_stop = Lock()
		self.nb_NC = 0
		self.nb_parsed = 0
		self.current_file = 0
		self.current_path = ""
		self.parser = ParserFunctions()

################################################################################
################################################################################


	def run(self):
		self.isRunning = True
		
		# resetting time
		start_time = time.time()

		print("Parsing Started...")
		parsing_choice = " >> ".join(self.path_choice.split('/')[2:])
		self.log_signal.emit("Parsing of " + parsing_choice + " started...")

		self.pool=Pool()

		args =  self.organism_df.loc[self.organism_df['path'].str.startswith(self.path_choice + '/')]
		self.nb_NC = len(args['NC'])
		self.progress_signal.emit(self.nb_NC)
		
		self.pool.map(partial(self.parser.parse_NC, region_choice = self.regions_choice, log_signal = self.log_signal, progress_signal = self.progress_signal, organism_df = self.organism_df, mutex = self.mutex, mutex_fetch = self.mutex_fetch, mutex_count = self.mutex_count, mutex_stop = self.mutex_stop), args.itertuples())

		# About 157421 files to parse in total, we test with the first 10
		# for (index_df, organism, path, NC_LIST, file_features) in args.itertuples():

		# 	parse_NC(index_df, organism, path, NC_LIST, file_features, self.regions_choice, self.log_signal, self.organism_df)
		# 	self.progress_signal.emit(self.nb_NC)

		# new dataframe with file features added
		with open("../pickle/organism_df", 'wb') as f:
			pickle.dump(self.organism_df, f)

		try: os.remove(bdd_path)
		except: pass

		msg = "Parsing finished in: "
		self.log_signal.emit(msg)
		self.log_signal.emit(str(time.time() - start_time))
		print(str(time.time() - start_time))
		self.end_signal.emit("----------- End -----------")


################################################################################
################################################################################


	def stop(self):
		self.isRunning = False
		self.parser.set_stop(self.mutex_count)
		while(self.parser.get_count(self.mutex_count)):
			print("Il reste {} thread".format(self.parser.get_count(self.mutex_count)))
			time.sleep(2)

		# while(len(self.pool)):
		# 	time.sleep(0.5)
		print("Stopping thread...",self.index)
		msg = "Stopping thread..." + str(self.index)
		self.log_signal.emit(msg)
		self.parser.reinit()
		self.terminate()


########################################################################################################################
##############################################  Parsing Functions ######################################################

		