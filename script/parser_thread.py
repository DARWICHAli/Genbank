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
green = [0,255,0,255]
white = [255,255,255,255]
purple = [255,0,255,255]
red = [255,0,0,255]

from threading import Thread, Lock


class ParserThread(QtCore.QThread):

	progress_signal = QtCore.pyqtSignal(int)
	log_signal = QtCore.pyqtSignal(str, list)
	dataframe_result = QtCore.pyqtSignal(pd.core.frame.DataFrame)
	end_signal = QtCore.pyqtSignal(str)
	quick_info_signal = QtCore.pyqtSignal(str, list)
	fin_parsing_signal = QtCore.pyqtSignal()

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
		self.parent.terminate_parsing_signal.connect(self.terminate)
		self.end_parsing = False

################################################################################
################################################################################


	def run(self):
		self.isRunning = True

		# resetting time
		start_time = time.time()

		parsing_choice = " >> ".join(self.path_choice.split('/')[2:])
		self.log_signal.emit("PARSING des " + parsing_choice + " ...", [255,255,255,255])

		self.pool=Pool()

		args =  self.organism_df.loc[self.organism_df['path'].str.startswith(self.path_choice + '/')]
		self.nb_NC = sum( len(i) for i in args['NC'])
		self.progress_signal.emit(self.nb_NC)

		self.pool.map(partial(self.parser.parse_NC, region_choice = self.regions_choice, log_signal = self.log_signal, progress_signal = self.progress_signal, organism_df = self.organism_df, mutex = self.mutex, mutex_fetch = self.mutex_fetch, mutex_count = self.mutex_count, mutex_stop = self.mutex_stop), args.itertuples())


		# new dataframe with file features added
		with open("../pickle/organism_df", 'wb') as f:
			pickle.dump(self.organism_df, f)

		try: os.remove(bdd_path)
		except: pass

		msg = "Parsing terminée en: "
		self.log_signal.emit(msg, green)
		self.log_signal.emit(str(round(time.time() - start_time,3)), green)

		if VERBOSE:
			print(str(time.time() - start_time))

		self.end_signal.emit("----------- FIN -----------")



################################################################################
################################################################################


	def stop(self):
		self.isRunning = False
		self.pool.close()
		self.parser.set_stop(self.mutex_count)
		print("Veuillez attendre la fin des threads.")
		while(self.parser.get_count(self.mutex_count)):
			print("Il reste {} thread".format(self.parser.get_count(self.mutex_count)))
			self.quick_info_signal.emit("Veuillez attendre la fin des threads.\nIl reste {} thread".format(self.parser.get_count(self.mutex_count)), purple)
			time.sleep(2)
		self.quick_info_signal.emit("Threads Terminées.", green)
		print("Stopping thread...",self.index)
		msg = "Thread " + str(self.index) + " terminée"
		self.log_signal.emit(msg, white)
		self.fin_parsing_signal.emit()
		self.parser.reinit()
		self.terminate()


########################################################################################################################
##############################################  Parsing Functions ######################################################
