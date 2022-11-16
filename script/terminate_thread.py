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


class TerminateThread(QtCore.QThread):
	
	def __init__(self, parent, thread):
		super(TerminateThread, self).__init__(parent)
		self.parent = parent
		self.thread = thread
		self.stop_thread = False

################################################################################
################################################################################


	def run(self):
		while(True):
			if(self.stop_thread):
				self.thread.stop()
				break


################################################################################
################################################################################
		

	def stop(self):
		self.stop_thread = True

########################################################################################################################
##############################################  Parsing Functions ######################################################

		