from signal import signal
from PyQt5 import QtCore
import time
import pickle
import pandas as pd
from ftp_downloader import *
import time
from parser_functions import ParserFunctions, bdd_path
from functools import partial
from multiprocessing.pool import Pool
from multiprocessing import Queue, Process, Event
import psutil

save_pickle = False
DEBUG = False
VERBOSE = False
green = [0,255,0,255]
white = [255,255,255,255]
purple = [255,0,255,255]
red = [255,0,0,255]

from multiprocessing import Lock,active_children
from functools import partial

def init_pool_processes(m, m_f, rc, e, q):
	global mutex
	global mutex_fetch
	global region_choice
	global queue
	global stop_event
	mutex = m
	mutex_fetch = m_f
	region_choice = rc
	stop_event = e
	queue = q

def my_parser(args):
    parser = ParserFunctions(region_choice, mutex, mutex_fetch, queue, stop_event)
    parser.parse_NC(args)

def handle_parser(mutex, mutex_fetch, regions_choice, stop_event, queue, args):
	pool = Pool(initializer=init_pool_processes, initargs=(mutex, mutex_fetch, regions_choice, stop_event, queue))
	pool.map(partial(my_parser), args)
	pool.close()
	pool.join()
	queue.put({"type":0})


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
		self.nb_NC = 0
		self.nb_parsed = 0
		self.current_file = 0
		self.current_path = ""
		self.parent.terminate_parsing_signal.connect(self.terminate)
		self.end_parsing = False

################################################################################
################################################################################


	def run(self):
		self.isRunning = True

		# resetting time
		start_time = time.time()
		parsing_choice = " >> ".join(self.path_choice.split('/')[2:])
		self.log_signal.emit(f"PARSING des {parsing_choice} ...", [255,255,255,255])
		queue = Queue()
		self.stop_event = Event()
		rows =  self.organism_df.loc[self.organism_df['path'].str.startswith(self.path_choice + '/')]
		args = [tuple(r) for r in rows.itertuples()] #if "NC_060948" in r[3]]
		self.nb_NC = len(args)
		self.progress_signal.emit(self.nb_NC)
		self.parser_handler = Process(target=handle_parser, args=(self.mutex, self.mutex_fetch, self.regions_choice, self.stop_event, queue, args))
		self.parser_handler.start()
		while(True):
			data = queue.get(block=True)
			if data['type'] == 1:
				self.log_signal.emit(data['msg'], data['color'])
			elif data['type'] == 2:
				self.progress_signal.emit(0)
			elif data['type'] == 3:
				if("intron" in self.organism_df[data['col']][data['row']]):
					#self.organism_df[data['col']][data['row']].remove("intron")
					self.organism_df[data['col']][data['row']] = list(filter(('intron').__ne__, self.organism_df[data['col']][data['row']]))
			elif data['type'] == 4:
				self.organism_df[data['col']][data['row']].append("intron")
			elif data['type'] == 5:
				self.organism_df[data['col']][data['row']] = data['value']
			elif data['type'] == 6:
				self.organism_df[data['col']][data['row']] = list(filter((data['value']).__ne__, self.organism_df[data['col']][data['row']]))
			else:
				break
		# new dataframe with file features added
		with open("../pickle/organism_df", 'wb') as f:
			pickle.dump(self.organism_df, f)

		# Deleting Update Path ... à voir
		try: os.remove(bdd_path)
		except: pass

		msg = "Parsing terminée en: "
		self.log_signal.emit(msg, green)
		self.log_signal.emit(str(round(time.time() - start_time,3)), green)
		if VERBOSE:
			print(time.time() - start_time)
		self.end_signal.emit("----------- FIN -----------")



################################################################################
################################################################################


	def stop(self):
     
		self.stop_event.set()
		parent = psutil.Process(self.parser_handler.pid)
		childs = parent.children(recursive=True)
		for i,child in enumerate(childs):
			self.quick_info_signal.emit(f"Veuillez attendre la fin des threads.\nIl reste {len(childs)-i} thread", purple)
			child.kill()
		parent.kill()
		with open("../pickle/organism_df", 'wb') as f:
			pickle.dump(self.organism_df, f)

		try: os.remove(bdd_path)
		except: pass
		self.quick_info_signal.emit("Threads Terminées.", green)
		msg = f"Thread {str(self.index)} terminée"
		self.log_signal.emit(msg, white)
		self.fin_parsing_signal.emit()
		self.end_signal.emit("----------- FIN -----------")
		self.terminate()


########################################################################################################################
##############################################  Parsing Functions ######################################################
