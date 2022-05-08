from importlib.resources import path
import shutil
from threading import Thread
from PyQt5 import QtWidgets, QtCore, QtGui
from ui import Ui_MainWindow
from downloader_thread import DownloaderThread
import os
import pandas as pd
from parser_thread import ParserThread
import time
from Timer import Timer

green = [0,255,0,255]
white = [255,255,255,255]
purple = [255,0,255,255] 
red = [255,0,0,255] 

class Genbank(QtWidgets.QMainWindow, QtCore.QObject):

	region_signal = QtCore.pyqtSignal(list)
	path_signal = QtCore.pyqtSignal(str)
	pause_signal = QtCore.pyqtSignal(bool)
	
	def __init__(self, parent = None, index = 0):
		super(Genbank, self).__init__(parent)
		QtWidgets.QMainWindow.__init__(self)
		self.isRunning = True
		self.parsing = False
		self.index = index
		self.MainWindow = QtWidgets.QMainWindow()
		self.mainwindow = Ui_MainWindow()
		self.mainwindow.setupUi(self.MainWindow)
		self.thread={}
		self.mutex_log = QtCore.QMutex()
		self.mutex_progress = QtCore.QMutex()
		self.mutex_end = QtCore.QMutex()
		self.mutex_info = QtCore.QMutex()
		self.mainwindow.connect_ui(self)
		self.region_choice = []
		self.path_choice = ""
		self.nb_NC = 0
		self.nb_parsed = 0
		self.nbThreads = 0
		self.thread[3] = Timer(parent = self)
		self.thread[3].time_signal.connect(self.update_time)
		self.thread[3].start()

		self.thread[1] = DownloaderThread(parent = self, index=1)
		self.thread[1].log_signal.connect(self.log)
		self.thread[1].dataframe_result.connect(self.get_result)
		self.thread[1].end_signal.connect(self.end)
		self.thread[1].start()
		self.quickInfo("Téléchargement de l'arborescence, veuillez patienter.\n Cela ne prendra pas plus de 3 minutes.", white)

################################################################################
################################################################################

	def update_time(self, time):
		self.mainwindow.timer.setText(time)

	def get_path_choice(self):
		self.path_choice = []

		if(self.isRunning):
			return

		try:
			index = self.mainwindow.treeView.selectedIndexes()[0]
			path = self.mainwindow.model.filePath(index)
			path = path.split('/')[::-1]

			for f in path:
				if(f == "Results"): break
				self.path_choice.append(f)
			self.path_choice.append("Results")
			self.path_choice.append("..")
			self.path_choice = self.path_choice[::-1]
			self.path_choice = "/".join(self.path_choice)
		except:
			self.path_choice = "../Results"

		parsing_choice = (" >> ".join(self.path_choice.split('/')[2:])).upper()
		
		self.clear_log()
		self.log("Selected organisms to parse:", white)
		self.log("	" + parsing_choice, green)
		parsing_choice = self.path_choice.split('/')[2:]

		if(len(parsing_choice) == 4 and len(parsing_choice[len(parsing_choice)-1]) > 30):
			parsing_choice[len(parsing_choice)-1] = " ".join(parsing_choice[len(parsing_choice)-1].split("_")[0:2])

		parsing_choice = "\n >> ".join(parsing_choice).upper().replace("_"," ")

		self.quickInfo("Organismes selectionnés:\n" + parsing_choice, green )

################################################################################
################################################################################


	def get_region_choice(self):
		selected_regions = []
		if(self.mainwindow.checkBox_rrna.isChecked()):
			selected_regions = selected_regions + ["rRNA"]
		if(self.mainwindow.checkBox_cds.isChecked()):
			selected_regions = selected_regions + ["CDS"]
		if(self.mainwindow.checkBox_trna.isChecked()):
			selected_regions = selected_regions + ["tRNA"]
		if(self.mainwindow.checkBox_centromere.isChecked()):
			selected_regions = selected_regions + ["centromere"]
		if(self.mainwindow.checkBox_telomere.isChecked()):
			selected_regions = selected_regions + ["telomere"]
		if(self.mainwindow.checkBox_3utr.isChecked()):
			selected_regions = selected_regions + ["3'UTR"]
		if(self.mainwindow.checkBox_5utr.isChecked()):
			selected_regions = selected_regions + ["5'UTR"]
		if(self.mainwindow.checkBox_mobile_element.isChecked()):
			selected_regions = selected_regions + ["mobile_element"]
		if(self.mainwindow.checkBox_mobile_ncrna.isChecked()):
			selected_regions = selected_regions + ["ncRNA"]
		if(self.mainwindow.checkBox_mobile_intron.isChecked()):
			selected_regions = selected_regions + ["intron"]
		if(len(self.mainwindow.inputRegion.toPlainText())):
			selected_regions = selected_regions + [self.mainwindow.inputRegion.toPlainText()]

		self.region_choice = selected_regions



################################################################################
################################################################################

	def get_result(self, organism_df):
		index = self.sender().index 
		if index == 1:
			self.organism_df = organism_df
		


################################################################################
################################################################################

	def end(self, msg):

		self.mutex_end.lock()
		index = self.sender().index
		self.pause_signal.emit(True)

		if(index == 1):
			self.log(str(msg), green)
			self.mainwindow.buttonStart.setEnabled(True)
			self.mainwindow.buttonStart.setText("Start Parsing")
			self.mainwindow.buttonStart.setStyleSheet("background-color: rgb(0, 250, 125);\n" "color:rgb(0, 4, 38);")
			self.thread[1].stop()
			self.isRunning = False
			self.quickInfo("Téléchargement de l'arborescence terminée.\nVeuillez selectionner les organismes et régions à parser.\nSi pas de selection d'organismes, tous seront parsés.\n\nVous pouvez arrêter le parsing à tout moment.\nIl faudra attendre la fin des threads.", white)
			
		elif(index == 2):
			self.log(str(msg), green)
			self.mainwindow.buttonStart.setEnabled(True)
			self.mainwindow.buttonStart.setText("Start Parsing")
			self.mainwindow.buttonStart.setStyleSheet("background-color: rgb(0, 250, 125);\n" "color:rgb(0, 4, 38);")
			self.thread[2].stop()
			self.isRunning = False
			self.mainwindow.progressBar.setValue(0)
			self.mainwindow.progressBar.setFormat("0/0")
			
		self.mutex_end.unlock()

################################################################################
################################################################################

	# Logger
	def log(self, str, color):
		self.mutex_log.lock()
		self.log_color(color)
		self.mainwindow.logOutput.insertPlainText(str + '\n')
		sb =self.mainwindow.logOutput.verticalScrollBar()
		sb.setValue(sb.maximum())
		self.mutex_log.unlock()
		
	def quickInfo(self, str, color):
		self.mutex_info.lock()
		self.quick_info_color(color)
		self.mainwindow.quickInfo.setText(str)
		self.mutex_info.unlock()

	def get_nb_thread(self, val):
		self.nbThreads = val
	
	def update_progress_bar(self, value):
		self.mutex_progress.lock()
		if value > 0:
			self.nb_NC = value
		else:
			self.nb_parsed += 1
			self.mainwindow.progressBar.setProperty("value", (float(self.nb_parsed)/float(self.nb_NC))*100. )
		self.mainwindow.progressBar.setFormat(str(self.nb_parsed) + " / " +str(self.nb_NC) + " NC")
		self.mutex_progress.unlock()

	def worker(self):

		if( self.isRunning == False):
			self.mainwindow.progressBar.setValue(0)
			self.mainwindow.progressBar.setFormat("")
			self.get_path_choice()
			self.get_region_choice()
			self.pause_signal.emit(False)
			if not len(self.region_choice):
				
				self.log("Erreur User: Il faut choisir au moins une région fonctionnelle!", [255,0,0,255])
				self.quickInfo("Veuillez choisir au moins une région fonctionnelle.", red)
				return

			if(self.path_choice == "../Results"):
				self.log("Parsing de toutes les kingdoms...", green)
				self.quickInfo("Parsing de toutes les kingdoms en cours...", green)
			else:
				parsing_choice = ("\n >> ".join(self.path_choice.split('/')[2:])).upper()
				self.quickInfo("Parsing de " + parsing_choice + " en cours...", green)

			self.isRunning = True
			self.thread[2] = ParserThread(parent = self, index=2, organism_df = self.organism_df, path_choice=self.path_choice, regions_choice=self.region_choice)
			self.thread[2].end_signal.connect(self.end)
			self.thread[2].log_signal.connect(self.log)
			self.thread[2].progress_signal.connect(self.update_progress_bar)
			self.thread[2].quick_info_signal.connect(self.quickInfo)
			self.thread[2].start()

			self.mainwindow.buttonStart.setText("Stop Parsing")
			self.mainwindow.buttonStart.setStyleSheet("background-color: rgb(100, 20, 15);\n" "color:rgb(255, 255, 255);")
			self.mainwindow.logOutput.insertPlainText('Program Started\n')

		else:
			self.mainwindow.buttonStart.setText("Start Parsing")
			self.mainwindow.buttonStart.setStyleSheet("background-color: rgb(0, 250, 125);\n" "color:rgb(0, 4, 38);")
			self.isRunning = False
			self.thread[2].stop()
			self.clear_log()
			self.mainwindow.logOutput.insertPlainText('------- Program stopped ------\n')
			self.mainwindow.progressBar.setValue(0)
			self.mainwindow.progressBar.setFormat("0/0")
		
	def log_color(self, color):
		try:
			self.mainwindow.logOutput.setTextColor(QtGui.QColor(color[0], color[1], color[2], color[3]))
		except: pass

	def quick_info_color(self, color):
		color_str = "color:rgb(" + str(color[0]) + "," + str(color[1]) + "," + str(color[2]) + ")"
		style = "font: 8.5pt Futura;\n background-color: rgb(0, 4, 38);\n" + color_str + ";"
		self.mainwindow.quickInfo.setStyleSheet( style )


	def clear_log(self):
		self.mainwindow.logOutput.clear()