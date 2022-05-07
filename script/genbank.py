from importlib.resources import path
import shutil
from threading import Thread
from PyQt5 import QtWidgets, QtCore, QtGui
from ui import Ui_MainWindow
from downloader_thread import DownloaderThread
import os
import pandas as pd
from parser_thread import ParserThread

class Genbank(QtWidgets.QMainWindow, QtCore.QObject):

	region_signal = QtCore.pyqtSignal(list)
	path_signal = QtCore.pyqtSignal(str)
	
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
		self.mainwindow.connect_ui(self)
		self.region_choice = []
		self.path_choice = []

		self.thread[1] = DownloaderThread(parent = self, index=1)
		self.thread[1].log_signal.connect(self.log)
		self.thread[1].dataframe_result.connect(self.get_result)
		self.thread[1].end_signal.connect(self.end)
		self.thread[1].start()
		

################################################################################
################################################################################

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

		print("Parsing Started...")
		parsing_choice = " >> ".join(self.path_choice.split('/')[2:])

		print("Selected organisms to parse: " + parsing_choice)

		self.log_color(255,255,255,255)
		self.clear_log()
		self.log("Selected organisms to parse: \n	" + parsing_choice)
		self.log_color(0, 250, 125, 255)

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
			self.log("Organism dataframe received")
			self.organism_df = organism_df
		


################################################################################
################################################################################

	def end(self, msg):

		self.mutex_end.lock()
		index = self.sender().index
		if(index == 1):
			self.log(str(msg))
			self.mainwindow.buttonStart.setEnabled(True)
			self.mainwindow.buttonStart.setText("Start Parsing")
			self.mainwindow.buttonStart.setStyleSheet("background-color: rgb(0, 250, 125);\n" "color:rgb(0, 4, 38);")
			self.thread[1].stop()
			self.isRunning = False
		elif(index == 2):
			self.log(str(msg))
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


	def reset(self):

		if os.getcwd().endswith("script"):
				os.chdir('../')
		try:
			if os.getcwd().endswith("GENOME_REPORTS"):
				os.chdir('../')
			shutil.rmtree('./GENOME_REPORTS')
		except: print("cannot delete GENOME_REPORTS")
		try:
			if os.getcwd().endswith("pickle"):
				os.chdir('../')
			shutil.rmtree('./pickle')
		except: print("cannot delete pickle")
		try:
			if os.getcwd().endswith("Results"):
				os.chdir('../')
			shutil.rmtree('./Results')
		except: print("cannot delete ../Results")

		self.mainwindow.logOutput.clear()
		self.mainwindow.progressBar.setValue(0)


################################################################################
################################################################################

	# Logger
	def log(self, str):
		self.mutex_log.lock()
		self.mainwindow.logOutput.insertPlainText(str + '\n')
		sb =self.mainwindow.logOutput.verticalScrollBar()
		sb.setValue(sb.maximum())
		self.mutex_log.unlock()
		
	def update_progress_bar(self, value):
		self.mutex_progress.lock()
		self.mainwindow.progressBar.setProperty("value", self.mainwindow.progressBar.value() + 100/value)
		self.mainwindow.progressBar.setFormat(str(int(self.mainwindow.progressBar.value()*value/100)) + " / " +str(value) + " NC")
		self.mutex_progress.unlock()

	def worker(self):

		if( self.isRunning == False):
			self.mainwindow.progressBar.setValue(0)
			self.mainwindow.progressBar.setFormat("")
			self.get_path_choice()
			self.get_region_choice()

			if not len(self.region_choice):
				self.log_color(255,0,0,255)
				self.log("Erreur User: Il faut choisir au moins une région fonctionnelle!")
				self.log_color(0, 250, 125, 255)
				return

			self.isRunning = True
			self.thread[2] = ParserThread(parent = self, index=2, organism_df = self.organism_df, path_choice=self.path_choice, regions_choice=self.region_choice)
			self.thread[2].end_signal.connect(self.end)
			self.thread[2].log_signal.connect(self.log)
			self.thread[2].progress_signal.connect(self.update_progress_bar)
			self.thread[2].start()

			self.mainwindow.buttonStart.setText("Stop Parsing")
			self.mainwindow.buttonStart.setStyleSheet("background-color: rgb(100, 20, 15);\n" "color:rgb(255, 255, 255);")
			self.mainwindow.logOutput.insertPlainText('Program Started\n')

		else:
			self.mainwindow.buttonStart.setText("Start Parsing")
			self.mainwindow.buttonStart.setStyleSheet("background-color: rgb(0, 250, 125);\n" "color:rgb(0, 4, 38);")
			self.isRunning = False
			self.clear_log()
			self.mainwindow.logOutput.insertPlainText('Program stopped\n')
			self.thread[2].stop()
			self.mainwindow.progressBar.setValue(0)
			self.mainwindow.progressBar.setFormat("0/0")
		
	def log_color(self, r, g, b, a):
		self.mainwindow.logOutput.setTextColor(QtGui.QColor(r, g, b, a))

	def clear_log(self):
		self.mainwindow.logOutput.clear()