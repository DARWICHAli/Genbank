#!/usr/bin/env python3

#from Bio import Entrez
import pandas as pd
#from Bio import SeqIO

from functools import partial
import shutil
from PyQt5 import QtWidgets
from pip import main
from ui import Ui_MainWindow
from downloader_thread import ThreadClass, kingdoms_choice, regions_choice
import os


class Genbank(QtWidgets.QMainWindow):
	def __init__(self):
		QtWidgets.QMainWindow.__init__(self)
		self.MainWindow = QtWidgets.QMainWindow()
		self.mainwindow = Ui_MainWindow()
		self.mainwindow.setupUi(self.MainWindow)
		self.thread={}
		self.mainwindow.connect_ui(self)
		

################################################################################
################################################################################

	def get_kingdom_choice(self):
		selected_kingdoms = []
		if(self.mainwindow.checkBox_prokaryota.isChecked()):
			selected_kingdoms = selected_kingdoms + ["prokaryota"]
		if(self.mainwindow.checkBox_archaea.isChecked()):
			selected_kingdoms = selected_kingdoms + ["archaea"]
		if(self.mainwindow.checkBox_bacteria.isChecked()):
			selected_kingdoms = selected_kingdoms + ["bacteria"]
		if(self.mainwindow.checkBox_eukaryota.isChecked()):
			selected_kingdoms = selected_kingdoms + ["eukaryota"]
		if(self.mainwindow.inputKingdom.toPlainText() != ""):
			selected_kingdoms = selected_kingdoms + [self.mainwindow.inputKingdom.toPlainText()]
		return selected_kingdoms


################################################################################
################################################################################


	def get_region_choice(self):
		selected_regions = []
		if(self.mainwindow.checkBox_rrna.isChecked()):
			selected_regions = selected_regions + ["rrna"]
		if(self.mainwindow.checkBox_cds.isChecked()):
			selected_regions = selected_regions + ["cds"]
		if(self.mainwindow.checkBox_trna.isChecked()):
			selected_regions = selected_regions + ["trna"]
		if(self.mainwindow.checkBox_centromere.isChecked()):
			selected_regions = selected_regions + ["centromere"]
		if(self.mainwindow.checkBox_telomere.isChecked()):
			selected_regions = selected_regions + ["telomere"]
		if(self.mainwindow.checkBox_3utr.isChecked()):
			selected_regions = selected_regions + ["3utr"]
		if(self.mainwindow.checkBox_5utr.isChecked()):
			selected_regions = selected_regions + ["5utr"]
		if(self.mainwindow.checkBox_mobile_element.isChecked()):
			selected_regions = selected_regions + ["mobile element"]
		if(self.mainwindow.checkBox_mobile_ncrna.isChecked()):
			selected_regions = selected_regions + ["ncrna"]
		if(self.mainwindow.checkBox_mobile_intron.isChecked()):
			selected_regions = selected_regions + ["intron"]
		if(self.mainwindow.inputRegion.toPlainText() != ""):
			selected_regions = selected_regions + [self.mainwindow.inputRegion.toPlainText()]
		return selected_regions


################################################################################
################################################################################


	def start(self, msg):
		index = self.sender().index
		if(index == 1):
			self.log(str(msg))

	def get_result(self, organism_df):
		self.log("Organism dataframe received")
		self.organism_df = organism_df

################################################################################
################################################################################

	def pause(self, str):
		index = self.sender().index
		if(index == 1):
			print("pause")

################################################################################
################################################################################


	def reset(self):
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
			shutil.rmtree('./Results/.')
		except: print("cannot delete Results/.")
		self.mainwindow.logOutput.clear()
		self.mainwindow.progressBar.setValue(0)
		print("reset")


################################################################################
################################################################################

	# Logger
	def log(self, str):
		self.mainwindow.logOutput.insertPlainText(str + '\n')
		sb =self.mainwindow.logOutput.verticalScrollBar()
		sb.setValue(sb.maximum())
		
	def update_progress_bar(self, value):
		index = self.sender().index
		if index == 1:
			self.mainwindow.progressBar.setProperty("value", self.mainwindow.progressBar.value() + value)


	def worker(self):
		kingdoms_choice = self.get_kingdom_choice()
		regions_choice = self.get_region_choice()
		self.log(str(kingdoms_choice))
		self.log(str(regions_choice))
		self.thread[1] = ThreadClass(parent = None, index=1)
		self.thread[1].start()
		self.thread[1].any_signal.connect(self.start)
		self.thread[1].dataframe_result.connect(self.get_result)
		self.thread[1].progress_signal.connect(self.update_progress_bar)
		self.thread[1].time_signal.connect(self.start)
		


	def stop_worker(self):
		self.thread[1].stop()
		self.thread[1].stop_signal.connect(self.pause)
