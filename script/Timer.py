import time
from PyQt5 import QtCore


class Timer(QtCore.QThread):

    time_signal = QtCore.pyqtSignal(str)

    def __init__(self, parent, time=0):
        super(Timer, self).__init__(parent)
        self.isRunning = True
        self.isPaused = False
        self.parent = parent
        self.time = time
        self.hours = int (self.time / 3600)
        self.minutes = int ((self.time % 3600) / 60)
        self.seconds = int ((self.time % 3600) % 60)
        self.time_str = str(self.hours) + ":" + str(self.minutes) + ":" + str(self.seconds)
        self.parent.pause_signal.connect(self.get_pause)

    def run(self):
        while(True):
            if(self.isRunning and not self.isPaused):
                time.sleep(1)
                self.time += 1
                self.hours = int(self.time / 3600)
                self.minutes = int ((self.time % 3600) / 60)
                self.seconds = int ((self.time % 3600) % 60)
                str_hours = str(self.hours) if int(self.hours / 10) > 0 else "0" + str(self.hours) 
                str_minutes = str(self.minutes) if int(self.minutes / 10) > 0 else "0" + str(self.minutes) 
                str_seconds = str(self.seconds) if int(self.seconds / 10) > 0 else "0" + str(self.seconds) 
                self.time_str = "   Temps écoulé: " + str_hours + ":" + str_minutes + ":" + str_seconds
                self.time_signal.emit(self.time_str)
            elif(not self.isRunning):
                break
            

    def get_pause(self, val):
        self.isPaused = val

    def stop(self):
        self.terminate()