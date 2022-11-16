from PyQt5.QtWidgets import QApplication
import genbank as genbank
from PyQt5 import QtGui
from multiprocessing import freeze_support

def main():
    import sys
    app = QApplication(sys.argv)
    pal = app.palette()
    pal.setColor(QtGui.QPalette.Window, QtGui.QColor(0, 4, 38,255))
    app.setPalette(pal)
    GenBank = genbank.Genbank()

    GenBank.MainWindow.show()
    sys.exit(app.exec_())
 
if __name__ == "__main__":
    freeze_support()
    main()

