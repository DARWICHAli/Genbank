from ui import *
from PyQt5.QtWidgets import QApplication

if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    ui.connect_ui()
    MainWindow.show()
    sys.exit(app.exec_())
 