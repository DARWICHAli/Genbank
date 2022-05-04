from PyQt5.QtWidgets import QApplication
import genbank as genbank

def main():
    import sys
    app = QApplication(sys.argv)

    GenBank = genbank.Genbank()

    GenBank.MainWindow.show()
    sys.exit(app.exec_())
 
if __name__ == "__main__":
    main()

