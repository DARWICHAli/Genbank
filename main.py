from PyQt5.QtWidgets import QApplication
import asyncio
import genbank as genbank

def main():
    import sys
    app = QApplication(sys.argv)

    GenBank = genbank.Genbank()

    GenBank.MainWindow.show()
    sys.exit(app.exec_())
 
asyncio.run(main())