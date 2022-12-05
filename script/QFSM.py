import typing
from PyQt5 import QtCore
from PyQt5.QtWidgets import QFileSystemModel
from PyQt5.QtGui import QBrush, QIcon, QPixmap

class QFSM(QFileSystemModel):

    data_changed = QtCore.pyqtSignal(QtCore.QModelIndex, QtCore.QModelIndex)

    def __init__(self):
        super().__init__()
        self.indexList = QtCore.QPersistentModelIndex

    def data(self, index: QtCore.QModelIndex, role: int = ... ) -> QtCore.QVariant:
        try:
            if(role == QtCore.Qt.DecorationRole):
                return QtCore.QVariant(QIcon(QPixmap("../images/white_folder.png")))
            return super().data(index, role)
        except:
            return super().data(index, role)

    def setData(self, index: QtCore.QModelIndex, value: typing.Any, role: int = ...) -> bool:
        try:
            if (role == QtCore.Qt.DecorationRole):

                if (value == "../images/green_icon.png"):
                    self.indexList.insert(index)
                else:
                    self.indexList.remove(index)
                self.data_changed.emit(index, index, {QtCore.Qt.DecorationRole})
                return True
            return super().setData(index, value, role)
        except:
            return super().setData(index, value, role)
