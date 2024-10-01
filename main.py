

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import (QMainWindow, QApplication, QWidget,QHeaderView,
    QPushButton, QAction, QLineEdit, QMessageBox, QTableWidgetItem, QAbstractItemView)
from PyQt5.QtWidgets import QDialog, QApplication, QTableWidgetItem
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QTableWidget, QStyledItemDelegate

import sys
def legendre_symbol(a, p):
    ls = pow(a, (p - 1) // 2, p)
    return -1 if ls == p - 1 else ls

def modular_sqrt(a, p):
    if legendre_symbol(a, p) != 1:
        return 0
    elif a == 0:
        return 0
    elif p == 2:
        return p
    elif p % 4 == 3:
        return pow(a, (p + 1) // 4, p)

    s = p - 1
    e = 0
    while s % 2 == 0:
        s //= 2
        e += 1

    n = 2
    while legendre_symbol(n, p) != -1:
        n += 1

    x = pow(a, (s + 1) // 2, p)
    b = pow(a, s, p)
    g = pow(n, s, p)
    r = e

    while True:
        t = b
        m = 0
        for m in range(r):
            if t == 1:
                break
            t = pow(t, 2, p)

        if m == 0:
            return x

        gs = pow(g, 2 ** (r - m - 1), p)
        g = (gs * gs) % p
        x = (x * gs) % p
        b = (b * g) % p
        r = m


def summa(a, b):  # функция нахождения суммы по модулю
    sum = (a + b) % p
    if (sum < 0):
        sum += p
        return sum
    else:
        return sum


def proiz(a, b):  # функция нахождения произведения по модулю
    sum = (a * b) % p
    if (sum < 0):
        sum += p
        return sum
    else:
        return sum


def step(a, b):  # функция возведения в степень числа !mod|a по простому модулю mod
    if (b < 0):
        z = b
        while z < 0:
            z = z + (p - 1)
        return (a ** z) % p
    else:
        return (a ** b) % p
def sum_tochek(A,x1,y1,x2,y2):
    if (x1 == x2) and (y1 == y2):
        x3 = int(summa(step(proiz(summa(proiz(3, step(x1, 2)), A), step(proiz(2, y1), (-1))), 2), proiz(-2, x1)))


        y3 = summa(-y1, proiz(proiz(summa(proiz(step(x1, 2), 3), A), step(proiz(2, y1), -1)), summa(x1, -x3)))
        return x3,y3
    else:
        x3 = summa(summa(step(proiz(summa(y2, -y1), step(summa(x2, -x1), -1)), 2), -x1), -x2)

        y3 = summa(proiz(summa(x1, -x3), proiz(summa(y2, -y1), step(summa(x2, -x1), -1))), -y1)
        return x3,y3

def por(a, b, mod, x1, y1):
    x2 = x1
    y2 = y1
    por = 1
    if y1 == 0:
        por = 1

    else:

        while por != -1:
            if (x2 == x1) and ((y2 - mod) == ((y1) * (-1))):
                break

            x2,y2=sum_tochek(a,x1,y1,x2,y2)
            por = por + 1
    return por + 1


def vse_t(a, b, p):
    q = [[0], 1]

    for i in range(p):

        s = summa(summa(step(i, 3), proiz(a, i)), b)
        e = legendre_symbol(s, p)
        if e == -1:

            continue
        elif e == 1:

            c = modular_sqrt(s, p)
            if c == 0:
                q.append([i, c])

                q.append(por(a, b, p, i, c))
            else:
                q.append([i, c])

                q7 = por(a, b, p, i, c)
                q.append(q7)
                q.append([i, c * (-1)])

                q.append(q7)
        elif e == 0:

            # Символ вычислять не нужно сразу записываем sj
            q.append([i, s])

            q.append(por(a, b, p, i, s))
    return q
def str_group(a,b,p):
    q = vse_t(a, b, p)
    n = int(len(q) / 2)

    t7 = 3
    max_p = q[t7 - 2]
    max_t = q[0]

    while t7 <= len(q):

        if q[t7] > max_p:
            max_p = q[t7]
            max_t = q[t7 - 1]
        t7 = t7 + 2

    if max_p == n:

        return "E(Fp)=Z", max_p, "=<", max_t[0], ":", max_t[1], ">"

    else:
        m = int(n / max_p)
        q2 = [max_t]
        x1 = max_t[0]
        y1 = max_t[1]
        x2 = x1
        y2 = y1
        for i in range(max_p - 1):
            x2, y2 = sum_tochek(a, x1, y1, x2, y2)
            q2.append([x2, y2])


        i2 = 0

        while i2 <= n:

            for i in range(max_p - 1):
                if q[i2] == q2[i]:
                    del q[i2:(i2 + 2)]
                    continue
            i2 = i2 + 2
        i = 1

        while i <= len(q) - 1:
            if q[i] == m:
                max_t2 = q[i - 1]
                break
            i = i + 2
        if len(max_t2) == 1 and max_t2[0] == 0:

            return "E(Fp)=Z", max_p, "+Z", m, "=<", max_t[0], ":", max_t[1], ">+<", max_t2[0], ">"

        return "E(Fp)=Z", max_p, "+Z", m, "=<", max_t[0], ":", max_t[1], ">+<", max_t2[0], ":",max_t2[1], ">"
def izomorf_kriv(a,b,p):
    count = (p - 1) // 2
    mas1=[]
    for i in range(1,count+1):
        mas1.append([proiz(a,step(i,4)),proiz(step(i,6),b)])
    return mas1



class ReadOnlyDelegate(QStyledItemDelegate):
    def createEditor(self, parent, option, index):
        
        return
def remove_commas(string):
    trans_table = {ord("'"): None,ord(","): None,ord(" "): None}
    return string.translate(trans_table)

def insert_row_in_the_end_of_table(table):
    row_index_to_isert = table.rowCount() + 1
    table.setRowCount(row_index_to_isert)

    for column in range(table.columnCount()):
        item = QTableWidgetItem('')
        item.setFlags(item.flags() & ~Qt.ItemIsEditable)
        table.setItem(row_index_to_isert - 1, column, item)
class Ui_MainWindow(object):
    def summa(self,a, b):  # функция нахождения суммы по модулю
        sum = (a + b) % p
        if (sum < 0):
            sum += p
            return sum
        else:
            return sum

    def proiz(self, a, b):  # функция нахождения произведения по модулю
        sum = (a * b) % p
        if (sum < 0):
            sum += p
            return sum
        else:
            return sum

    def step(self, a, b):  # функция возведения в степень числа !mod|a по простому модулю mod
        if (b < 0):
            z = b
            while z < 0:
                z = z + (p - 1)
            return (a ** z) % p
        else:
            return (a ** b) % p

    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 600)
        MainWindow.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        MainWindow.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.p_label = QtWidgets.QLabel(self.centralwidget)
        self.p_label.setGeometry(QtCore.QRect(0, 0, 300, 40))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.p_label.setFont(font)
        self.p_label.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.p_label.setTextFormat(QtCore.Qt.AutoText)
        self.p_label.setAlignment(QtCore.Qt.AlignCenter)
        self.p_label.setObjectName("p_label")
        self.p_btn = QtWidgets.QPushButton(self.centralwidget)
        self.p_btn.setGeometry(QtCore.QRect(370, 0, 80, 40))
        self.p_btn.setStyleSheet("font: 16pt \"MS Shell Dlg 2\";\n"
"background-color: rgb(223, 223, 223);")
        self.p_btn.setObjectName("p_btn")
        self.p_line = QtWidgets.QLineEdit(self.centralwidget)
        self.p_line.setGeometry(QtCore.QRect(300, 0, 70, 40))
        self.p_line.setStyleSheet("font: 14pt \"MS Shell Dlg 2\";")
        self.p_line.setObjectName("p_line")
        self.pred_table = QtWidgets.QTableWidget(self.centralwidget)
        self.pred_table.setGeometry(QtCore.QRect(0, 80, 500, 520))
        self.pred_table.setObjectName("pred_table")
        self.pred_table.setColumnCount(3)
        self.pred_table.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setStyleStrategy(QtGui.QFont.PreferAntialias)
        item.setFont(font)
        self.pred_table.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        font = QtGui.QFont()
        font.setPointSize(14)
        item.setFont(font)
        self.pred_table.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        font = QtGui.QFont()
        font.setPointSize(14)
        item.setFont(font)
        self.pred_table.setHorizontalHeaderItem(2, item)
        self.pred_table.horizontalHeader().setVisible(True)
        self.pred_table.horizontalHeader().setCascadingSectionResizes(False)
        self.pred_table.horizontalHeader().setDefaultSectionSize(166)
        self.pred_table.horizontalHeader().setHighlightSections(False)
        self.pred_table.horizontalHeader().setMinimumSectionSize(25)
        self.pred_table.horizontalHeader().setSortIndicatorShown(False)
        self.pred_table.horizontalHeader().setStretchLastSection(False)
        header = self.pred_table.horizontalHeader()
        header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        self.pred_label = QtWidgets.QLabel(self.centralwidget)
        self.pred_label.setGeometry(QtCore.QRect(0, 40, 500, 40))
        self.pred_label.setStyleSheet("font: 18pt \"MS Shell Dlg 2\";")
        self.pred_label.setAlignment(QtCore.Qt.AlignCenter)
        self.pred_label.setObjectName("pred_label")

        self.izomorf_table = QtWidgets.QTableWidget(self.centralwidget)
        self.izomorf_table.setGeometry(QtCore.QRect(500, 80, 300, 520))
        self.izomorf_table.setObjectName("izomorf_table")
        self.izomorf_table.setColumnCount(1)
        self.izomorf_table.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignLeft)
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setStyleStrategy(QtGui.QFont.PreferAntialias)
        item.setFont(font)
        self.izomorf_table.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        font = QtGui.QFont()
        font.setPointSize(14)
        item.setFont(font)
        self.izomorf_table.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        font = QtGui.QFont()
        font.setPointSize(14)
        item.setFont(font)
        self.izomorf_table.setHorizontalHeaderItem(2, item)
        self.izomorf_table.horizontalHeader().setVisible(True)
        self.izomorf_table.horizontalHeader().setCascadingSectionResizes(False)
        self.izomorf_table.horizontalHeader().setDefaultSectionSize(166)
        self.izomorf_table.horizontalHeader().setHighlightSections(False)
        self.izomorf_table.horizontalHeader().setMinimumSectionSize(25)
        self.izomorf_table.horizontalHeader().setSortIndicatorShown(False)
        self.izomorf_table.horizontalHeader().setStretchLastSection(False)
        header = self.izomorf_table.horizontalHeader()
        header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)



        self.izomorf_btn = QtWidgets.QPushButton(self.centralwidget)
        self.izomorf_btn.setGeometry(QtCore.QRect(450, 0, 350, 40))
        self.izomorf_btn.setStyleSheet("font: 16pt \"MS Shell Dlg 2\";\n"
"background-color: rgb(223, 223, 223);")

        self.izomorf_btn.setObjectName("izomorf_btn")
        self.str_group_btn = QtWidgets.QPushButton(self.centralwidget)
        self.str_group_btn.setGeometry(QtCore.QRect(500, 40, 300, 40))
        self.str_group_btn.setStyleSheet("font: 12pt \"MS Shell Dlg 2\";\n"
"background-color: rgb(223, 223, 223);")
        self.str_group_btn.setObjectName("str_group_btn")
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        delegate = ReadOnlyDelegate(self.pred_table)
        self.pred_table.setItemDelegateForRow(1, delegate)
        self.pred_table.setItemDelegateForColumn(0, delegate)
        self.pred_table.setItemDelegateForColumn(1, delegate)
        self.pred_table.setItemDelegateForColumn(2, delegate)
        delegate = ReadOnlyDelegate(self.izomorf_table)
        self.izomorf_table.setItemDelegateForColumn(0, delegate)





    def on_click(self):
        row = self.pred_table.rowCount()
        
        if row!=0:
            for i in range(row+1):
                self.pred_table.removeRow(0)
        row = self.izomorf_table.rowCount()
        
        if row != 0:
            for i in range(row + 1):
                self.izomorf_table.removeRow(0)

        textboxValue = self.p_line.text()
        if textboxValue.isdigit()!=True:
            infoBox = QMessageBox()
            infoBox.setText("Ошибка")
            infoBox.setInformativeText("Введите число")
            infoBox.setWindowTitle("Внимание")
            infoBox.setStandardButtons(QMessageBox.Ok)
            infoBox.setIcon(QMessageBox.Warning)
            infoBox.setWindowIcon(QtGui.QIcon('icon.png'))
            infoBox.exec_()
            self.p_line.setText("")
            return 0

        if int(textboxValue)<=3:
            infoBox = QMessageBox()
            infoBox.setText("Ошибка")
            infoBox.setInformativeText("Введите модуль>3")
            infoBox.setWindowTitle("Внимание")
            infoBox.setStandardButtons(QMessageBox.Ok)
            infoBox.setIcon(QMessageBox.Warning)
            infoBox.setWindowIcon(QtGui.QIcon('icon.png'))
            infoBox.exec_()
            self.p_line.setText("")
            return 0
        temp3 = 0
        global pervobr
        with open("первообразные_корни.txt") as file:
            pervobr=''
            for item in file:
                temp = str("")
                stroka = item

                for j in range(len(stroka)):
                    stroka = list(stroka)

                    if stroka[j] == " ":
                        break
                    temp = temp + (stroka[j])
                    temp2 = j

                if int(temp) == int(textboxValue):
                    for i in range(temp2 + 1, len(stroka)):
                        pervobr = pervobr + stroka[i]
                        temp3 = 1
                    break


        zero="0"
        if temp3 != 1:
            infoBox = QMessageBox()
            infoBox.setText("Ошибка")
            infoBox.setInformativeText("Введите простой модуль <1000")
            infoBox.setWindowTitle("Внимание")
            infoBox.setStandardButtons(QMessageBox.Ok)
            infoBox.setIcon(QMessageBox.Warning)
            infoBox.setWindowIcon(QtGui.QIcon('icon.png'))
            infoBox.exec_()
            self.p_line.setText("")
            return 0
        mas=[]
        j_inv = 0
        global p
        p = int(textboxValue)
        global j0
        j0=1728%p
        j0str=str(j0)
        pervobr = int(pervobr)
        if (p % 12 == 1) or (p % 12 == 7):
            
            for i in range(6):

                mas.append(0)
                mas.append(0)
                mas.append((pervobr**i)%p)
                temp=((pervobr**i)%p)
                temp2 = str(temp)

                temp2_2 = "E(Fp): Y^2=X^3+"+temp2
                self.zero = QtWidgets.QTableWidgetItem(zero)
                self.zero.setTextAlignment(QtCore.Qt.AlignHCenter)
                rowPosition = self.pred_table.rowCount()
                self.pred_table.insertRow(rowPosition)
                self.pred_table.setItem(rowPosition, 0, QtWidgets.QTableWidgetItem(self.zero))
                self.pred_table.setItem(rowPosition, 1, QtWidgets.QTableWidgetItem(temp2_2))



        if p % 12 == 5 or p % 12 == 11:
            
            for i in range(2):
                mas.append(0)
                mas.append(0)
                mas.append((pervobr ** i) % p)
                temp = ((pervobr ** i) % p)
                temp2 = str(temp)

                temp2_2 = "E(Fp): Y^2=X^3+" + temp2
                self.zero = QtWidgets.QTableWidgetItem(zero)
                self.zero.setTextAlignment(QtCore.Qt.AlignHCenter)
                rowPosition = self.pred_table.rowCount()
                self.pred_table.insertRow(rowPosition)
                self.pred_table.setItem(rowPosition, 0, QtWidgets.QTableWidgetItem(self.zero))
                self.pred_table.setItem(rowPosition, 1, QtWidgets.QTableWidgetItem(temp2_2))

        for i in range(0, p-1):
                j_inv += 1
                if j_inv == j0:
                    if (p % 12 == 1) or (p % 12 == 5):
                        
                        for i in range(4):
                            mas.append(j0)

                            mas.append((pervobr ** i) % p)
                            mas.append(0)
                            temp = ((pervobr ** i) % p)
                            temp2 = str(temp)

                            temp2_2 = "E(Fp): Y^2=X^3+" + temp2+"X"
                            self.j0str = QtWidgets.QTableWidgetItem(j0str)
                            self.j0str.setTextAlignment(QtCore.Qt.AlignHCenter)
                            rowPosition = self.pred_table.rowCount()
                            self.pred_table.insertRow(rowPosition)
                            self.pred_table.setItem(rowPosition, 0, QtWidgets.QTableWidgetItem(self.j0str))
                            self.pred_table.setItem(rowPosition, 1, QtWidgets.QTableWidgetItem(temp2_2))
                    if p % 12 == 7 or p % 12 == 11:
                        
                        for i in range(2):
                            mas.append(j0)

                            mas.append((pervobr ** i) % p)
                            mas.append(0)
                            temp = ((pervobr ** i) % p)
                            temp2 = str(temp)
                            self.j0str = QtWidgets.QTableWidgetItem(j0str)
                            self.j0str.setTextAlignment(QtCore.Qt.AlignHCenter)
                            temp2_2 = "E(Fp): Y^2=X^3+" + temp2+"X"
                            rowPosition = self.pred_table.rowCount()
                            self.pred_table.insertRow(rowPosition)
                            self.pred_table.setItem(rowPosition, 0, QtWidgets.QTableWidgetItem(self.j0str))
                            self.pred_table.setItem(rowPosition, 1, QtWidgets.QTableWidgetItem(temp2_2))
                    continue

                mas.append(j_inv)
                temp=proiz(proiz(step(4,-1),-27),proiz(j_inv,step(j_inv-j0,-1)))
                mas.append(temp)
                temp1=str(j_inv)
                temp2=str(temp)
                pervobr=int(pervobr)
                temp2_2="E(Fp): Y^2=X^3+"+(temp2)+"X+"+(temp2)
                temp3_1=(temp*(pervobr**2))%p
                temp3_2 = (temp * (pervobr ** 3)) % p
                mas.append(temp3_1)
                mas.append(temp3_2)
                temp3_11=str(temp3_1)
                temp3_22=str(temp3_2)
                temp3_3 = "E(Fp): Y^2=X^3+" + (temp3_11) + "X+" + (temp3_22)
                self.temp1 = QtWidgets.QTableWidgetItem(temp1)
                self.temp1.setTextAlignment(QtCore.Qt.AlignHCenter)
                rowPosition = self.pred_table.rowCount()
                self.pred_table.insertRow(rowPosition)
                self.pred_table.setItem(rowPosition, 0, QtWidgets.QTableWidgetItem(self.temp1))
                self.pred_table.setItem(rowPosition, 1, QtWidgets.QTableWidgetItem(temp2_2))
                self.pred_table.setItem(rowPosition, 2, QtWidgets.QTableWidgetItem(temp3_3))

        self.pred_table.setSelectionMode(QAbstractItemView.SingleSelection)
    def on_click_strgroup(self):

        row = self.pred_table.currentIndex().row()
        col = self.pred_table.currentIndex().column()


        temp=self.pred_table.model().data(self.pred_table.currentIndex())
        
        
        
        if (row and col)<0 or temp==None or col==0:
            infoBox = QMessageBox()
            infoBox.setText("Выберите эллиптическую кривую")
            infoBox.setInformativeText("")
            infoBox.setWindowTitle("Внимание")
            infoBox.setStandardButtons(QMessageBox.Ok)
            infoBox.setIcon(QMessageBox.Warning)
            infoBox.setWindowIcon(QtGui.QIcon('icon.png'))
            infoBox.exec_()

            return 0
        if temp.count("X")==2:
            if temp.count("+")==2:
                temp_pl=temp.find("+")
                temp_X = temp.find("X",temp_pl)
                temp_a=int(temp[temp_pl:temp_X])
                temp_b=int(temp[temp_X+1:])
                
                
            else:
                temp_pl = temp.find("+")
                temp_X = temp.find("X", temp_pl)
                temp_a = int(temp[temp_pl:temp_X])
                temp_b=0
                
        else:
            temp_pl = temp.find("+")

            temp_b = int(temp[temp_pl:])
            temp_a = 0
            

        temp_c=str_group(temp_a,temp_b,p)
        temp_c=str(temp_c)
        temp_c=remove_commas(temp_c)

        temp_c.lstrip(), temp_c.rstrip()
        infoBox = QMessageBox()
        

        infoBox.setText("Структура группы:")
        infoBox.setInformativeText(temp_c)
        infoBox.setWindowTitle("Внимание")
        infoBox.setStandardButtons(QMessageBox.Ok)
        infoBox.setIcon(QMessageBox.Information)
        infoBox.setWindowIcon(QtGui.QIcon('icon.png'))
        infoBox.exec_()

        return 0


    def on_click_izomorf(self):
        row = self.izomorf_table.rowCount()
        
        if row != 0:
            for i in range(row + 1):
                self.izomorf_table.removeRow(0)
        row = self.pred_table.currentIndex().row()
        col = self.pred_table.currentIndex().column()



        temp = self.pred_table.model().data(self.pred_table.currentIndex())
        
        
        
        
        if (row and col) < 0 or temp == None or col == 0:
            infoBox = QMessageBox()
            infoBox.setText("Выберите эллиптическую кривую")
            infoBox.setInformativeText("")
            infoBox.setWindowTitle("Внимание")
            infoBox.setStandardButtons(QMessageBox.Ok)
            infoBox.setIcon(QMessageBox.Warning)
            infoBox.setWindowIcon(QtGui.QIcon('icon.png'))
            infoBox.exec_()

            return 0
        if temp.count("X")==2:
            if temp.count("+")==2:
                temp_pl=temp.find("+")
                temp_X = temp.find("X",temp_pl)
                temp_a=int(temp[temp_pl:temp_X])
                temp_b=int(temp[temp_X+1:])
                
                
                
            else:
                temp_pl = temp.find("+")
                temp_X = temp.find("X", temp_pl)
                temp_a = int(temp[temp_pl:temp_X])
                temp_b=0
                
                
        else:
            temp_pl = temp.find("+")

            temp_b = int(temp[temp_pl:])
            temp_a = 0
            
            
        row = self.pred_table.currentIndex().row()
        col = self.pred_table.currentIndex().column()
        temp0 = self.pred_table.item(row, 0).text()
        temp0=str(temp0)
        

        
        
        j1 = str(j0)
        if temp0=="0":
            if (p % 12 == 1) or (p % 12 == 7):
                

                for j in range(0, p-1, 6):

                    temp1=((pervobr**j)*temp_b)%p
                    temp1 = str(temp1)
                    temp1="E:y^2=X^3+" + temp1

                    
                    rowPosition = self.izomorf_table.rowCount()
                    self.izomorf_table.insertRow(rowPosition)

                    self.izomorf_table.setItem(rowPosition, 0, QtWidgets.QTableWidgetItem(temp1))
                    self.temp1 = QtWidgets.QTableWidgetItem(temp1)
                    self.temp1.setTextAlignment(QtCore.Qt.AlignHCenter)
            if p % 12 == 5 or p % 12 == 11:
                
                for j in range(0, p-1 , 2):
                    temp1 = ((pervobr**j)*temp_b)%p
                    temp1 = str(temp1)
                    temp1 = "E:y^2=X^3+"+temp1

                    
                    rowPosition = self.izomorf_table.rowCount()
                    self.izomorf_table.insertRow(rowPosition)

                    self.izomorf_table.setItem(rowPosition, 0, QtWidgets.QTableWidgetItem(temp1))
                    self.temp1 = QtWidgets.QTableWidgetItem(temp1)
                    self.temp1.setTextAlignment(QtCore.Qt.AlignHCenter)

        elif temp0==j1:

            if (p % 12 == 1) or (p % 12 == 5):
                
                for j in range(0, p-1, 4):
                    temp2=((pervobr**j)*temp_a)%p
                    temp2 = str(temp2)
                    temp2=("E:y^2=X^3+"+temp2+"X")

                    
                    rowPosition = self.izomorf_table.rowCount()
                    self.izomorf_table.insertRow(rowPosition)

                    self.izomorf_table.setItem(rowPosition, 0, QtWidgets.QTableWidgetItem(temp2))
                    self.temp2 = QtWidgets.QTableWidgetItem(temp2)
                    self.temp2.setTextAlignment(QtCore.Qt.AlignHCenter)
            if p % 12 == 7 or p % 12 == 11:
                
                for j in range(0, p-1, 2):
                    temp2=((pervobr**j)*temp_a)%p
                    temp2 = str(temp2)
                    temp2="E:y^2=X^3+"+temp2+"X"

                    
                    rowPosition = self.izomorf_table.rowCount()
                    self.izomorf_table.insertRow(rowPosition)

                    self.izomorf_table.setItem(rowPosition, 0, QtWidgets.QTableWidgetItem(temp2))
                    self.temp2 = QtWidgets.QTableWidgetItem(temp2)
                    self.temp2.setTextAlignment(QtCore.Qt.AlignHCenter)
        else:
            mas1=izomorf_kriv(temp_a,temp_b,p)

            for i in range(0,(p-1)//2):
                a=mas1[i][0]
                a=str(a)
                b=mas1[i][1]
                b=str(b)
                temp3_3 = "E(Fp):Y^2=X^3+" + a + "X+" + b
                temp3_3=str(temp3_3)
                

                rowPosition = self.izomorf_table.rowCount()
                self.izomorf_table.insertRow(rowPosition)

                self.izomorf_table.setItem(rowPosition, 0, QtWidgets.QTableWidgetItem(temp3_3))
                self.temp3_3 = QtWidgets.QTableWidgetItem(temp3_3)
                self.temp3_3.setTextAlignment(QtCore.Qt.AlignHCenter)
        self.izomorf_table.setSelectionMode(QAbstractItemView.SingleSelection)




        return 0





    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Изоморфизм эллиптических кривых"))
        self.p_label.setText(_translate("MainWindow", "Введите модуль:"))
        self.p_btn.setText(_translate("MainWindow", "OK"))
        item = self.pred_table.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "j"))
        item = self.pred_table.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Ej,1"))
        item = self.pred_table.horizontalHeaderItem(2)
        item.setText(_translate("MainWindow", "Ej,2"))
        item = self.izomorf_table.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Изоморфные кривые↓"))

        self.pred_label.setText(_translate("MainWindow", "Представители классов↓"))
        self.izomorf_btn.setText(_translate("MainWindow", "Найти изоморфные кривые"))
        self.str_group_btn.setText(_translate("MainWindow", "Определить структуру группы"))

        self.p_btn.clicked.connect(self.on_click)
        self.str_group_btn.clicked.connect(self.on_click_strgroup)
        self.izomorf_btn.clicked.connect(self.on_click_izomorf)
        self.izomorf_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.izomorf_table.horizontalHeader().setMinimumSectionSize(0)

class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)

        layout = QtWidgets.QGridLayout(self.centralwidget)
        layout.addWidget(self.tableWidget)

        self.pred_table.setColumnCount(4)

        self.pred_table.setAlternatingRowColors(True)
        self.pred_table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)  # !!!

        self.pred_table.verticalHeader().setDefaultSectionSize(20)
        self.pred_table.horizontalHeader().setDefaultSectionSize(40)

        self.pred_table.horizontalHeader().setSectionResizeMode(
            0, QtWidgets.QHeaderView.Fixed)
        self.pred_table.horizontalHeader().setSectionResizeMode(
            1, QtWidgets.QHeaderView.Fixed)




if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.setFixedSize(800, 600)
    MainWindow.show()

    sys.exit(app.exec_())

