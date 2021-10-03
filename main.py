import sys
from types import MemberDescriptorType
from PyQt5.QtCore import flush
from PyQt5.QtWidgets import QApplication, QMainWindow
import Ui_mygui
import os
from draw import drawing
from rndalloysgenerator3D_rev2 import calculation

A1 = '' 
A2 = ''
A3 = ''
A4 = ''
A5 = ''
P1 = 20
P2 = 20
P3 = 20
P4 = 20
P5 = 20
Mode = 2  
Le = 5
Calen=1
Vac= 0
Vacp = 20
Seed=12126
textname = "MyOutput.txt"
entropy = 0

metal_list = {'Li','Na','K','Rb','Cs','Be','Mg','Ca','Sr','Ba',
              'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
              'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',    
              'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Al',
              'Al','Ga','In','Tl','Ge','Sn','Pb','Sb','Bi','Po',
              'C','Si','Se','Te','Po' }

def read():
    dataerror = 0
    global A1 
    A1 = ui.LEI1.text()
    if A1 not in metal_list:
        ui.Error_message.setText( "Inappropriate atomic label" )
        dataerror = 1
    global A2 
    A2 = ui.LEI2.text()
    if A2 not in metal_list:
        ui.Error_message.setText( "Inappropriate atomic label" )
        dataerror = 1
    global A3 
    A3 = ui.LEI3.text()
    if A3 not in metal_list:
        ui.Error_message.setText( "Inappropriate atomic label" )
        dataerror = 1
    global A4
    A4 = ui.LEI4.text()
    if A4 not in metal_list:
        ui.Error_message.setText( "Inappropriate atomic label" )
        dataerror = 1
    global A5 
    A5 = ui.LEI5.text()
    if A5 not in metal_list:
        ui.Error_message.setText( "Inappropriate atomic label" )
        dataerror = 1

    global P1 
    if str.isnumeric(ui.LEP1.text()) :
        P1 = float(ui.LEP1.text())
    else :
        ui.Error_message.setText( "Inappropriate atom proportion" )
        dataerror = 1
    
    global P2 
    if str.isnumeric(ui.LEP2.text()) :
        P2 = float(ui.LEP2.text())
    else :
        ui.Error_message.setText( "Inappropriate atom proportion" )
        dataerror = 1

    global P3 
    if str.isnumeric(ui.LEP3.text()):
        P3 = float(ui.LEP3.text())
    else :
        ui.Error_message.setText( "Inappropriate atom proportion" )
        dataerror = 1

    global P4 
    if str.isnumeric(ui.LEP4.text()):
        P4 = float(ui.LEP4.text())
    else :
        ui.Error_message.setText( "Inappropriate atom proportion" )
        dataerror = 1

    global P5 
    if str.isnumeric(ui.LEP5.text()):
        P5 = float(ui.LEP5.text())
    else :
        ui.Error_message.setText( "Inappropriate atom proportion" )
        dataerror = 1

    global Le
    if str.isnumeric(ui.Edge.text()) : 
        Le = int(ui.Edge.text())
    else :
        ui.Error_message.setText( "Inappropriate atom number on edge" )
        dataerror = 1

    global Vacp

    global Vac

    if Vac == 0:
        Vacp = 0 
    elif str.isnumeric(ui.proportionedit.text() ) and Vac == 1:
        Vacp = float(ui.proportionedit.text())
    else :
        ui.Error_message.setText( "Inappropriate random vacancy proportion" )
        dataerror = 1
    
    
    global Seed

    if str.isnumeric(ui.seeding.text() ):
        Seed = int(ui.seeding.text())
    else :
        ui.Error_message.setText( "Inappropriate Seeding number" )
        dataerror = 1
    
    

    global textname

    if  ui.textn.text()!="":
        textname = ui.textn.text() + ".txt"
    else :
        ui.Error_message.setText( "Inappropriate Textname" )
        dataerror = 1
    
    return dataerror


def Button_random():
    global Vac
    Vac = ui.random_but.isChecked() 

def Button_Entro():
    global Calen
    Calen = ui.CalE_but.isChecked() 
    print("Cal Entropy")
  
def Button_struct():
    global Mode
    if ui.radioButton_fcc.isChecked():
        print("fcc")
        Mode = 2
    elif ui.radioButton_bcc.isChecked():
        print("bcc")
        Mode = 1
    

def Click1():
    global entropy
    ui.label_10.setText("OUTPUT")
    ui.OUTPUT.setText("")
    ui.Error_message.setText( "" )
    ui.Error_message_2.setText( "" )
    dataerror = read()
    if Le % 2 == 0 : 
        ui.Error_message.setText( "Atom number on edge must be odd number !!!" )
    elif dataerror == 1:
        ui.Error_message_2.setText( "Wrong input arguments !!!" )
    else :   
        ui.Error_message.setText( "" )
        entropy = calculation(A1,A2,A3,A4,A5,P1,P2,P3,P4,P5,Mode,Le,Calen,Vac,Vacp,Seed,textname)
        entropy = round(entropy, 2)
        drawing(textname)


def Click2():

    os.system("start notepad " + textname)


def Click3():
    ui.label_10.setText("OUTPUT")
    ui.OUTPUT.setText("")
    if Mode == 1:
        struct = "bcc"
    elif Mode == 2:
        struct = "fcc"
    A_list = '{'+A1+A2+A3+A4+A5+'}' 
    proplist = str(int(P1))+":"+str(int(P2))+":"+str(int(P3))+":"+str(int(P4))+":"+str(int(P5))
    if Calen == 1:
        Calenstring = "Yes"
    elif Calen == 0:
        Calenstring = "No"
    ui.OUTPUT.setText(
                        'Textname: '+textname+'\n'
                        'Structure: ' + struct+'\n'
                        'Cal. Entropy: '+Calenstring+'\n'
                        'Atom list: '+'{'+A1+' '+A2+' '+A3+' '+A4+' '+A5+'}\n'
                        'Atom proportion:  '+proplist+'\n'
                        'Random vacancy: '+str(Vac)+'\n'
                        'Entropy of this system: '+str(entropy)+'\n'
                    )

def Click_est():
    global Le
    
    if str.isnumeric(ui.Edge.text()) : 
        Le = int(ui.Edge.text())
        if Le % 2 == 0 : 
            ui.Error_message.setText( "Atom number on edge must be odd number !!!" )
            ui.estimate_text.setText( "NaN" )
        else :   
            ui.Error_message.setText( "" )
            if (Mode==1):
                result = pow(Le,3)+pow(Le-1,3)
            elif (Mode==2):
                result = Le*(pow(Le,2)+pow(Le-1,2))+(Le-1)*Le*(Le-1)*2
            result = int(result)
            ui.estimate_text.setText(str(result))
    else :
        ui.Error_message.setText( "Inappropriate atom number on edge" )

def Click_exp():
    ui.LEI1.setText("Ni")
    ui.LEI2.setText("Co")
    ui.LEI3.setText("Mn")
    ui.LEI4.setText("Fe")
    ui.LEI5.setText("Ti")
    ui.LEP1.setText("2")
    ui.LEP2.setText("1")
    ui.LEP3.setText("3")
    ui.LEP4.setText("1")
    ui.LEP5.setText("4")
    ui.random_but.setChecked(True)
    global Vac
    Vac = True
    ui.proportionedit.setText("20")
    ui.Edge.setText("7")
    
    ui.CalE_but.setChecked(True)
    global Calen
    Calen = 1
    
    ui.radioButton_bcc.setChecked(True)
    global Mode 
    Mode = 1
    ui.textn.setText("MyOutput")
    global Seed 
    ui.seeding.setText("12126")

    ui.label_10.setText("Available atoms")
    ui.OUTPUT.setText(" Li   Na    K   Rb   Cs   Be   Mg   Ca   Sr  \n"
                      " Sc   Ti    V   Cr   Mn   Fe   Co   Ni   Cu  \n" 
                      " Zn    Y   Zr   Nb   Mo   Tc   Ru   Rh   Pd  \n" 
                      " Ag   Cd   Hf   Ta    W   Re   Os   Ir   Pt  \n"
                      " Au   Hg   Al   Al   Ga   In   Tl   Ge   Sn  Ba\n"
                      " Pb   Sb   Bi   Po    C   Si   Se   Te   Po  ")
        






    

if __name__ == '__main__':

    app = QApplication(sys.argv)
    
    MainWindow = QMainWindow()

    ui = Ui_mygui.Ui_MainWindow()
    ui.setupUi(MainWindow)


    MainWindow.show()
    MainWindow.setWindowTitle("Random Alloy Generator")

    ui.radioButton_fcc.setChecked(True)
    ui.Button1.clicked.connect(Click1)
    ui.Button2.clicked.connect(Click2)
    ui.Button3.clicked.connect(Click3)
    ui.Button_estimate.clicked.connect(Click_est)
    ui.CalE_but.clicked.connect(Button_Entro)
    ui.Button_example.clicked.connect(Click_exp)
    ui.radioButton_bcc.clicked.connect(Button_struct)
    ui.radioButton_fcc.clicked.connect(Button_struct)
    ui.random_but.clicked.connect(Button_random)



    sys.exit(app.exec_())


   
