from ase.visualize import view
from ase import Atoms


def drawing(file_path):
    f = open(file_path, 'r') # 開啟並讀取檔案
    lines = f.readlines() # 讀取檔案內容的每一行文字為陣列 是一堆strings的list
    atom_num = len(lines) # 原子數
    print("Atom_numbers: ",atom_num)
    
    for i in range(atom_num):  #去掉結尾換行
        lines[i] = lines[i].replace('\n', '').replace('\r', '')  
    f.close() # 關閉檔案

    Cell=[30, 30, 30]
    name_list = []
    position_list = []

    for i in range(atom_num):
        lines[i] = lines[i].split(' ')  #現在每個 lines[i] 是一個 4 個字串組成的 list
        #print(lines[i], end = '')        
        x = float(lines[i][1])*30
        y = float(lines[i][2])*30
        z = float(lines[i][3])*30    
        name_list.append(lines[i][0])
        position_list.append((x,y,z))
    All_atoms = Atoms(name_list,  positions=position_list,  cell=Cell) 
    view(All_atoms)
