# Library
import random
import numpy as np
import math

# mode=1:bcc / mode=2: fcc


def calculation(atom1,atom2,atom3,atom4,atom5,cpst1,cpst2,cpst3,cpst4,cpst5,mode,L,calentropy,vacancy,vacanpropor,seedingnumber,filename):
   
    entropy = 0
    vacanporporreduce=vacanpropor/100
    if mode == 1:
    # parameters setting
    # only odd number is valid in this code
        randq = 0       # 0 to 1;
        eventnumber=5   # numbers of atom

        cpstt=[cpst1, cpst2, cpst3, cpst4, cpst5]

        token1=1
        token2=2
        token3=3
        token4=4
        token5=5

        #====calculation====
        random.seed(seedingnumber)
        # re-distribution composition
        rcpst = []
        for i in range(eventnumber):
            rcpst.append(cpstt[i]/sum(cpstt))
        an_real=L**2*(L+1)/2+(L-1)**2*(L-1)/2

        # 3-D space matrix
        amap_odd=np.zeros((L,L,L))      # odd layer
        amap_even=np.zeros((L-1,L-1,L-1))  # even layer
        rmap_odd=np.random.random((L,L,L))
        rmap_even=np.random.random((L-1,L-1,L-1))

        x=0
        y=0
        z=0
        # i,j,k distribution mode 1, naturally random
        if mode == 1:
            atom1p=[]
            atom2p=[]
            atom3p=[]
            atom4p=[]
            atom5p=[]

            # odd layer loop
            for i in range(L):
                for j in range(L):
                    for k in range(L):
                        if rmap_odd[i,j,k] >= 0 and rmap_odd[i,j,k] < sum(rcpst[0:1]):  # re-distribute 1
                            amap_odd[i,j,k]=token1
                            # deciding first atom position
                            x=(i)/(L-1)
                            y=(j)/(L-1)
                            z=(k)/(L-1)
                            atom1p.append([x, y, z])

                        elif rmap_odd[i,j,k] >= sum(rcpst[0:1]) and rmap_odd[i,j,k] < sum(rcpst[0:2]):  # re-dist2
                            amap_odd[i,j,k]=token2
                            x=(i)/(L-1)
                            y=(j)/(L-1)
                            z=(k)/(L-1)
                            atom2p.append([x, y, z])

                        elif rmap_odd[i,j,k] >= sum(rcpst[0:2]) and rmap_odd[i,j,k] < sum(rcpst[0:3]):  # re-dist3
                            amap_odd[i,j,k]=token3
                            x=(i)/(L-1)
                            y=(j)/(L-1)
                            z=(k)/(L-1)
                            atom3p.append([x, y, z])

                        elif rmap_odd[i,j,k] >= sum(rcpst[0:3]) and rmap_odd[i,j,k] < sum(rcpst[0:4]):  # re-dist4
                            amap_odd[i,j,k]=token4
                            x=(i)/(L-1)
                            y=(j)/(L-1)
                            z=(k)/(L-1)
                            atom4p.append([x, y, z])

                        elif rmap_odd[i,j,k] >= sum(rcpst[0:4]) and rmap_odd[i,j,k] <= sum(rcpst[0:5]):  # re-dist5
                            amap_odd[i,j,k]=token5
                            x=(i)/(L-1)
                            y=(j)/(L-1)
                            z=(k)/(L-1)
                            atom5p.append([x, y, z])
                        else:
                            print(i)
                            print(j)
                            print(k)

            for i in range(L-1):
                for j in range(L-1):
                    for k in range((L-1)):
                        if rmap_even[i,j,k] >= 0 and rmap_even[i,j,k] < sum(rcpst[0:1]): #re-distribute 1
                            amap_even[i,j,k]=token1
                            x=(i+1/2)/(L-1)
                            y=(j+1/2)/(L-1)
                            z=(k+1/2)/(L-1)
                            atom1p.append([x, y, z])               
                            #%deciding first atom position
                        elif rmap_even[i,j,k] >= sum(rcpst[0:1]) and rmap_even[i,j,k] <sum(rcpst[0:2]):  #re-dist2
                            amap_even[i,j,k]=token2
                            x=(i+1/2)/(L-1)
                            y=(j+1/2)/(L-1)
                            z=(k+1/2)/(L-1)
                            atom2p.append([x, y, z])
                        elif rmap_even[i,j,k] >= sum(rcpst[0:2]) and rmap_even[i,j,k] <sum(rcpst[0:3]): #re-dist3
                            amap_even[i,j,k]=token3
                            x=(i+1/2)/(L-1)
                            y=(j+1/2)/(L-1)
                            z=(k+1/2)/(L-1)
                            atom3p.append([x, y, z])                
                        elif rmap_even[i,j,k] >= sum(rcpst[0:3]) and rmap_even[i,j,k] <sum(rcpst[0:4]): #re-dist4
                            amap_even[i,j,k]=token4
                            x=(i+1/2)/(L-1)
                            y=(j+1/2)/(L-1)
                            z=(k+1/2)/(L-1)
                            atom4p.append([x, y, z])          
                        elif rmap_even[i,j,k] >= sum(rcpst[0:4]) and rmap_even[i,j,k] <= sum(rcpst[0:5]): #re-dist5
                            amap_even[i,j,k]=token5
                            x=(i+1/2)/(L-1)
                            y=(j+1/2)/(L-1)
                            z=(k+1/2)/(L-1)
                            atom5p.append([x, y, z])
                        else:
                            print(i)
                            print(j)
                            print(k)


        # position saving
        # k=1 -> z=0's plane and growing toward top

        # i,j,k distribution conditionally random which based on randomness

        B=[]
        if calentropy == 1:
            count=0
            for i in range(eventnumber):
                for j in range(eventnumber):
                    for k in range(eventnumber):
                        for l in range(eventnumber):      
                            B.append([[i, j], [k, l]])    # define state
                            count=count+1;   

            # calculate possibility
            P=np.zeros(len(B))

            amap_odd_e=amap_odd
            amap_even_e=amap_even
            amap_odd_e=amap_odd
            amap_even_e=amap_even

            amap_odd_e=np.concatenate((amap_odd_e, [amap_odd_e[0]]), axis=0)

            a=np.split(amap_odd_e,[1],axis=1)
            amap_odd_e=np.concatenate((amap_odd_e, a[0]), axis=1)

            amap_even_e=np.concatenate((amap_even_e, [amap_even_e[0]]), axis=0)

            b=np.split(amap_even_e,[1],axis=1)
            amap_even_e=np.concatenate((amap_even_e, b[0]), axis=1)

            # odd
            for i in range(L-1):
                for j in range(L-1):
                    for k in range(L):
                        for kkk in range(len(B)):
                            posi=[[amap_odd_e[i,j,k], amap_odd_e[i+1,j,k]], [amap_odd_e[i,j+1,k], amap_odd_e[i+1,j+1,k]]]
                            #check=posi==B[kkk]
                            if posi==B[kkk]:
                                P[kkk]=P[kkk]+1
                            else:
                                P[kkk]=P[kkk]+0

            for i in range(L-2):
                for j in range(L-2):
                    for k in range(L-1):
                        for kkk in range(len(B)):
                            posi=[[amap_even_e[i,j,k], amap_even_e[i+1,j,k]], [amap_even_e[i,j+1,k], amap_even_e[i+1,j+1,k]]]
                            if posi==B[kkk]:
                                P[kkk]=P[kkk]+1
                            else:
                                P[kkk]=P[kkk]+0


            P=P/sum(P)  #/length(B)
            entropy=0
            for i in range(len(B)):
                if P[i]!=0:
                    entropy=entropy+ -P[i]*math.log(P[i], 2)
            print(entropy)

        an_eval=len(atom1p)+len(atom2p)+len(atom3p)+len(atom4p)+len(atom5p)
        print(an_eval)

        # ###output###
        # fid = open(filename, 'w')   # path
        # for i in range(len(atom1p)):
        #     fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom1, atom1p[i][0], atom1p[i][1], atom1p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))

        # for i in range(len(atom2p)):
        #     fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom2, atom2p[i][0], atom2p[i][1], atom2p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))

        # for i in range(len(atom3p)):
        #     fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom3, atom3p[i][0], atom3p[i][1], atom3p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))

        # for i in range(len(atom4p)):
        #     fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom4, atom4p[i][0], atom4p[i][1], atom4p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))

        # for i in range(len(atom5p)):
        #     fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom5, atom5p[i][0], atom5p[i][1], atom5p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))

        # fid.close()

    elif mode == 2:
    # for fcc
        randq = 0       # 0 to 1;
        seedingnumber=126
        eventnumber=5   # numbers of atom
        cpstt=[cpst1, cpst2, cpst3, cpst4, cpst5]
        token1=1
        token2=2
        token3=3
        token4=4
        token5=5
        mode=1

        #====calculation====
        random.seed(seedingnumber)
        # re-distribution composition
        rcpst = []
        for i in range(eventnumber):
            rcpst.append(cpstt[i]/sum(cpstt))
        an_real=L**2*(L+1)/2+(L-1)**2*(L-1)/2

        # 3-D space matrix
        amap_odd=np.zeros((L,L,L))          # odd layer 1
        amap_even=np.zeros((L-1,L-1,L))     # odd layer 2
        amap_2=np.zeros((2*L-1,2*L-1,L-1))  #fcc even
        rmap_odd=np.random.random((L,L,L))
        rmap_even=np.random.random((L-1,L-1,L))
        rmap_2=np.random.random((2*L-1,2*L-1,L-1))

        x=0
        y=0
        z=0

        for i in range(0, 2*L-1, 2):
            for j in range(0, 2*L-1, 2):
                rmap_2[i,j,:]=-5

        for i in range(1, 2*L-2, 2):
            for j in range(1, 2*L-2, 2):
                rmap_2[i,j,:]=-5

        # i,j,k distribution mode 1, naturally random
        atom1p=[]
        atom2p=[]
        atom3p=[]
        atom4p=[]
        atom5p=[]

        # odd layer loop
        for i in range(L):
            for j in range(L):
                for k in range(L):
                    if rmap_odd[i,j,k] >= 0 and rmap_odd[i,j,k] < sum(rcpst[0:1]):  # re-distribute 1
                        amap_odd[i,j,k]=token1
                        # deciding first atom position
                        x=(i)/(L-1)
                        y=(j)/(L-1)
                        z=(k)/(L-1)
                        atom1p.append([x, y, z])

                    elif rmap_odd[i,j,k] >= sum(rcpst[0:1]) and rmap_odd[i,j,k] < sum(rcpst[0:2]):  # re-dist2
                        amap_odd[i,j,k]=token2
                        x=(i)/(L-1)
                        y=(j)/(L-1)
                        z=(k)/(L-1)
                        atom2p.append([x, y, z])

                    elif rmap_odd[i,j,k] >= sum(rcpst[0:2]) and rmap_odd[i,j,k] < sum(rcpst[0:3]):  # re-dist3
                        amap_odd[i,j,k]=token3
                        x=(i)/(L-1)
                        y=(j)/(L-1)
                        z=(k)/(L-1)
                        atom3p.append([x, y, z])

                    elif rmap_odd[i,j,k] >= sum(rcpst[0:3]) and rmap_odd[i,j,k] < sum(rcpst[0:4]):  # re-dist4
                        amap_odd[i,j,k]=token4
                        x=(i)/(L-1)
                        y=(j)/(L-1)
                        z=(k)/(L-1)
                        atom4p.append([x, y, z])

                    elif rmap_odd[i,j,k] >= sum(rcpst[0:4]) and rmap_odd[i,j,k] <= sum(rcpst[0:5]):  # re-dist5
                        amap_odd[i,j,k]=token5
                        x=(i)/(L-1)
                        y=(j)/(L-1)
                        z=(k)/(L-1)
                        atom5p.append([x, y, z])
                    else:
                        print(i)
                        print(j)
                        print(k)

        # even layer loop
        for i in range(L-1):
            for j in range(L-1):
                for k in range(L):
                    if rmap_even[i,j,k] >= 0 and rmap_even[i,j,k] < sum(rcpst[0:1]): #re-distribute 1
                        amap_even[i,j,k]=token1
                        x=(i+1/2)/(L-1)
                        y=(j+1/2)/(L-1)
                        z=(k)/(L-1)
                        atom1p.append([x, y, z])               
                        #%deciding first atom position
                    elif rmap_even[i,j,k] >= sum(rcpst[0:1]) and rmap_even[i,j,k] <sum(rcpst[0:2]):  #re-dist2
                        amap_even[i,j,k]=token2
                        x=(i+1/2)/(L-1)
                        y=(j+1/2)/(L-1)
                        z=(k)/(L-1)
                        atom2p.append([x, y, z])
                    elif rmap_even[i,j,k] >= sum(rcpst[0:2]) and rmap_even[i,j,k] <sum(rcpst[0:3]): #re-dist3
                        amap_even[i,j,k]=token3
                        x=(i+1/2)/(L-1)
                        y=(j+1/2)/(L-1)
                        z=(k)/(L-1)
                        atom3p.append([x, y, z])                
                    elif rmap_even[i,j,k] >= sum(rcpst[0:3]) and rmap_even[i,j,k] <sum(rcpst[0:4]): #re-dist4
                        amap_even[i,j,k]=token4
                        x=(i+1/2)/(L-1)
                        y=(j+1/2)/(L-1)
                        z=(k)/(L-1)
                        atom4p.append([x, y, z])          
                    elif rmap_even[i,j,k] >= sum(rcpst[0:4]) and rmap_even[i,j,k] <= sum(rcpst[0:5]): #re-dist5
                        amap_even[i,j,k]=token5
                        x=(i+1/2)/(L-1)
                        y=(j+1/2)/(L-1)
                        z=(k)/(L-1)
                        atom5p.append([x, y, z])
                    else:
                        print(i)
                        print(j)
                        print(k)

        # fcc
        for i in range(2*L-1):
            for j in range(2*L-1):
                for k in range(L-1):
                    if rmap_2[i,j,k] >= 0 and rmap_2[i,j,k] < sum(rcpst[0:1]): #re-distribute 1
                        amap_2[i,j,k]=token1
                        x=(i)/(2*L-2)
                        y=(j)/(2*L-2)
                        z=(k+1/2)/(L-1)
                        atom1p.append([x, y, z])               
                        #%deciding first atom position
                    elif rmap_2[i,j,k] >= sum(rcpst[0:1]) and rmap_2[i,j,k] <sum(rcpst[0:2]):  #re-dist2
                        amap_2[i,j,k]=token2
                        x=(i)/(2*L-2)
                        y=(j)/(2*L-2)
                        z=(k+1/2)/(L-1)
                        atom2p.append([x, y, z])
                    elif rmap_2[i,j,k] >= sum(rcpst[0:2]) and rmap_2[i,j,k] <sum(rcpst[0:3]): #re-dist3
                        amap_2[i,j,k]=token3
                        x=(i)/(2*L-2)
                        y=(j)/(2*L-2)
                        z=(k+1/2)/(L-1)
                        atom3p.append([x, y, z])                
                    elif rmap_2[i,j,k] >= sum(rcpst[0:3]) and rmap_2[i,j,k] <sum(rcpst[0:4]): #re-dist4
                        amap_2[i,j,k]=token4
                        x=(i)/(2*L-2)
                        y=(j)/(2*L-2)
                        z=(k+1/2)/(L-1)
                        atom4p.append([x, y, z])          
                    elif rmap_2[i,j,k] >= sum(rcpst[0:4]) and rmap_2[i,j,k] <= sum(rcpst[0:5]): #re-dist5
                        amap_2[i,j,k]=token5
                        x=(i)/(2*L-2)
                        y=(j)/(2*L-2)
                        z=(k+1/2)/(L-1)
                        atom5p.append([x, y, z])


        # position saving
        # k=1 -> z=0's plane and growing toward top

        # i,j,k distribution conditionally random which based on randomness

        B=[]
        if calentropy == 1:
            count=0
            for i in range(eventnumber):
                for j in range(eventnumber):
                    for k in range(eventnumber):
                        for l in range(eventnumber):      
                            B.append([[i, j], [k, l]])    # define state
                            count=count+1;   

            # calculate possibility
            P=np.zeros(len(B))

            amap_odd_e=amap_odd
            amap_even_e=amap_even
            amap_2_e=amap_2

            amap_odd_e=np.concatenate((amap_odd_e, [amap_odd_e[0]]), axis=0)

            a=np.split(amap_odd_e,[1],axis=1)
            amap_odd_e=np.concatenate((amap_odd_e, a[0]), axis=1)

            amap_even_e=np.concatenate((amap_even_e, [amap_even_e[0]]), axis=0)

            b=np.split(amap_even_e,[1],axis=1)
            amap_even_e=np.concatenate((amap_even_e, b[0]), axis=1)

            amap_2_e=np.concatenate((amap_2_e, [amap_2_e[0]]), axis=0)

            b=np.split(amap_2_e,[1],axis=1)
            amap_2_e=np.concatenate((amap_2_e, b[0]), axis=1)  

            # odd
            for i in range(L-1):
                for j in range(L-1):
                    for k in range(L):
                        for kkk in range(len(B)):
                            posi=[[amap_odd_e[i,j,k], amap_odd_e[i+1,j,k]], [amap_odd_e[i,j+1,k], amap_odd_e[i+1,j+1,k]]]
                            #check=posi==B[kkk]
                            if posi==B[kkk]:
                                P[kkk]=P[kkk]+1
                            else:
                                P[kkk]=P[kkk]+0

            for i in range(L-2):
                for j in range(L-2):
                    for k in range(L):
                        for kkk in range(len(B)):
                            posi=[[amap_even_e[i,j,k], amap_even_e[i+1,j,k]], [amap_even_e[i,j+1,k], amap_even_e[i+1,j+1,k]]]
                            if posi==B[kkk]:
                                P[kkk]=P[kkk]+1
                            else:
                                P[kkk]=P[kkk]+0
            # fcc layer 2
            for i in range(2*L-2):
                for j in range(2*L-2):
                    for k in range(L-1):
                        if amap_2[i,j,k]!=0:
                            for kkk in range(len(B)):
                                posi=[[amap_2_e[i,j,k], amap_2_e[i+1,j,k]], [amap_2_e[i,j+1,k], amap_2_e[i+1,j+1,k]]]
                                if posi==B[kkk]:
                                    P[kkk]=P[kkk]+1
                                else:
                                    P[kkk]=P[kkk]+0

            P=P/sum(P)  #/length(B)
            entropy=0
            for i in range(len(B)):
                if P[i]!=0:
                    entropy=entropy+ -P[i]*math.log(P[i], 2)
            print(entropy)


        an_eval=len(atom1p)+len(atom2p)+len(atom3p)+len(atom4p)+len(atom5p)
        print(an_eval)

    ###output###
    if vacancy==1:
        fid = open(filename, 'w')   # path
        
        for i in range(len(atom1p)):
            rndnum=random.random()
            if rndnum > vacanporporreduce:
                fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom1, atom1p[i][0], atom1p[i][1], atom1p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))
        for i in range(len(atom2p)):
            rndnum=random.random()
            if rndnum > vacanporporreduce:
                fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom2, atom2p[i][0], atom2p[i][1], atom2p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))
        for i in range(len(atom3p)):
            rndnum=random.random()
            if rndnum > vacanporporreduce:
                fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom3, atom3p[i][0], atom3p[i][1], atom3p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))
        for i in range(len(atom4p)):
            rndnum=random.random()
            if rndnum > vacanporporreduce:
                fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom4, atom4p[i][0], atom4p[i][1], atom4p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))
        for i in range(len(atom5p)):
            rndnum=random.random()
            if rndnum > vacanporporreduce:
                fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom5, atom5p[i][0], atom5p[i][1], atom5p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))

        fid.close()

    else:
        fid = open(filename, 'w')   # path
        for i in range(len(atom1p)):
            fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom1, atom1p[i][0], atom1p[i][1], atom1p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))

        for i in range(len(atom2p)):
            fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom2, atom2p[i][0], atom2p[i][1], atom2p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))

        for i in range(len(atom3p)):
            fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom3, atom3p[i][0], atom3p[i][1], atom3p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))

        for i in range(len(atom4p)):
            fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom4, atom4p[i][0], atom4p[i][1], atom4p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))

        for i in range(len(atom5p)):
            fid.write('{:2s} {:8.6f} {:8.6f} {:8.6f}\n'.format(atom5, atom5p[i][0], atom5p[i][1], atom5p[i][2]))   #按列輸出，若要按行輸出：fprintf(fid,'%.4\t',A(jj))

        fid.close()
    return entropy    