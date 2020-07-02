#!/bin/env python
import os
import matplotlib.pyplot as plt
import numpy as np
import sys
from copy import deepcopy
from sklearn.linear_model import LinearRegression

def max_list(lt):
    temp = 0
    for i in lt:
        if lt.count(i) > temp:
            max_str = i
            temp = lt.count(i)
    return max_str

def dcd_rotation_dna_protein(read_path, out_path):
    print('dcd_rotation_dna_protien')
    cwd = os.getcwd()
    os.chdir(read_path)
    all_files = os.listdir('.')
    all_dcd = [f for f in all_files if 'dcd' in f]
    for dcd_file in all_dcd:
        print(dcd_file, '  reading...')
        top_file = cwd + "\\" + "top.top"
        inp_file = cwd + '\\' + "rotation_group.inp"
        out_file = out_path + '\\'+ 'rotation_' + dcd_file.replace('dcd','txt')
        cmd     = "p_dcd_rotation_dna_protein  " + \
                    ' -f ' + dcd_file +\
                    ' -s ' + top_file +\
                    ' -i ' + inp_file +\
                    ' -o ' + out_file
        temp = os.popen(cmd).readlines()

def dcd_search_mode_dna_protein(read_path,out_path):
    print('dcd_search_mode_dna_protein')
    cwd = os.getcwd()
    os.chdir(read_path)
    all_files = os.listdir('.')
    all_dcd = [f for f in all_files if 'dcd' in f]
    for dcd_file in all_dcd:
        print(dcd_file, '  reading...')
        top_file = cwd + '\\'+ "top.top"
        inp_file = cwd + '\\' + "mode_group.inp"
        out_file = out_path + '\\'+ 'mode_' + dcd_file.replace('dcd','txt')
        cmd     = "p_dcd_search_mode_dna_protein  " + \
                    ' -f ' + dcd_file +\
                    ' -s ' + top_file +\
                    ' -i ' + inp_file +\
                    ' -o ' + out_file
        temp = os.popen(cmd).readlines()

def dcd_interface_dna_protein_full(read_path,out_path):
    print('\n')
    print('dcd_interface_dna_protein_full')
    cwd = os.getcwd()
    os.chdir(read_path)
    all_files = os.listdir('.')
    all_dcd = [f for f in all_files if 'dcd' in f]
    for dcd_file in all_dcd:
        print(dcd_file, '  reading...')
        top_file = cwd + '\\'+ "top.top"
        inp_file = cwd + '\\' + "interface_group_full.inp"
        out_file = out_path + '\\'+ 'fulll_interface_' + dcd_file.replace('dcd','txt')
        cmd     = "p_dcd_interface_lig_rec  " + \
                    ' -f ' + dcd_file +\
                    ' -s ' + top_file +\
                    ' -i ' + inp_file +\
                    ' -o ' + out_file
        temp = os.popen(cmd).readlines()

def dcd_interface_dna_protein_truncated(read_path,out_path):
    print('\n')
    print('dcd_interface_dna_protein_truncated')
    cwd = os.getcwd()
    os.chdir(read_path)
    all_files = os.listdir('.')
    all_dcd = [f for f in all_files if 'dcd' in f]
    for dcd_file in all_dcd:
        print(dcd_file, '  reading...')
        top_file = cwd + '\\'+ "top.top"
        inp_file = cwd + '\\' + "interface_group_truncated.inp"
        out_file = out_path + '\\'+ 'truncated_interface_' + dcd_file.replace('dcd','txt')
        cmd     = "p_dcd_interface_lig_rec  " + \
                    ' -f ' + dcd_file +\
                    ' -s ' + top_file +\
                    ' -i ' + inp_file +\
                    ' -o ' + out_file
        temp = os.popen(cmd).readlines()

def dcd_distances_com(read_path,out_path):
    print('\n')
    print('dcd_distances_com')
    cwd = os.getcwd()
    os.chdir(read_path)
    all_files = os.listdir('.')
    all_dcd = [f for f in all_files if 'dcd' in f]
    for dcd_file in all_dcd:
        print(dcd_file, '  reading...')
        top_file = cwd + '\\'+ "top.top"
        inp_file = cwd + '\\' + "com_dis_group.inp"
        out_file = out_path + '\\'+ 'com_dis_' + dcd_file.replace('dcd','txt')
        cmd     = "p_dcd_distances_com  " + \
                    ' -f ' + dcd_file +\
                    ' -s ' + top_file +\
                    ' -i ' + inp_file +\
                    ' -o ' + out_file
        temp = os.popen(cmd).readlines()

def dcd_distance_rec_lig_n(read_path,out_path):
    print('\n')
    print('dcd_distance_lig_rec_N')
    cwd = os.getcwd()
    os.chdir(read_path)
    all_files = os.listdir('.')
    all_dcd = [f for f in all_files if 'dcd' in f]
    for dcd_file in all_dcd:
        print(dcd_file, '  reading...')
        top_file = cwd + '\\'+ "top.top"
        inp_file = cwd + '\\' + "dis_rec_lig_N.inp"
        out_file = out_path + '\\'+ 'rl_n_dis_' + dcd_file.replace('dcd','txt')
        cmd     = "p_dcd_distance_lig_rec  " + \
                    ' -f ' + dcd_file +\
                    ' -s ' + top_file +\
                    ' -i ' + inp_file +\
                    ' -o ' + out_file
        temp = os.popen(cmd).readlines()

def dcd_distance_rec_lig_dbd(read_path,out_path):
    print('\n')
    print('dcd_distance_lig_rec_dbd')
    cwd = os.getcwd()
    os.chdir(read_path)
    all_files = os.listdir('.')
    all_dcd = [f for f in all_files if 'dcd' in f]
    for dcd_file in all_dcd:
        print(dcd_file, '  reading...')
        top_file = cwd + '\\'+ "top.top"
        inp_file = cwd + '\\' + "dis_rec_lig_DBD.inp"
        out_file = out_path + '\\'+ 'rl_dbd_dis_' + dcd_file.replace('dcd','txt')
        cmd     = "p_dcd_distance_lig_rec  " + \
                    ' -f ' + dcd_file +\
                    ' -s ' + top_file +\
                    ' -i ' + inp_file +\
                    ' -o ' + out_file
        temp = os.popen(cmd).readlines()

def dcd_distance_rec_lig_c(read_path,out_path):
    print('\n')
    print('dcd_distance_lig_rec_C')
    cwd = os.getcwd()
    os.chdir(read_path)
    all_files = os.listdir('.')
    all_dcd = [f for f in all_files if 'dcd' in f]
    for dcd_file in all_dcd:
        print(dcd_file, '  reading...')
        top_file = cwd + '\\'+ "top.top"
        inp_file = cwd + '\\' + "dis_rec_lig_C.inp"
        out_file = out_path + '\\'+ 'rl_c_dis_' + dcd_file.replace('dcd','txt')
        cmd     = "p_dcd_distance_lig_rec  " + \
                    ' -f ' + dcd_file +\
                    ' -s ' + top_file +\
                    ' -i ' + inp_file +\
                    ' -o ' + out_file
        temp = os.popen(cmd).readlines()

def dcd_contact_number_n(read_path,out_path):
    print('\n')
    print('dcd_contact_number_n')
    cwd = os.getcwd()
    os.chdir(read_path)
    all_files = os.listdir('.')
    all_dcd = [f for f in all_files if 'dcd' in f]
    for dcd_file in all_dcd:
        print(dcd_file, '  reading...')
        top_file = cwd + '\\'+ "top.top"
        inp_file = cwd + '\\' + "contact_group_N.inp"
        out_file = out_path + '\\'+ 'contact_n_' + dcd_file.replace('dcd','txt')
        cmd     = "p_dcd_contact_number  " + \
                    ' -f ' + dcd_file +\
                    ' -s ' + top_file +\
                    ' -i ' + inp_file +\
                    ' -o ' + out_file
        temp = os.popen(cmd).readlines()

def dcd_contact_number_dbd(read_path,out_path):
    print('\n')
    print('dcd_contact_number_dbd')
    cwd = os.getcwd()
    os.chdir(read_path)
    all_files = os.listdir('.')
    all_dcd = [f for f in all_files if 'dcd' in f]
    for dcd_file in all_dcd:
        print(dcd_file, '  reading...')
        top_file = cwd + '\\'+ "top.top"
        inp_file = cwd + '\\' + "contact_group_DBD.inp"
        out_file = out_path + '\\'+ 'contact_dbd_' + dcd_file.replace('dcd','txt')
        cmd     = "p_dcd_contact_number  " + \
                    ' -f ' + dcd_file +\
                    ' -s ' + top_file +\
                    ' -i ' + inp_file +\
                    ' -o ' + out_file
        temp = os.popen(cmd).readlines()

def dcd_contact_number_c(read_path,out_path):
    print('\n')
    print('dcd_contact_number_c')
    cwd = os.getcwd()
    os.chdir(read_path)
    all_files = os.listdir('.')
    all_dcd = [f for f in all_files if 'dcd' in f]
    for dcd_file in all_dcd:
        print(dcd_file, '  reading...')
        top_file = cwd + '\\'+ "top.top"
        inp_file = cwd + '\\' + "contact_group_C.inp"
        out_file = out_path + '\\'+ 'contact_c_' + dcd_file.replace('dcd','txt')
        cmd     = "p_dcd_contact_number  " + \
                    ' -f ' + dcd_file +\
                    ' -s ' + top_file +\
                    ' -i ' + inp_file +\
                    ' -o ' + out_file
        temp = os.popen(cmd).readlines()

def dcd_dna_curvature(read_path,out_path):
    print('\n')
    print('dcd_dna_curvature')
    cwd = os.getcwd()
    os.chdir(read_path)
    all_files = os.listdir('.')
    all_dcd = [f for f in all_files if 'dcd' in f]
    for dcd_file in all_dcd:
        print(dcd_file, '  reading...')
        top_file = cwd + '\\'+ "top.top"
        inp_file = cwd + '\\' + "dna_curvature.inp"
        out_file = out_path + '\\'+ 'curvature_' + dcd_file.replace('dcd','txt')
        cmd     = "p_dcd_dna_curvature  " + \
                    ' -f ' + dcd_file +\
                    ' -s ' + top_file +\
                    ' -i ' + inp_file +\
                    ' -o ' + out_file
        temp = os.popen(cmd).readlines()

def fun_001():
    fr = open('top.top','r')
    all_lines = fr.readlines()
    fr.close()
    grp1 = []
    grp2 = []
    i = 0
    for each in all_lines:
        if 'DS' in each:
            i = i + 1
            if i < 101:
                grp1.append(each.split()[0])
            else:
                grp2.append(each.split()[0])
         
    fw = open('rotation_group.inp','w')
    fw.write('GROUP1:')
    for g in grp1[:-1]:
        fw.write(str(g)+',')
    fw.write(str(grp1[-1]) + '\n')
   
    grp2.reverse() 
    fw.write('GROUP2:')
    for g in grp2[:-1]:
        fw.write(str(g)+',')
    fw.write(str(grp2[-1]) + '\n')

    fw.write('GROUP3:')
    for g in range(599,871):
        fw.write(str(g)+',')
    fw.write('871\n')
    fw.close()
#fun_001()
#exit(0)

def hist_residue_interface_dna(out_path,figure_path):
    print('hist_residue_interface_dna')
    discard_line =  2000
    total_length = 732
    initi_length = 97
    os.chdir(out_path)
    all_files = os.listdir('.')
    inpfile_str=["interface_xpa_dna"]
    outfile_str={}
    for s in inpfile_str:
        for f in all_files:
            if not( s in f):
                continue
            con = f.split('_')[-2]
            name = s + '_' + con
            if name in outfile_str.keys():
                outfile_str[name].append(f)
            else:
                outfile_str[name] = []
                outfile_str[name].append(f)
              
    i = 0  
    for name in outfile_str.keys():
        i = i + 1
        residues_list = [0]*(total_length-598)
        line_num = 0
        for f in outfile_str[name]:
            print(f,' ...reading...')
            fr = open(f,'r')
            all_lines = fr.readlines()
            fr.close()
            lig_lines = [ line[8:] for line in all_lines[1+discard_line*4::4] if line[8:] != " \n"]           
            line_num += len(lig_lines)           
            for line in lig_lines:
                for each in line.split():
                    residues_list[int(each)-1 -598] += 1

        residues_list = np.array(residues_list)
        p_list = residues_list*1.0/line_num
        fig = plt.figure(i)
        ax = fig.add_subplot(1,1,1)
        ax.plot(range(1+initi_length,total_length-598+initi_length+1),p_list)
        ax.set_title(name)
        ax.set_ylabel('probability')
        ax.set_xlabel('residue id')
        if i == 1:
            t_p = sorted(p_list)[-20]
            ind = np.where(p_list>t_p)
            print(np.array(ind)+initi_length+1)
        plt.savefig(figure_path + r'\hist_residue_interface_dna_'+ name.split('_')[-1] + '.png')
    plt.show()
    
def traj_distance_rec_lig(out_path,figure_path):
    print('traj_distance_rec_lig')
    discard_line =  2000
    os.chdir(out_path)
    all_files = os.listdir('.')
    inpfile_str=["rl_dis_xpa_dna"]
    outfile_str={}
    for s in inpfile_str:
        for f in all_files:
            if not( s in f):
                continue
            con = f.split('_')[-2]
            name = s + '_' + con
            if name in outfile_str.keys():
                outfile_str[name].append(f)
            else:
                outfile_str[name] = []
                outfile_str[name].append(f)
              
    i = 0  
    for name in outfile_str.keys():
        i = i + 1
        j = 0
        for f in outfile_str[name]:
            j = j + 1
            if j != 4 :
                continue
            
            
            print(f,' ...reading...')
            fr = open(f,'r')
            all_lines = fr.readlines()
            fr.close()
            dis = [ float(line.split()[1]) for line in all_lines[discard_line:]]    
            dis = np.array(dis)       

            fig = plt.figure(i)
            ax = fig.add_subplot(1,1,1)
            ax.plot(range(0,len(dis)),dis)
            ax.set_title(f)
            ax.set_ylabel('distance(A)')
            ax.set_xlabel('md_step')
            plt.savefig(figure_path + r'\traj_' + f.replace('txt','png'))
    plt.show()
    
def hist_distance_rec_lig(out_path,figure_path):
    print('hist_distance_rec_lig')
    discard_line =  2000
    
    os.chdir(out_path)
    all_files = os.listdir('.')
    inpfile_str=["rl_dis_xpa_dna"]
    outfile_str={}
    for s in inpfile_str:
        for f in all_files:
            if not( s in f):
                continue
            con = f.split('_')[-2]
            name = s + '_' + con
            if name in outfile_str.keys():
                outfile_str[name].append(f)
            else:
                outfile_str[name] = []
                outfile_str[name].append(f)
       
    _p = {}
    i = 0  
    for name in outfile_str.keys():
        i = i + 1
        dis = []
        for f in outfile_str[name]:

            print(f,' ...reading...')
            fr = open(f,'r')
            all_lines = fr.readlines()
            fr.close()
            dis.extend( [ float(line.split()[1]) for line in all_lines[discard_line:]] )
            
        dis = np.array(dis)    
        _p[name] = np.sum(dis<50)/len(dis)
        fig = plt.figure(i)
        ax = fig.add_subplot(1,1,1)
        n, bins, patches = ax.hist(dis, bins=100, normed=1,edgecolor='None',facecolor='red') 
        ax.set_title(name)
        ax.set_xlabel('distance(A)')
        ax.set_ylabel('probability')
        plt.savefig(figure_path + '\hist_distance_rec_lig_' + name.split('_')[-1] + '.png')
    print(_p)
    #plt.show()
            
# Analysis of the electrostatic surface potential revealed that three positively charged patches created by Lys141/179, Arg207/211, and Lys151 interact with the backbones of DNA duplex
#(F.-M. Lian, X. Yang, Y.-L. Jiang, et al @International Journal of Biological Macromolecules(2020))
# 141, 179,207,211,151 used to define the search mode.....
 
def hist_distance_rec_lig_98_105(out_path,figure_path):
    print("\nhist_distance_rec_lig_98_105")
    discard_line =  2000
    os.chdir(out_path)
    all_files = os.listdir('.')
    inpfile_str=["rl_n_dis_xpa_dna","rl_dbd_dis_xpa_dna","rl_c_dis_xpa_dna"]
    inpfile_str=["98_105_rl_dis_xpa_dna"]
    outfile_str={}
    for s in inpfile_str:
        for f in all_files:
            if not( s in f):
                continue
            con = f.split('_')[-2]
            name = s + '_' + con
            if name in outfile_str.keys():
                outfile_str[name].append(f)
            else:
                outfile_str[name] = []
                outfile_str[name].append(f)
             
    _p = {}
    i = 0  
    for name in outfile_str.keys():
        i = i + 1
        dis = []
        for f in outfile_str[name]:

            print(f,' ...reading...')
            fr = open(f,'r')
            all_lines = fr.readlines()
            fr.close()
            dis.extend( [ float(line.split()[1]) for line in all_lines[discard_line:]] )
            
        dis = np.array(dis)       
        _p[name] = np.sum(dis<50)/len(dis)
        fig = plt.figure(i)
        ax = fig.add_subplot(1,1,1)
        n, bins, patches = ax.hist(dis, bins=100, normed=1,edgecolor='None',facecolor='red')
        ax.set_title(name)
        ax.set_xlabel('distance(A)')
        ax.set_ylabel('probability')
        plt.savefig(figure_path + '\hist_distance_rec_lig_98_105_' + name.split('_')[-1] + '.png')
    print(_p)
##<<<<<<<<<-------------------------
##<<<<<<<<<-------------------------
## the snapshots with the contact between DNA and protein is considered. 
def diffusion_1D(out_path,figure_path):
    print('diffusion_1D')
    discard_line =  2000
    ss_dna_len = 299
    d_t = 100
    g_t = 2
    os.chdir(out_path)
    all_files = os.listdir('.')
    inpfile_str=["interface_xpa_dna"]
    outfile_str={}
    for s in inpfile_str:
        for f in all_files:
            if not( s in f):
                continue
            con = f.split('_')[-2]
            name = s + '_' + con
            if name in outfile_str.keys():
                outfile_str[name].append(f)
            else:
                outfile_str[name] = []
                outfile_str[name].append(f)
    
    con = [float(ele.split('_')[-1]) for ele in list(outfile_str.keys())]
    con = sorted(con)
    names = [ inpfile_str[0] + "_" + str('%.2f'%c) for c in con]
    fig_ii = 0
    Total_D1 = []
    for name in names:
        file_pos = []
        for f in outfile_str[name]:
            print(f,' ...reading...')
            fr = open(f,'r')
            all_lines = fr.readlines()
            fr.close()
            lig_lines = [ line[9:] for line in all_lines[2+discard_line*4::4]]            

            bp_all = []
            for line in lig_lines:
                if line == "\n":
                    bp_all.append(-99999)
                    continue
                
                bp = []
                for each in line.split():
                    if int(each) > ss_dna_len:
                        t =  int((ss_dna_len*2 - int(each))/3) + 1 
                    else:
                        t =  int(int(each)/3) + 1
                    if not (t in bp):
                        bp.append(t)
                        
                bp_all.append(np.mean(bp))
                
            bp_ind = []  
            ind_pos = {}
            start_flag = False
            for i in range(bp_all.__len__()):
                pos = bp_all[i]
                if pos != -99999:
                    if start_flag:
                        ind_pos[i+1] = pos
                    else:
                        start_flag = True
                        ind_pos[i+1] = pos
                else:
                    if not start_flag:
                        continue
                    if i + g_t + 1 < bp_all.__len__():
                        for ii in range(i+1, i+g_t+2):
                            if bp_all[ii] != -99999:
                                break
                        if ii == i + g_t + 1:
                            start_flag = False
                            if ind_pos.__len__() < d_t:
                                ind_pos = {}
                            else:
                                bp_ind.append(deepcopy(ind_pos))
                                ind_pos = {}
                                continue
                    else:
                        break
                
            if start_flag == True and ind_pos.__len__() > d_t:
                bp_ind.append(deepcopy(ind_pos))     
                
            if bp_ind.__len__() > 0:
                file_pos.append(deepcopy(bp_ind))
        
        D_t = [d_t/10,d_t/4,d_t/2,d_t, d_t+d_t/2, d_t+(d_t*2)/2, d_t+(d_t*3)/2, d_t+(d_t*4)/2, d_t+ (d_t*5)/2, d_t+ (d_t*6)/2  ] 
        all_msd = []
        D1 = []
        for d in D_t:
            SD = []
            for bp_ind in file_pos:
                for pos_key in bp_ind:
                    for k in sorted(list(pos_key.keys())):
                        if k + d in pos_key.keys():
                            SD.append(((pos_key[k+d] - pos_key[k])*3.4)**2)       #1bp ~ 3.4A
                    
            MSD = np.mean(SD)
            all_msd.append(MSD)   #10^8 step ~ 0.1ms
        
        fig = plt.figure(fig_ii)
        ax = fig.add_subplot(1,1,1)
        ax.set_title("truncated")
        #ax.set_xlabel(r'$dt 10^{-3} ms$')
        #ax.set_ylabel(r'$MSD(Å^{2})$')
        ax.set_xlabel(r'$log(t/10^{-3} ms)$')
        ax.set_ylabel(r'$log(MSD/Å^{2})$')
        
        
        D1 = np.array(all_msd)*1.0/(np.array(D_t)/d_t)*10**(6)  #10^8 step ~ 0.1ms, 10^6 step ~ 10^(-6)s
        Total_D1.append(np.mean(D1))
        
        #D_t.insert(0,0)
        #all_msd.insert(0,0)
        ax.plot(np.log(np.array(D_t)/d_t),np.log(np.array(all_msd)),label=name.split('_')[-1])
        ax.legend(loc="upper left")
    
    #plt.savefig(figure_path + '\diffusion_1D.png')
    plt.savefig(figure_path + "\diffusion_1D_log"+'.png')
    fig = plt.figure(fig_ii+1)
    ax = fig.add_subplot(1,1,1)
    ax.plot(np.array(con),np.array(Total_D1))
    ax.set_title("truncated")
    ax.set_xlabel(r'ion concentration(M)')
    ax.set_ylabel(r'$D1(Å^{2}/s)$')
    plt.savefig(figure_path + "\diffusion_1D_hist.png")
    print(con)
    print(Total_D1)
      
    #plt.show()

def hist_rotation(out_path,figure_path):
    print('hist_rotation')
    discard_line =  2000
    os.chdir(out_path)
    all_files = os.listdir('.')
    
    inpfile_str = "rotation_xpa_dna"
    str_contact = "contact_xpa_dna"
    outfile_str={}
    
    for f in all_files:
        if not( inpfile_str in f):
            continue
        con = f.split('_')[-2]
        name = inpfile_str + '_' + con
        if name in outfile_str.keys():
            outfile_str[name].append(f)
        else:
            outfile_str[name] = []
            outfile_str[name].append(f)
    
    i = 0  
    for name in outfile_str.keys():
        i = i + 1
        rotation = []
        y_projection = []
        for f in outfile_str[name]:

            print(f,' ...reading...')
            fr = open(f.replace(inpfile_str,str_contact),'r')
            all_lines = fr.readlines()
            fr.close()
            contact_num = [ int(line.split()[1]) for line in all_lines[discard_line:] ]
   
            fr = open(f,'r')
            all_lines = fr.readlines()
            fr.close()
            t_a = [ int(line.split()[3]) for line in all_lines[discard_line+1:]]
                           
            con_start = False
            start_ind = []
            end_ind = []
            for num_ii in range(contact_num.__len__()):
                if contact_num[num_ii] > 0 and not con_start :
                    con_start = True
                    start_ind.append(num_ii)
                if con_start and contact_num[num_ii] == 0:
                    con_start = False
                    end_ind.append(num_ii-1)
            
            if con_start:
                end_ind.append(num_ii)    
            

            
            t_p = [ float(line.split()[1]) for line in all_lines[discard_line+1:]]
            t_r = [ float(line.split()[2]) for line in all_lines[discard_line+1:]]
            
            for t_ind in range(end_ind.__len__()):
                axis_flag = False
                for x in range(start_ind[t_ind],end_ind[t_ind]+1):
                    if t_a[x] != 3:
                        continue
                    if not axis_flag:
                        axis_flag = True
                        x_start = x
                    y_projection.append(t_p[x] - t_p[x_start])
                    rotation.append(t_r[x] - t_r[x_start])
                    #y_projection.append(t_p[x] - t_p[start_ind[t_ind]])
                    #rotation.append(t_r[x] - t_r[start_ind[t_ind]])
                    
        
        fig = plt.figure(i)
        ax = fig.add_subplot(1,1,1)
        ax.scatter(rotation[::100],y_projection[::100])
        ax.set_title(name)
        ax.set_xlabel(r'$degree$')
        ax.set_ylabel('y_projection(A)')
        
        
        x = np.array(rotation[::100]).reshape((-1, 1))
        y = np.array(y_projection[::100])
        model = LinearRegression()
        model = model.fit(x,y)
        
        '''get result
        y = b0 + b1x
        '''
        
        r_sq = model.score(x, y)
        print('coefficient of determination(R^2):',r_sq)
        print('intercept:', model.intercept_)
        # scalar b0 intercept: this will be an array when y is also 2-dimensional
        print('slope:', model.coef_)
        # array slope b1 slope: ---------this will be 2-d array when y is also 2-dimensional
        x_a = np.arange(x.min(),x.max())
        y_a = model.coef_[0] * x_a + model.intercept_
        ax.plot(x_a,y_a,'r')
        ax.text(x.min()-30,y.max(),s='R: '+str('%.3f'%np.sqrt(r_sq)))
        ax.text(x.min()-30,y.max()-5,s='k: '+str('%.3f'%model.coef_[0]))
        plt.savefig(figure_path + '\hist_rotation_'+ name.split('_')[-1] + '.png') 
        
        
    #plt.show()
##<<<<<<<<<-------------------------
##<<<<<<<<<-------------------------
## the snapshots with the contact between DNA and protein is considered.
    
    
    

#<!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#only certain snapshots is considered, where the contact between DNA and protein is formed at certain cumstance. . 
def diffusion_1D_specific(out_path,figure_path):
    print('diffusion_1D_specific')
    d_t = 100
    g_t = 4
    discard_line =  2000
    ss_dna_len = 299
    os.chdir(out_path)
    all_files = os.listdir('.')
    inpfile_str=["interface_xpa_dna"]
    outfile_str={}
    for s in inpfile_str:
        for f in all_files:
            if not( s in f):
                continue
            con = f.split('_')[-2]
            name = s + '_' + con
            if name in outfile_str.keys():
                outfile_str[name].append(f)
            else:
                outfile_str[name] = []
                outfile_str[name].append(f)
    
    con = [float(ele.split('_')[-1]) for ele in list(outfile_str.keys())]
    con = sorted(con)
    names = [ inpfile_str[0] + "_" + str('%.2f'%c) for c in con]
    fig_ii = 0
    Total_D1 = []
    for name in names:

        file_pos = []
        for f in outfile_str[name]:
            print(f,' ...reading...')
            fr = open(f,'r')
            all_lines = fr.readlines()
            fr.close()
            lig_lines = [ line[9:] for line in all_lines[2+discard_line*4::4]]

            fr = open(f.replace(inpfile_str[0], "mode_xpa_dna"),'r')
            all_lines = fr.readlines()
            fr.close()
            modes = [ int(line.split()[1]) for line in all_lines[discard_line:]]
        
            bp_all = []
            line_ii = -1
            for line in lig_lines:
                line_ii += 1
                
                if line == "\n" or modes[line_ii] == 2:
                    bp_all.append(-99999)
                    continue
                
                bp = []
                for each in line.split():
                    if int(each) > ss_dna_len:
                        t =  int((ss_dna_len*2 - int(each))/3) + 1 
                    else:
                        t =  int(int(each)/3) + 1
                    if not (t in bp):
                        bp.append(t)
                        
                bp_all.append(np.mean(bp))
                
            bp_ind = []  
            ind_pos = {}
            start_flag = False
            for i in range(bp_all.__len__()):
                pos = bp_all[i]
                if pos != -99999:
                    if start_flag:
                        ind_pos[i+1] = pos
                    else:
                        start_flag = True
                        ind_pos[i+1] = pos
                else:
                    if not start_flag:
                        continue
                    if i + g_t + 1 < bp_all.__len__():
                        for ii in range(i+1, i+g_t+2):
                            if bp_all[ii] != -99999:
                                break
                        if ii == i + g_t + 1:
                            start_flag = False
                            if ind_pos.__len__() < d_t:
                                ind_pos = {}
                            else:
                                bp_ind.append(deepcopy(ind_pos))
                                ind_pos = {}
                                continue
                    else:
                        break
                
            if start_flag == True and ind_pos.__len__() > d_t:
                bp_ind.append(deepcopy(ind_pos))     
                
            if bp_ind.__len__() > 0:
                file_pos.append(deepcopy(bp_ind))
        
        D_t = [d_t/10,d_t/4,d_t/2, d_t, d_t+d_t/2, d_t+(d_t*2)/2, d_t+(d_t*3)/2, d_t+(d_t*4)/2, d_t+ (d_t*5)/2, d_t+ (d_t*6)/2  ]
        all_msd = []
        D1 = []
        for d in D_t:
            SD = []
            for bp_ind in file_pos:
                for pos_key in bp_ind:
                    for k in sorted(list(pos_key.keys())):
                        if k + d in pos_key.keys():
                            SD.append(((pos_key[k+d] - pos_key[k])*3.4)**2)       #1bp ~ 3.4A
                    
            MSD = np.mean(SD)
            all_msd.append(MSD)   #10^8 step ~ 0.1ms
        
        fig = plt.figure(fig_ii)
        ax = fig.add_subplot(1,1,1)
        ax.set_title("truncated")
        #ax.set_xlabel(r'$dt 10^{-3} ms$')
        #ax.set_ylabel(r'$MSD(Å^{2})$')
        ax.set_xlabel(r'$log(t/10^{-3} ms)$')
        ax.set_ylabel(r'$log(MSD/Å^{2})$')
        
        
        D1 = np.array(all_msd)*1.0/(np.array(D_t)/d_t)*10**(6)  #10^8 step ~ 0.1ms, 10^6 step ~ 10^(-6)s
        Total_D1.append(np.mean(D1))
        
        #D_t.insert(0,0)
        #all_msd.insert(0,0)
        ax.plot(np.log(np.array(D_t)/d_t),np.log(np.array(all_msd)),label=name.split('_')[-1])
        ax.legend(loc="upper left")
    
    #plt.savefig(figure_path + '\diffusion_1D_specific.png')
    plt.savefig(figure_path + '\diffusion_1D_specific_log_2.png')
 
    fig = plt.figure(fig_ii+1)
    ax = fig.add_subplot(1,1,1)
    ax.plot(np.array(con),np.array(Total_D1))
    ax.set_title("DBD")
    ax.set_xlabel(r'ion concentration(M)')
    ax.set_ylabel(r'$D1(A^{2}/s)$')
    plt.savefig(figure_path + '\diffusion_1D_specific_hist_2.png')   
    
    print(con)
    print(Total_D1)
      
    plt.show()

def hist_rotation_specific(out_path,figure_path):
    print('hist_rotation_specific')
    d_t = 10
    ddr = 10
    ddy = 2
    
    discard_line =  2000
    os.chdir(out_path)
    all_files = os.listdir('.')
    inpfile_str = "rotation_xpa_dna"
    str_contact = "contact_xpa_dna"
    str_mode    = "mode_xpa_dna"
    outfile_str={}
    
    for f in all_files:
        if not( inpfile_str in f):
            continue
        con = f.split('_')[-2]
        name = inpfile_str + '_' + con
        if name in outfile_str.keys():
            outfile_str[name].append(f)
        else:
            outfile_str[name] = []
            outfile_str[name].append(f)
    
    i = 0  
    for name in outfile_str.keys():
        if not '0.05' in name:
            continue
        
        i = i + 1
        rotation = []
        y_projection = []
        for f in outfile_str[name]:
            if not '0.05_001' in f:
                continue
            
            print(f,' ...reading...')
            fr = open(f.replace(inpfile_str,str_contact),'r')
            all_lines = fr.readlines()
            fr.close()
            contact_num = [ int(line.split()[1]) for line in all_lines[discard_line:] ]
   
            fr = open(f,'r')
            all_lines = fr.readlines()
            fr.close()
            t_a = [ int(line.split()[3]) for line in all_lines[discard_line+1:]]            
            t_p = [ float(line.split()[1]) for line in all_lines[discard_line+1:]]
            t_r = [ float(line.split()[2]) for line in all_lines[discard_line+1:]]
                   
            fr = open(f.replace(inpfile_str,str_mode),'r')
            all_lines = fr.readlines()
            fr.close()
            modes = [ int(line.split()[1]) for line in all_lines[discard_line:]]
                           
            con_start = False
            start_ind = []
            end_ind = []
            for num_ii in range(contact_num.__len__()):
                if contact_num[num_ii] > 0 and (not con_start) and modes[num_ii] != 2 :
                    con_start = True
                    start_ind.append(num_ii)
                if con_start and ( contact_num[num_ii] == 0 or modes[num_ii] == 2 ):
                    con_start = False       
                    if num_ii - start_ind[-1] < d_t:
                        start_ind.pop(-1)
                    else:
                        end_ind.append(num_ii-1)
            
            if con_start:
                if num_ii - start_ind[-1] < d_t:
                    start_ind.pop(-1)
                else:
                    end_ind.append(num_ii)    
            
            print(start_ind)
            print(end_ind)
            T_axis = []
            for t_ind in range(end_ind.__len__()):
                for x in range(start_ind[t_ind],end_ind[t_ind]+1):
                    T_axis.append(t_a[x])
            m = max_list(T_axis)
            
            for t_ind in range(end_ind.__len__()):
                
            #for t_ind in range(5,6):
                print('....')
                axis_flag = False
                last_a = -100
                for x in range(start_ind[t_ind],end_ind[t_ind]+1):
                    if t_a[x] != m:
                        last_a = t_a[x]
                        continue
                    #    print('end: ',x)
                    #    break
                    
                    if t_a[x] != last_a:
                        axis_flag = False
                        #continue
                    
                    if not axis_flag:
                        axis_flag = True
                        x_start = x
                        last_a = t_a[x]
                        continue
                        #print('start: ',x_start)
                            
                    y_projection.append(t_p[x] - t_p[x_start])
                    rotation.append(t_r[x] - t_r[x_start])
                    
                    last_a = t_a[x]
                    #y_projection.append(t_p[x] - t_p[start_ind[t_ind]])
                    #rotation.append(t_r[x] - t_r[start_ind[t_ind]])
        
            
                fig = plt.figure(t_ind)
                ax = fig.add_subplot(1,1,1)
                ax.scatter(rotation[::1],y_projection[::1])
                ax.set_title(name)
                ax.set_xlabel(r'$degree$')
                ax.set_ylabel('y_projection(A)')
        
        
        """
        x = np.array(rotation[::1]).reshape((-1, 1))
        y = np.array(y_projection[::1])
        model = LinearRegression()
        model = model.fit(x,y)
        
        '''get result
        y = b0 + b1x
        '''
        
        r_sq = model.score(x, y)
        print('coefficient of determination(R^2):',r_sq)
        print('intercept:', model.intercept_)
        # （标量） 系数b0 intercept: this will be an array when y is also 2-dimensional
        print('slope:', model.coef_)
        # （数组）斜率b1 slope: ---------this will be 2-d array when y is also 2-dimensional
        x_a = np.arange(x.min(),x.max())
        y_a = model.coef_[0] * x_a + model.intercept_
        ax.plot(x_a,y_a,'r')
        #ax.text(x.min()-30,y.max(),s='R: '+str('%.3f'%np.sqrt(r_sq)))
        #ax.text(x.min()-30,y.max()-5,s='k: '+str('%.3f'%model.coef_[0]))
        """
        #plt.savefig(figure_path + '\hist_rotation_specific_'+ name.split('_')[-1] + '.png') 
    plt.show()

def hist_mode_specific(out_path,figure_path):
    print('hist_mode_specific')
    discard_line =  2000
    os.chdir(out_path)
    all_files = os.listdir('.')
    
    inpfile_str = "mode_xpa_dna"
    outfile_str={}
    
    for f in all_files:
        if not( inpfile_str in f):
            continue
        con = f.split('_')[-2]
        name = inpfile_str + '_' + con
        if name in outfile_str.keys():
            outfile_str[name].append(f)
        else:
            outfile_str[name] = []
            outfile_str[name].append(f)
    
    sliding = {}
    hopping = {}
    D3      = {}
    
    fig_s = plt.figure(1)
    fig_h = plt.figure(2)
    fig_d = plt.figure(3)
    ax_s = fig_s.add_subplot(1,1,1)
    ax_h = fig_h.add_subplot(1,1,1)
    ax_d = fig_d.add_subplot(1,1,1)
    i = 0  
    for name in outfile_str.keys():
        i = i + 1
        s = 0
        h = 0
        d = 0
        for f in outfile_str[name]:
            print(f,' ..reading....')
            fr = open(f,'r')
            all_lines = fr.readlines()
            fr.close()
            
            for line in all_lines[discard_line:]:
                if int(line.split()[1]) == 2:
                    d += 1
                if int(line.split()[1]) == 1:
                    h += 1
                if int(line.split()[1]) == 0:
                    s += 1
        
        sliding[float(name.split('_')[-1])] = s*1.0/(s+d+h)
        hopping[float(name.split('_')[-1])] = h*1.0/(s+d+h)
        D3[float(name.split('_')[-1])]      = d*1.0/(s+d+h)
    
    x = np.array(sorted(sliding.keys()))
    y = [sliding[a] for a in x]
    ax_s.plot(x,y)
    ax_s.set_title('sliding')
    ax_s.set_xlabel('salt concentration(M)')
    ax_s.set_ylabel('probability')
    ax_s.set_ylim(0, 1.0)
    plt.figure(1)
    plt.savefig(figure_path + r'\hist_mode_specific_s.png')
    
    x = np.array(sorted(hopping.keys()))
    y = [hopping[a] for a in x]
    ax_h.plot(x,y)
    ax_h.set_title('hoping')
    ax_h.set_xlabel('salt concentration(M)')
    ax_h.set_ylabel('probability')
    ax_h.set_ylim(0, 0.2)
    plt.figure(2)
    plt.savefig(figure_path + r'\hist_mode_specific_h.png')
    
    x = np.array(sorted(D3.keys()))
    y = [D3[a] for a in x]
    ax_d.plot(x,y)
    ax_d.set_title('3d diffusion')
    ax_d.set_xlabel('salt concentration(M)')
    ax_d.set_ylabel('probability')
    ax_d.set_ylim(0, 1.0)
    plt.figure(3)
    plt.savefig(figure_path + r'\hist_mode_specific_d.png')
    #plt.show()

def diffusion_along_DNA(out_path,figure_path):
    print("\ndiffusion_along_DNA")
    discard_line =  2000
    ss_dna_len = 299
    os.chdir(out_path)
    all_files = os.listdir('.')
    inpfile_str=["interface_xpa_dna"]
    outfile_str={}
    for s in inpfile_str:
        for f in all_files:
            if not( s in f):
                continue
            con = f.split('_')[-2]
            name = s + '_' + con
            if name in outfile_str.keys():
                outfile_str[name].append(f)
            else:
                outfile_str[name] = []
                outfile_str[name].append(f)
    
    con = [float(ele.split('_')[-1]) for ele in list(outfile_str.keys())]
    con = sorted(con)
    names = [ inpfile_str[0] + "_" + str('%.2f'%c) for c in con]
    fig_ii = 0
    
    for name in names:
        d_t = 100
        g_t = 5
        file_pos = []
        files = []
        for f in outfile_str[name]:
            print(f,' ...reading...')
            fr = open(f,'r')
            all_lines = fr.readlines()
            fr.close()
            lig_lines = [ line[9:] for line in all_lines[2+discard_line*4::4]]            

            bp_all = []
            for line in lig_lines:
                if line == "\n":
                    bp_all.append(-99999)
                    continue
                
                bp = []
                for each in line.split():
                    if int(each) > ss_dna_len:
                        t =  int((ss_dna_len*2 - int(each))/3) + 1 
                    else:
                        t =  int(int(each)/3) + 1
                    if not (t in bp):
                        bp.append(t)
                        
                bp_all.append(np.mean(bp))
                
            bp_ind = []  
            ind_pos = {}
            start_flag = False
            for i in range(bp_all.__len__()):
                pos = bp_all[i]
                if pos != -99999:
                    if start_flag:
                        ind_pos[i+1] = pos
                    else:
                        start_flag = True
                        ind_pos[i+1] = pos
                else:
                    if not start_flag:
                        continue
                    if i + g_t + 1 < bp_all.__len__():
                        for ii in range(i+1, i+g_t+2):
                            if bp_all[ii] != -99999:
                                break
                        if ii == i + g_t + 1:
                            start_flag = False
                            if ind_pos.__len__() < d_t:
                                ind_pos = {}
                            else:
                                bp_ind.append(deepcopy(ind_pos))
                                ind_pos = {}
                                continue
                    else:
                        break
                
            if start_flag == True and ind_pos.__len__() > d_t:
                bp_ind.append(deepcopy(ind_pos))     
                
            if bp_ind.__len__() > 0:
                file_pos.append(deepcopy(bp_ind))
                files.append(f)
        
        max_size = -1  
        pos_max = 0
        
        temp_file = ''
        xx = 0
        for bp_ind in file_pos:
            for pos_key in bp_ind:
                if max_size < max(pos_key.values()) - min(pos_key.values()):
                    pos_max = deepcopy(pos_key)
                    max_size = max(pos_key.values()) - min(pos_key.values())
                    temp_file = files[xx]
            
            xx = xx + 1
        
        step = sorted(list(pos_max.keys())) 
        id_array = []
        for s in step:
            id_array.append(pos_max[s])
        fig = plt.figure(fig_ii)
        ax = fig.add_subplot(1,1,1)
        ax.plot(np.array(step),np.array(id_array))
        ax.set_title(temp_file.replace('txt','png'))
        ax.set_xlabel('MD step(10^4)')
        ax.set_ylabel('DNA bp index')
        plt.savefig(figure_path + "\\" + temp_file.replace('txt','png') )
        fig_ii = fig_ii + 1
        
    #plt.show()
#<!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#only certain snapshots is considered, where the contact between DNA and protein is formed at certain cumstance. .


def hist_rate_bound(out_path,figure_path):
    print("\hist_rate_bound")
    d_t = 5
    discard_line =  2000
    os.chdir(out_path)
    all_files = os.listdir('.')
    inpfile_str=["contact_xpa_dna"]
    outfile_str={}
    for s in inpfile_str:
        for f in all_files:
            if not( s in f):
                continue
            con = f.split('_')[-2]
            name = s + '_' + con
            if name in outfile_str.keys():
                outfile_str[name].append(f)
            else:
                outfile_str[name] = []
                outfile_str[name].append(f)
             
    i = 0  
    time_bound = {}
    for name in outfile_str.keys():
        i = i + 1 
        T = []
        for f in outfile_str[name]:
            print(f,' ...reading...')
            fr = open(f,'r')
            all_lines = fr.readlines()
            fr.close()
            dis = [ float(line.split()[1]) for line in all_lines[discard_line:]] 
            
            unbound_flag = False
            start_flag = False
            for ind_ in range(len(dis)):
                if dis[ind_] < 10 and not start_flag:
                    continue
                else:
                    start_flag = True
                
                if dis[ind_] == 0 and not unbound_flag:
                    print('unbound  ',ind_)
                    unbound_flag = True  
                    ii_un = ind_
                
                if dis[ind_] > 10 and unbound_flag :
                    
                    if ind_ - ii_un > d_t:
                        print('bound ',ind_)
                        T.append(ind_ - ii_un)
                    unbound_flag = False
        
        if T.__len__() == 0:
            print(name,'  is empty')
        else:
            time_bound[float(name.split('_')[-1])] = np.mean(T)
    rate = []
    for k in sorted(list(time_bound.keys())):
        rate.append(1.0/time_bound[k]) 
    print(sorted(list(time_bound.keys())))
    print(rate)
    
cwd = os.getcwd()
read_path = '/scratch/wfli1/wwzhang/data/XPA/Full/outputfile'
out_path  = r'D:\temp\XPA\XPA_Truncated_out'
figure_path = r'D:\temp\XPA\truncated_figure'

if False:
    os.chdir(cwd)
    dcd_dna_curvature(read_path,out_path)
    sys.exit()
    
    os.chdir(cwd)
    dcd_contact_number_n(read_path,out_path)
    
    os.chdir(cwd)
    dcd_contact_number_dbd(read_path,out_path)
    
    os.chdir(cwd)
    dcd_contact_number_c(read_path,out_path)
    
    os.chdir(cwd)
    dcd_search_mode_dna_protein(read_path,out_path)
    
    os.chdir(cwd)
    dcd_rotation_dna_protein(read_path,out_path)
    
    os.chdir(cwd)
    dcd_distances_com(read_path,out_path)
    
    os.chdir(cwd)
    dcd_distance_rec_lig_n(read_path,out_path)
    
    os.chdir(cwd)
    dcd_distance_rec_lig_dbd(read_path,out_path)
    
    os.chdir(cwd)
    dcd_distance_rec_lig_c(read_path,out_path)
    
    os.chdir(cwd)
    dcd_interface_dna_protein_full(read_path,out_path)
    
    os.chdir(cwd)
    dcd_interface_dna_protein_truncated(read_path,out_path)
    
if True:
    hist_rate_bound(out_path,figure_path)
    sys.exit()
    #hist_mode_specific(out_path,figure_path)
    #sys.exit()
    
    hist_rotation_specific(out_path,figure_path)
    sys.exit()
    
    diffusion_1D_specific(out_path,figure_path)
    sys.exit()
    
    #hist_rotation(out_path,figure_path)
    #sys.exit()
    
    diffusion_1D(out_path,figure_path)
    sys.exit()
    
    #hist_distance_rec_lig(out_path,figure_path)
    #sys.exit()
    
    #traj_distance_rec_lig(out_path,figure_path)
    #sys.exit()
    
    #hist_residue_interface_dna(out_path,figure_path)
    #sys.exit()
    
    hist_distance_rec_lig_98_105(out_path,figure_path)
    sys.exit()
    
    diffusion_along_DNA(out_path,figure_path)
    sys.exit()
 

    
