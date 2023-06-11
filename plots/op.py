import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

plt.rc('font', family='Times')
plt.rc('mathtext', fontset='cm')

# Read 'input.txt'
with open("input.txt") as f:
    f.readline()
    N = int(f.readline().split(',')[0]) #line 2
    f.readline()
    R = int(f.readline()) #line 4
    for i in range(2):
        f.readline()
    TEMP_list = f.readline().split() #line 7
    for i in range(2):
        f.readline()
    H_list = f.readline().split() #line 10
    for i in range(2):
        f.readline()
    p_list = f.readline().split() #line 13
    f.readline()
    C = int(f.readline()) #line 15
    f.readline()
    MCINI = int(f.readline()) #line 17
    f.readline()
    NSEEDS = int(f.readline()) #line 19
    f.readline()
    SC = int(f.readline()) #line 21

MCTOT = int(MCINI + C*SC/2)

# Dictionary of colors for each p value
colors = {'0.00': '#E50D00', '0.15': '#E6038A', '0.20': '#AF06E7', '0.50': '#1F09E8', '0.7': '#0D87E9', '0.85': '#10EABC', '1.00': '#14EB31'}

os.system(f"mkdir -p results/data/")
f = open("results/data/q.dat", "w")
f.write('#N     TEMP    H       p       <|q|>\n')

q_id = [[[[] for i in range(len(p_list))] for i in range(len(H_list))] for i in range(len(TEMP_list))] # Order parameter integrating disorder
q_id_ini = [[[[] for i in range(len(p_list))] for i in range(len(H_list))] for i in range(len(TEMP_list))] # Order parameter integrating disorder after MCINI
q_avg = [[[[] for i in range(len(p_list))] for i in range(len(H_list))] for i in range(len(TEMP_list))] # Average order parameter integrating disorder after MCINI

q_id_self = [[[[] for i in range(len(p_list))] for i in range(len(H_list))] for i in range(len(TEMP_list))] # Order parameter integrating disorder

for i in range(len(TEMP_list)):
    for j in range(len(H_list)):

        TEMP = TEMP_list[i][:1] + TEMP_list[i][2:]
        H = H_list[j][:1] + H_list[j][2:]

        for k in range(len(p_list)):

            p = p_list[k][:1] + p_list[k][2:]

            for SEED in range(100,100+NSEEDS-1):

                if SEED == 100:
                    data = np.loadtxt(f"results/order/T{TEMP}_Γ{H}/q_{p}_{SEED}.dat") 
                    MCS = data[:,0]
                    q = data[:,1]
                    q_self = data[:,2]
                else:
                    data = np.loadtxt(f"results/order/T{TEMP}_Γ{H}/q_{p}_{SEED}.dat") 
                    q += data[:,1]
                    q_self += data[:,2]

            q_id[i][j][k] = q/NSEEDS
            q_id_self[i][j][k] = q_self/NSEEDS
            index = np.where(MCS>=MCINI)[0][0]
            q_id_ini[i][j][k] = q[index:]/NSEEDS
            q_avg[i][j][k] = np.average(q_id_ini[i][j][k])
            f.write(f"{N:3d}    {TEMP_list[i]}    {H_list[j]}    {p_list[k]}    {q_avg[i][j][k]} \n")

        # Plot time series
        fig = plt.figure(figsize=(7,4))
        ax = plt.subplot()
        ax1 = ax.twinx()
        plt.subplots_adjust(wspace=0,hspace=0,left=0.1,top=0.85,right=0.9,bottom=0.1)

        ax.tick_params(direction='in', top=True, right=True)

        plt.title(f"$T={TEMP_list[i]}$  $\\Gamma={H_list[j]}$", loc = 'right')
        for k in range(len(p_list)):
            ax.plot(MCS, q_id[i][j][k], linestyle='-',linewidth=1, marker=' ',markersize=1, color=colors[p_list[k]], label=f'p = {p_list[k]}')
            ax.plot(MCS, q_id_self[i][j][k], linestyle='-',linewidth=1, marker=' ',markersize=1)
            ax.axhline(y = q_avg[i][j][k], linestyle = '--', alpha=0.5, color=colors[p_list[k]])
            ax.axvline(x = MCINI, linestyle = '-',alpha=0.2, color='black')
            ax1.plot(np.NaN, np.NaN, ls='--',alpha=0.5, color=colors[p_list[k]], label = f'$\\langle |q| \\rangle = {q_avg[i][j][k]:.3f}$')
            ax1.get_yaxis().set_visible(False)
            ax.set_xlim(0,MCTOT)

        ax.set_xlabel('MCS', fontfamily='Times')
        ax.set_ylabel('$|q|$', fontfamily='Times')
        ax.set_xlim(0,)
        ax.set_ylim(0,1)

        ax1.legend(loc=1, fancybox=False,shadow=False,frameon=False, ncol=4)

        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.18), fancybox=False,shadow=False,frameon=False, ncol=4)
        plt.savefig(f"results/op_T{TEMP}_Γ{H}.pdf")

    # Plot order parameter mean value vs H
    fig = plt.figure(figsize=(5,5))
    ax2 = plt.subplot()
    plt.subplots_adjust(wspace=0,hspace=0,left=0.1,top=0.85,right=0.9,bottom=0.1)

    ax2.tick_params(direction='in', top=True, right=True)

    plt.title(f"$T={TEMP_list[i]}$", loc = 'right')

    H_values = [float(x) for x in H_list]
    q_values = []

    for k in range(len(p_list)):
        for x in range(len(H_values)):
            q_values.append(q_avg[i][x][k])
        ax2.plot(H_values, q_values, linestyle='-',linewidth=1, marker=' ',markersize=1, color=colors[p_list[k]], label=f'p = {p_list[k]}')

    ax2.set_xlabel('$\\Gamma$', fontfamily='Times')
    ax2.set_ylabel('$\\langle |q| \\rangle$', fontfamily='Times')
    ax2.set_xlim(0,6)
    ax2.set_ylim(0,1)

    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.18), fancybox=False,shadow=False,frameon=False, ncol=4)        
    plt.savefig(f"results/q_mean_T{TEMP}.pdf")