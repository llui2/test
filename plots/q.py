import matplotlib.pyplot as plt
import numpy as np

plt.rc('font', family='Times')
plt.rc('mathtext', fontset='cm')

# Read data from file
data = np.loadtxt('results/data/q.dat')

if data.ndim == 1:       
       data = np.array([data])

# Set up arrays
N = data[:, 0]
TEMP = data[:, 1]
H = data[:, 2]
p = data[:, 3]
q = data[:, 4]

# Check unique values of H
TEMP_values = np.unique(TEMP)
N_values = np.unique(N)

# Set p value
p_value = 0.5

# Dictionary of colors for each p value
# colors = {40: '#E50D00', 60: '#E6038A', 70: '#AF06E7', 80: '#1F09E8', 100: '#0D87E9', 120: '#10EABC', 140: '#14EB31'}

for i in range(len(TEMP_values)):

       p_H = [[] for x in range(len(N_values))]
       p_q = [[] for x in range(len(N_values))]
       p_N = [[] for x in range(len(N_values))]

       for j in range(len(H)):
              for k in range(len(N_values)):
                     if TEMP[j] == TEMP_values[i] and N[j] == N_values[k] and p[j] == p_value:
                            p_H[k].append(H[j])
                            p_q[k].append(q[j])
                            p_N[k] = N_values[k]

       # Set up figure and subplots
       fig = plt.figure(figsize=(5,5))
       ax = plt.subplot()
       plt.subplots_adjust(wspace=0,hspace=0,left=0.1,top=0.85,right=0.9,bottom=0.1)

       ax.tick_params(direction='in', top=True, right=True)

       # Plot
       plt.title(f"$p={p_value} \\quad T={TEMP_values[i]}$", loc = 'right')

       for j in range(len(N_values)):
              ax.plot(p_H[j], p_q[j], linestyle='-',linewidth=1, marker='o',markersize=2,markeredgewidth=1, label=f'N = {int(N_values[j]):3d}') #, color=colors[N_values[j]]
       ax.set_xlabel('$\Gamma$', fontfamily='Times')
       ax.set_ylabel('$\\langle |q| \\rangle$', fontfamily='Times')
       ax.set_xticks(np.arange(0, 5.5, step=0.5))
       ax.set_yticks(np.arange(0, 1.2, step=0.2))
       ax.set_xlim(0, 5)
       ax.set_ylim(0, 1)

       plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.18), fancybox=False,shadow=False,frameon=False, ncol=4)

       plt.margins(0,10)

       dummy = f'{TEMP_values[i]:.2f}'

       # Show plot
       plt.savefig(f'results/q_T{dummy[:1]+dummy[2:]}.pdf')
