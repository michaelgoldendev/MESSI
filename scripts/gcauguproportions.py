#!/usr/bin/env python
# a stacked bar plot with errorbars
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import lines

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.figure(figsize=(6, 4))

mat = np.ones((6,14))/6.0

names = ["WDF (r)", "TYLCV (r)", "MSV (r)", "Bocavirus (r)", "BCTV (r)", "Tobamovirus (r)", "Rhinovirus A (r)", "H. poliovirus 1 (r)", "Hepatitis A (r)", "FMDV (r)", "RF01846 (r)", "RF00379 (r)", "RF00010 (r)", "RF00003 (r)", "RF00001 (r)"]
mat = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.628499, 0.371501, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.35533, 0.0, 0.64467, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
medians = [[5.20114, 3.29335, 1.42303], [3.44638, 3.3834, 1.5268], [5.0632, 3.72829, 1.40195], [4.17431, 2.57663, 1.05141], [4.35757, 2.80693, 1.08446], [2.94822, 2.45572, 2.23185], [10.4561, 5.41171, 8.97726], [3.46465, 3.07641, 3.52301], [3.66289, 1.90069, 1.30974], [3.95973, 3.68022, 2.44305], [4.53752, 2.92626, 2.16112], [7.51311, 3.67866, 1.92968], [5.51895, 5.16185, 2.33751], [5.89402, 3.6172, 2.47293], [6.6818, 4.74913, 2.15773]]

#names = ["MSV (r)", "TYLCV (r)", "BCTV (r)", "Bocavirus (r)", "WDF (r)", "Hepatitis A (r)", "Rhinovirus A (r)", "Tobamovirus (r)", "H. poliovirus (r)", "FMDV (r)", "RF01846 (r)", "RF00379 (r)", "RF00010 (r)", "RF00003 (r)", "RF00001 (r)"]
#mat = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.628499, 0.371501, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.35533, 0.0, 0.64467, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
#medians = [[5.0632, 3.72829, 1.40195], [3.44638, 3.3834, 1.5268], [4.35757, 2.80693, 1.08446], [4.17431, 2.57663, 1.05141], [5.20114, 3.29335, 1.42303], [3.66289, 1.90069, 1.30974], [10.4561, 5.41171, 8.97726], [2.94822, 2.45572, 2.23185], [3.46465, 3.07641, 3.52301], [3.95973, 3.68022, 2.44305], [4.53752, 2.92626, 2.16112], [7.51311, 3.67866, 1.92968], [5.51895, 5.16185, 2.33751], [5.89402, 3.6172, 2.47293], [6.6818, 4.74913, 2.15773]]
#mat = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.998278, 0.00172216, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.998862, 0.00113766, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0, 0.0, 0.0], [0.451993, 0.0, 0.548007, 0.0, 0.0, 0.0], [0.0214693, 0.000167729, 0.814995, 0.0, 0.1632, 0.000167729], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
#mat = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]]


mat = np.transpose(np.array(mat))

#medians = [[9.4269, 6.19605, 1.8343], [8.74215, 5.41848, 1.54507], [5.63509, 4.48039, 1.17535], [5.16661, 3.03902, 1.12004], [7.16224, 5.76437, 1.74421], [6.95397, 4.79819, 2.30639], [8.30365, 4.50197, 3.14475], [9.14189, 7.30034, 3.05268], [13.1789, 7.31558, 3.57201], [7.78112, 5.0996, 2.68072], [10.2107, 3.5208, 1.00776], [22.6621, 12.7565, 20.4124], [6.8476, 5.6736, 5.61306], [6.25416, 5.29102, 5.9331], [8.11754, 7.14594, 4.61971]]
#medians = [[5.20114, 3.29335, 1.42303]]

#names = ["WDV","Bocavirus","BCTV","TYLCV","MSV","RF00001","RF00003","RF00010","RF00379","RF01846","Hepatitis A","Rhinovirus A","Tobamovirus","Human poliovirus 1","FMDV"]
#names = ["WDV (r)"]



colors = ['#60BD68','#FAA43A','#4D4D4D','#5DA5DA','#F15854','#DECF3F']
labels = [r'$\lambda_\textrm{GC} > \lambda_\textrm{AU} > \lambda_\textrm{GU}$', 
                     r'$\lambda_\textrm{AU} > \lambda_\textrm{GC} > \lambda_\textrm{GU}$',
                     r'$\lambda_\textrm{GC} > \lambda_\textrm{GU} > \lambda_\textrm{AU}$',
                     r'$\lambda_\textrm{AU} > \lambda_\textrm{GU} > \lambda_\textrm{GC}$',
                     r'$\lambda_\textrm{GU} > \lambda_\textrm{GC} > \lambda_\textrm{AU}$',
                     r'$\lambda_\textrm{GU} > \lambda_\textrm{AU} > \lambda_\textrm{GC}$']






for row in range(1, mat.shape[0]):
    mat[mat.shape[0]-row-1,:] += mat[mat.shape[0]-row,:]

plt.text(1.17, mat.shape[1]+1.1, "Medians", weight='bold', horizontalalignment='center',verticalalignment='center',color='black',fontsize=10)
plt.text(1.07, mat.shape[1], r'$\lambda_\textrm{GC}$', weight='bold', horizontalalignment='center',verticalalignment='center',color='black',fontsize=10)
plt.text(1.17, mat.shape[1], r'$\lambda_\textrm{AU}$', weight='bold', horizontalalignment='center',verticalalignment='center',color='black',fontsize=10)
plt.text(1.27, mat.shape[1], r'$\lambda_\textrm{GU}$', weight='bold', horizontalalignment='center',verticalalignment='center',color='black',fontsize=10)
for row in range(0, mat.shape[1]):
    for col in range(0, mat.shape[0]):
        v = mat[col,row]
        if col < 5:
            v = mat[col,row]-mat[col+1,row]
        if v > 0.05:
            fcolor = 'black'
            if col == 2:
                fcolor = 'white'
            plt.text(mat[col,row]-v/2, row-0.08, "%0.2f" % v, horizontalalignment='center',verticalalignment='center',color=fcolor,fontsize=10)
    plt.text(1.07, row-0.08, "%0.2f" % medians[row][0], horizontalalignment='center',verticalalignment='center',color='black',fontsize=10)
    plt.text(1.17, row-0.08, "%0.2f" % medians[row][1], horizontalalignment='center',verticalalignment='center',color='black',fontsize=10)
    plt.text(1.27, row-0.08, "%0.2f" % medians[row][2], horizontalalignment='center',verticalalignment='center',color='black',fontsize=10)
    #plt.plot([0.9, 1.08], [row,row], color="black", linewidth=1.0)
    #line = lines.Line2D([1.0, 1.08], [row,row], lw=1.0, color='black')
    #ax2.add_line(line)
	
    
N = mat.shape[1]
ind = np.arange(N)    # the x locations for the groups
width = 0.6       # the width of the bars: can also be len(x) sequence

pltrefs = []
for row in range(0, mat.shape[0]):
    p = plt.barh(ind, mat[row,:], width, color=colors[row])
    pltrefs.append(p)


plt.title(r'Relative magnitudes of $\lambda_{\textrm{GC}}$, $\lambda_{\textrm{AU}}$, and $\lambda_{\textrm{GU}}$', fontsize=14)
plt.xlabel('Posterior probabilities', fontsize=13)
#plt.xticks(ind+width/2., ('G1', 'G2', 'G3', 'G4', 'G5') )
plt.ylabel('Datasets', fontsize=13)
plt.yticks(ind, names)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.xlim(0.0,1.0)


#plt.legend(pltrefs, labels, loc='upper center', bbox_to_anchor=(0.5, -0.15),  ncol=2)
plt.legend(pltrefs, labels, loc='center left', bbox_to_anchor=(1.32, 0.5),  ncol=1)

plt.show()