
#Modelling_Nucleus_Accumbens
#3D view of the neurons and the electrodes considered in the computational model.

from brian2 import *
start_scope()

import numpy as np



from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


number_of_neurons_in_msnd1_core=100
number_of_neurons_in_msnd1_shell=100
xN_msnd1_core=5*rand(number_of_neurons_in_msnd1_core)
yN_msnd1_core=14*rand(number_of_neurons_in_msnd1_core)
zN_msnd1_core=7*rand(number_of_neurons_in_msnd1_core)
#xN_msnd2_core=5*rand(number_of_neurons_in_msnd2_core)
#yN_msnd2_core=10*rand(number_of_neurons_in_msnd2_core)
#zN_msnd2_core=10*rand(number_of_neurons_in_msnd2_core)
xN_msnd1_shell=5+5*rand(number_of_neurons_in_msnd1_shell)
yN_msnd1_shell=14*rand(number_of_neurons_in_msnd1_shell)
zN_msnd1_shell=7*rand(number_of_neurons_in_msnd1_shell)
#xN_msnd2_shell=5+5*rand(number_of_neurons_in_msnd2_shell)
#yN_msnd2_shell=10*rand(number_of_neurons_in_msnd2_shell)
#zN_msnd2_shell=10*rand(number_of_neurons_in_msnd2_shell)
#xN_nacc_in=10*rand(number_of_neurons_in_nacc_in)
#yN_nacc_in=10*rand(number_of_neurons_in_nacc_in)
#zN_nacc_in=10*rand(number_of_neurons_in_nacc_in)


array_dim=1



xE=[3]
yE=[5]
zE=[3+rand()]
for i in range(5):
    for j in range(5):
        xE.append(xE[0]+j*array_dim)
        yE.append(yE[0]+i*array_dim)   
        zE.append(zE[0])







print(xE)

print(yE)

print(zE)


#xE11=5*rand()
#yE11=10*rand()
#zE11=10*rand()
#
#
#
#
#
#
#
#xE21=5+5*rand()
#yE21=10*rand()
#zE21=10*rand()



#zline = zN_msnd1_core
#xline = xN_msnd1_core
#yline = yN_msnd1_core
#ax.plot3D(xline, yline, zline, 'gray')


# Data for a three-dimensional line
#zline = np.linspace(0, 15, 1000)
#xline = np.sin(zline)
#yline = np.cos(zline)
#ax.plot3D(xline, yline, zline, 'gray')

# Data for three-dimensional scattered points


xline=[[xE[0],xE[0]]]
yline=[[yE[0],yE[0]]]
zline=[[zE[0],xE[0]+10]]



for i in range(25):
    xline.append([xE[i+1],xE[i+1]])
    yline.append([yE[i+1],yE[i+1]])
    zline.append([zE[i+1],zE[i+1]+10])







fig = plt.figure()
plt.subplots_adjust(top = 0.98, bottom = 0.02, left = 0.02, right = 0.98, hspace = .1, wspace=0.06)

ax1 = plt.axes(projection='3d')
ax1.scatter3D(xN_msnd1_shell, yN_msnd1_shell, zN_msnd1_shell, c='Blue',s=100, label="Core");#xN_msnd1_shell, cmap='Blues');
ax1.scatter3D(xN_msnd1_core, yN_msnd1_core, zN_msnd1_core, c='Green',s=100, label="Shell");#xN_msnd1_core, cmap='Greens');
ax1.scatter3D(xE, yE, zE, c='red',s=100, label="Electrot");#, cmap='Reds');
for i in range(26):
    ax1.plot3D(xline[i], yline[i], zline[i], 'gray')
#ax.plot3D(xline5, yline5, zline5, 'gray')
#ax.plot3D(xline6, yline6, zline6, 'gray')
#ax.plot3D(xline10, yline10, zline10, 'gray')
#ax.plot3D(xline25, yline25, zline25, 'gray')

#ax.scatter3D(xE11, yE11, zE11, c='red');#, cmap='Reds');
#ax.scatter3D(xE21, yE21, zE21, c='red');#, cmap='Reds');





ax1.grid(True)
#ax1.view_init(elev=0., azim=90) ##### You can choose this view or 

ax1.view_init(elev=45., azim=45)   ##### You can choose this view
plt.legend(loc='best',prop={'family': 'Times New Roman', 'size':'25'},frameon=False)#,fontsize=50)

ax1.set_xlabel('X',visible=0)
ax1.set_ylabel('Y',visible=0)
ax1.set_zlabel('Z',visible=0)




figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
figurre = plt.gcf() # get current figure
figurre.set_size_inches(18, 10)
plt.savefig('3dview.pdf',dpi=600)

plt.show()


