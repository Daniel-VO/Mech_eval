"""
Created 14. November 2023 by Daniel Van Opdenbosch, Technical University of Munich

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed without any warranty or implied warranty of merchantability or fitness for a particular purpose. See the GNU general public license for more details: <http://www.gnu.org/licenses/>
"""

import numpy as np
import glob
import os
import quantities as pq
from quantities import UncertainQuantity as uq
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy import signal
from scipy import stats

def common_elements(list1, list2):
	return list(set(list1) & set(list2))

os.system('mv read.log read.alt')

data=np.load('data.npy',allow_pickle=True)
#[filename,R,E,A,W,Re,Ag,At,Wt]

nameslist=[]
for n,value in enumerate(data):
	# ~ print(value[0])
	nameslist.append(value[0].split('-')[0]+value[0].split('-')[1])
	# ~ print(value[0])

nameplot=[]
Rplot=[]
Eplot=[]
Aplot=[]
Wplot=[]
Replot=[]
Agplot=[]
Atplot=[]
Wtplot=[]
for name in sorted(common_elements(nameslist,nameslist)):
	R=np.array([])
	E=np.array([])
	A=np.array([])
	W=np.array([])
	Re=np.array([])
	Ag=np.array([])
	At=np.array([])
	Wt=np.array([])
	for n,value in enumerate(data):
		if value[0].split('-')[0]+value[0].split('-')[1]==name:
			R=np.append(R,float(value[1]))
			E=np.append(E,float(value[2]))
			A=np.append(A,float(value[3]))
			W=np.append(W,float(value[4]))
			Re=np.append(Re,float(value[5]))
			Ag=np.append(Ag,float(value[6]))
			At=np.append(At,float(value[7]))
			Wt=np.append(Wt,float(value[8]))
	nameplot.append(name.replace('_','-'))
	Rplot.append(uq(np.median(R),pq.Pa,stats.median_abs_deviation(R)))
	Eplot.append(uq(np.median(E),pq.Pa,stats.median_abs_deviation(E)))
	Aplot.append(uq(np.median(A),pq.dimensionless,stats.median_abs_deviation(A)))
	Wplot.append(uq(np.median(W),pq.J/pq.m**3,stats.median_abs_deviation(W)))
	Replot.append(uq(np.median(Re),pq.Pa,stats.median_abs_deviation(Re)))
	Agplot.append(uq(np.median(Ag),pq.dimensionless,stats.median_abs_deviation(Ag)))
	Atplot.append(uq(np.median(At),pq.dimensionless,stats.median_abs_deviation(At)))
	Wtplot.append(uq(np.median(Wt),pq.J/pq.m**3,stats.median_abs_deviation(Wt)))

with open('read.log','a') as e:
	for s,values in enumerate(nameplot):
		e.write(str([i for i in [nameplot[s],'R:',Rplot[s],'E:',Eplot[s],'A:',Aplot[s],'W:',Wplot[s],'Re:',Replot[s],'Ag:',Agplot[s],'At:',Atplot[s],'Wt:',Wtplot[s]]]).replace('UncertainQuantity','').replace('array','').replace('\n','').replace('[','').replace(']','').replace('(','').replace(')','').replace("'","")+'\n')

Rplot=uq([i.magnitude for i in Rplot],Rplot[0].units,[i.uncertainty for i in Rplot])
Eplot=uq([i.magnitude for i in Eplot],Eplot[0].units,[i.uncertainty for i in Eplot])
Aplot=uq([i.magnitude for i in Aplot],Aplot[0].units,[i.uncertainty for i in Aplot])
Wplot=uq([i.magnitude for i in Wplot],Wplot[0].units,[i.uncertainty for i in Wplot])
Replot=uq([i.magnitude for i in Replot],Replot[0].units,[i.uncertainty for i in Replot])
Agplot=uq([i.magnitude for i in Agplot],Agplot[0].units,[i.uncertainty for i in Agplot])
Atplot=uq([i.magnitude for i in Atplot],Atplot[0].units,[i.uncertainty for i in Atplot])
Wtplot=uq([i.magnitude for i in Wtplot],Wtplot[0].units,[i.uncertainty for i in Wtplot])

labels=nameplot
x=np.arange(len(labels))
xlabel=r'$\rm{Probe}$'

####

plt.close('all')
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
plt.subplots(figsize=(7.5/2.54,5.3/2.54))

plt.errorbar(nameplot,Rplot.magnitude,marker='s',color='k',yerr=np.array(Rplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)
plt.errorbar(nameplot,Replot.magnitude,marker='o',color='k',yerr=np.array(Replot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)

# ~ plt.xlabel(xlabel,fontsize=10)
# ~ plt.xticks(x,labels)
# ~ plt.yscale('log')
plt.setp(plt.gca().xaxis.get_majorticklabels(),rotation=45,ha='right',rotation_mode='anchor')
plt.ylabel(r'$R/\rm{Pa}$',fontsize=10)
plt.tick_params(axis='both',pad=2,labelsize=8)
plt.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
plt.gca().yaxis.get_offset_text().set_size(8)
plt.tight_layout(pad=0.1)
plt.savefig('R.pdf',transparent=True)
plt.savefig('R.png',dpi=300)

####

plt.close('all')
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
plt.subplots(figsize=(7.5/2.54,5.3/2.54))

plt.errorbar(nameplot,Eplot.magnitude,marker='s',color='k',yerr=np.array(Eplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)

# ~ plt.xlabel(xlabel,fontsize=10)
# ~ plt.xticks(x,labels)
# ~ plt.yscale('log')
plt.setp(plt.gca().xaxis.get_majorticklabels(),rotation=45,ha='right',rotation_mode='anchor')
plt.ylabel(r'$E/\rm{Pa}$',fontsize=10)
plt.tick_params(axis='both',pad=2,labelsize=8)
plt.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
plt.gca().yaxis.get_offset_text().set_size(8)
plt.tight_layout(pad=0.1)
plt.savefig('E.pdf',transparent=True)
plt.savefig('E.png',dpi=300)

####

plt.close('all')
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
plt.subplots(figsize=(7.5/2.54,5.3/2.54))

plt.errorbar(nameplot,Aplot.magnitude,marker='s',color='k',yerr=np.array(Aplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)
plt.errorbar(nameplot,Agplot.magnitude,marker='o',color='k',yerr=np.array(Agplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)
plt.errorbar(nameplot,Atplot.magnitude,marker='^',color='k',yerr=np.array(Atplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)

# ~ plt.xlabel(xlabel,fontsize=10)
# ~ plt.xticks(x,labels)
# ~ plt.yscale('log')
plt.setp(plt.gca().xaxis.get_majorticklabels(),rotation=45,ha='right',rotation_mode='anchor')
plt.ylabel(r'$A/1$',fontsize=10)
plt.tick_params(axis='both',pad=2,labelsize=8)
plt.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
plt.gca().yaxis.get_offset_text().set_size(8)
plt.tight_layout(pad=0.1)
plt.savefig('A.pdf',transparent=True)
plt.savefig('A.png',dpi=300)

####

plt.close('all')
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
plt.subplots(figsize=(7.5/2.54,5.3/2.54))

plt.errorbar(nameplot,Wplot.magnitude,marker='s',color='k',yerr=np.array(Wplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)
plt.errorbar(nameplot,Wtplot.magnitude,marker='^',color='k',yerr=np.array(Wtplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)

# ~ plt.xlabel(xlabel,fontsize=10)
# ~ plt.xticks(x,labels)
# ~ plt.yscale('log')
plt.setp(plt.gca().xaxis.get_majorticklabels(),rotation=45,ha='right',rotation_mode='anchor')
plt.ylabel(r'$W/(\rm{J\cdot m}^{-3})$',fontsize=10)
plt.tick_params(axis='both',pad=2,labelsize=8)
plt.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
plt.gca().yaxis.get_offset_text().set_size(8)
plt.tight_layout(pad=0.1)
plt.savefig('W.pdf',transparent=True)
plt.savefig('W.png',dpi=300)
