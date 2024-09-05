"""
Created 05. September 2024 by Daniel Van Opdenbosch, Technical University of Munich

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed without any warranty or implied warranty of merchantability or fitness for a particular purpose. See the GNU general public license for more details: <http://www.gnu.org/licenses/>
"""

import os
import numpy as np
import quantities as pq
from quantities import UncertainQuantity as uq
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy import stats

def namesplitter(d):
	return d[0].split('-')[0]+d[0].split('-')[1]

os.system('mv read.log read.alt')

data=np.load('mech.npy',allow_pickle=True)
#[filename,R,E,A,W,Rp,Ag,At,Wt]

nameslist=[]
for d in data:
	# ~ print(d[0])
	nameslist.append(namesplitter(d))
	# ~ print(d[0])

nameplot=[]
Rplot=[]
Eplot=[]
Aplot=[]
Wplot=[]
Rpplot=[]
Agplot=[]
Atplot=[]
Wtplot=[]
for name in np.unique(nameslist):
	R=[float(d[1]) for d in data if namesplitter(d)==name]
	E=[float(d[2]) for d in data if namesplitter(d)==name]
	A=[float(d[3]) for d in data if namesplitter(d)==name]
	W=[float(d[4]) for d in data if namesplitter(d)==name]
	Rp=[float(d[5]) for d in data if namesplitter(d)==name]
	Ag=[float(d[6]) for d in data if namesplitter(d)==name]
	At=[float(d[7]) for d in data if namesplitter(d)==name]
	Wt=[float(d[8]) for d in data if namesplitter(d)==name]
	nameplot.append(name.replace('_','-'))
	Rplot.append(uq(np.median(R),pq.Pa,stats.median_abs_deviation(R)))
	Eplot.append(uq(np.median(E),pq.Pa,stats.median_abs_deviation(E)))
	Aplot.append(uq(np.median(A),pq.dimensionless,stats.median_abs_deviation(A)))
	Wplot.append(uq(np.median(W),pq.J/pq.m**3,stats.median_abs_deviation(W)))
	Rpplot.append(uq(np.median(Rp),pq.Pa,stats.median_abs_deviation(Rp)))
	Agplot.append(uq(np.median(Ag),pq.dimensionless,stats.median_abs_deviation(Ag)))
	Atplot.append(uq(np.median(At),pq.dimensionless,stats.median_abs_deviation(At)))
	Wtplot.append(uq(np.median(Wt),pq.J/pq.m**3,stats.median_abs_deviation(Wt)))

for s,values in enumerate(nameplot):
	print(str([i for i in [nameplot[s],'R:',Rplot[s],'E:',Eplot[s],'A:',Aplot[s],'W:',Wplot[s],'Rp:',Rpplot[s],'Ag:',Agplot[s],'At:',Atplot[s],'Wt:',Wtplot[s]]]).replace('UncertainQuantity','').replace('array','').replace('\n','').replace('[','').replace(']','').replace('(','').replace(')','').replace("'",""),file=open('read.log','a'))

Rplot=uq([i.magnitude for i in Rplot],Rplot[0].units,[i.uncertainty for i in Rplot])
Eplot=uq([i.magnitude for i in Eplot],Eplot[0].units,[i.uncertainty for i in Eplot])
Aplot=uq([i.magnitude for i in Aplot],Aplot[0].units,[i.uncertainty for i in Aplot])
Wplot=uq([i.magnitude for i in Wplot],Wplot[0].units,[i.uncertainty for i in Wplot])
Rpplot=uq([i.magnitude for i in Rpplot],Rpplot[0].units,[i.uncertainty for i in Rpplot])
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
plt.figure(figsize=(7.5/2.54,5.3/2.54))

plt.errorbar(nameplot,Rplot.magnitude,marker='s',color='k',yerr=np.array(Rplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)
plt.errorbar(nameplot,Rpplot.magnitude,marker='o',color='k',yerr=np.array(Rpplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)

# ~ plt.xlabel(xlabel,fontsize=10)
# ~ plt.xticks(x,labels)
# ~ plt.yscale('log')
plt.setp(plt.gca().xaxis.get_majorticklabels(),rotation=45,ha='right',rotation_mode='anchor')
plt.ylabel(r'$R/\rm{Pa}$',fontsize=10)
plt.tick_params(axis='both',pad=2,labelsize=8)
plt.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
plt.gca().yaxis.get_offset_text().set_size(8)
plt.tight_layout(pad=0.1)
plt.savefig('R.pdf')
plt.savefig('R.png',dpi=300)

####

plt.close('all')
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
plt.figure(figsize=(7.5/2.54,5.3/2.54))

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
plt.savefig('E.pdf')
plt.savefig('E.png',dpi=300)

####

plt.close('all')
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
plt.figure(figsize=(7.5/2.54,5.3/2.54))

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
plt.savefig('A.pdf')
plt.savefig('A.png',dpi=300)

####

plt.close('all')
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
plt.figure(figsize=(7.5/2.54,5.3/2.54))

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
plt.savefig('W.pdf')
plt.savefig('W.png',dpi=300)
