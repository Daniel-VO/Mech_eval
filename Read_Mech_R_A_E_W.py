"""
Created 22. November by Daniel Van Opdenbosch, Technical University of Munich

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed without any warranty or implied warranty of merchantability or fitness for a particular purpose. See the GNU general public license for more details: <http://www.gnu.org/licenses/>
"""

import numpy
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import signal
from scipy import stats
import quantities as pq
from quantities import UncertainQuantity as uq

os.system('mv Read.log Read.alt')

data=numpy.load('data.npy')
#[filename,R,E,A,W,Ag,At,Wt]

nameslist=[]
for n,value in enumerate(data):
	# ~ print(value[0])
	nameslist.append(value[0].split('-')[0]+value[0].split('-')[1])
	# ~ print(value[0])

def common_elements(list1, list2):
	return list(set(list1) & set(list2))

samples=sorted(common_elements(nameslist,nameslist))

nameplot=[]
Rplot=[]
Eplot=[]
Aplot=[]
Wplot=[]
Agplot=[]
Atplot=[]
Wtplot=[]
for name in samples:
	R=numpy.array([])
	E=numpy.array([])
	A=numpy.array([])
	W=numpy.array([])
	Ag=numpy.array([])
	At=numpy.array([])
	Wt=numpy.array([])
	for n,value in enumerate(data):
		if value[0].split('-')[0]+value[0].split('-')[1]==name:
			R=numpy.append(R,float(value[1]))
			E=numpy.append(E,float(value[2]))
			A=numpy.append(A,float(value[3]))
			W=numpy.append(W,float(value[4]))
			Ag=numpy.append(Ag,float(value[5]))
			At=numpy.append(Ag,float(value[6]))
			Wt=numpy.append(Wt,float(value[7]))
	nameplot.append(name.replace('_','-'))
	Rplot.append(uq(numpy.average(R),pq.Pa,numpy.std(R)))
	Eplot.append(uq(numpy.average(E),pq.Pa,numpy.std(E)))
	Aplot.append(uq(numpy.average(A),pq.dimensionless,numpy.std(A)))
	Wplot.append(uq(numpy.average(W),pq.J/pq.m**3,numpy.std(W)))
	Agplot.append(uq(numpy.average(Ag),pq.dimensionless,numpy.std(Ag)))
	Atplot.append(uq(numpy.average(At),pq.dimensionless,numpy.std(At)))
	Wtplot.append(uq(numpy.average(Wt),pq.J/pq.m**3,numpy.std(Wt)))

with open('Read.log','a') as e:
	for s,values in enumerate(nameplot):
		e.write(str([i for i in [nameplot[s],'R:',Rplot[s],'E:',Eplot[s],'A:',Aplot[s],'W:',Wplot[s],'Ag:',Agplot[s],'At:',Atplot[s],'Wt:',Wtplot[s]]]).replace('UncertainQuantity','').replace('array','').replace('\n','').replace('[','').replace(']','').replace('(','').replace(')','').replace("'","")+'\n')

Rplot=uq([i.magnitude for i in Rplot],Rplot[0].units,[i.uncertainty for i in Rplot])
Eplot=uq([i.magnitude for i in Eplot],Eplot[0].units,[i.uncertainty for i in Eplot])
Aplot=uq([i.magnitude for i in Aplot],Aplot[0].units,[i.uncertainty for i in Aplot])
Wplot=uq([i.magnitude for i in Wplot],Wplot[0].units,[i.uncertainty for i in Wplot])
Agplot=uq([i.magnitude for i in Agplot],Agplot[0].units,[i.uncertainty for i in Agplot])
Atplot=uq([i.magnitude for i in Atplot],Atplot[0].units,[i.uncertainty for i in Atplot])
Wtplot=uq([i.magnitude for i in Wtplot],Wtplot[0].units,[i.uncertainty for i in Wtplot])

labels=nameplot
x=numpy.arange(len(labels))
xlabel=r'$\rm{Probe}$'

########################################################################

plt.clf()
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
fig,ax1=plt.subplots(figsize=(7.5/2.54,5.3/2.54))

ax1.errorbar(nameplot,Rplot.magnitude,marker='s',color='k',yerr=numpy.array(Rplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)

# ~ ax1.set_xlabel(xlabel,fontsize=10)
# ~ plt.xticks(x,labels)
# ~ ax1.set_yscale('log')
plt.setp(ax1.xaxis.get_majorticklabels(),rotation=45,ha='right',rotation_mode='anchor')
ax1.set_ylabel(r'$R/\rm{Pa}$',fontsize=10)
ax1.tick_params(direction='out')
ax1.tick_params(axis='x',pad=2,labelsize=8)
ax1.tick_params(axis='y',pad=2,labelsize=8)
ax1.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
ax1.xaxis.get_offset_text().set_size(8)
ax1.yaxis.get_offset_text().set_size(8)
plt.tight_layout(pad=0.5)
plt.savefig('R.pdf',transparent=True)
plt.savefig('R.png',dpi=600)
plt.close('all')


########################################################################

plt.clf()
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
fig,ax1=plt.subplots(figsize=(7.5/2.54,5.3/2.54))

ax1.errorbar(nameplot,Eplot.magnitude,marker='s',color='k',yerr=numpy.array(Eplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)

# ~ ax1.set_xlabel(xlabel,fontsize=10)
# ~ plt.xticks(x,labels)
# ~ ax1.set_yscale('log')
plt.setp(ax1.xaxis.get_majorticklabels(),rotation=45,ha='right',rotation_mode='anchor')
ax1.set_ylabel(r'$E/\rm{Pa}$',fontsize=10)
ax1.tick_params(direction='out')
ax1.tick_params(axis='x',pad=2,labelsize=8)
ax1.tick_params(axis='y',pad=2,labelsize=8)
ax1.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
ax1.xaxis.get_offset_text().set_size(8)
ax1.yaxis.get_offset_text().set_size(8)
plt.tight_layout(pad=0.5)
plt.savefig('E.pdf',transparent=True)
plt.savefig('E.png',dpi=600)
plt.close('all')

########################################################################

plt.clf()
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
fig,ax1=plt.subplots(figsize=(7.5/2.54,5.3/2.54))

ax1.errorbar(nameplot,Aplot.magnitude,marker='s',color='k',yerr=numpy.array(Aplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)
ax1.errorbar(nameplot,Agplot.magnitude,marker='o',color='k',yerr=numpy.array(Agplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)
ax1.errorbar(nameplot,Atplot.magnitude,marker='^',color='k',yerr=numpy.array(Atplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)

# ~ ax1.set_xlabel(xlabel,fontsize=10)
# ~ plt.xticks(x,labels)
# ~ ax1.set_yscale('log')
plt.setp(ax1.xaxis.get_majorticklabels(),rotation=45,ha='right',rotation_mode='anchor')
ax1.set_ylabel(r'$A/1$',fontsize=10)
ax1.tick_params(direction='out')
ax1.tick_params(axis='x',pad=2,labelsize=8)
ax1.tick_params(axis='y',pad=2,labelsize=8)
ax1.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
ax1.xaxis.get_offset_text().set_size(8)
ax1.yaxis.get_offset_text().set_size(8)
plt.tight_layout(pad=0.5)
plt.savefig('A.pdf',transparent=True)
plt.savefig('A.png',dpi=600)
plt.close('all')

########################################################################

plt.clf()
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
fig,ax1=plt.subplots(figsize=(7.5/2.54,5.3/2.54))

ax1.errorbar(nameplot,Wplot.magnitude,marker='s',color='k',yerr=numpy.array(Wplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)
ax1.errorbar(nameplot,Wtplot.magnitude,marker='^',color='k',yerr=numpy.array(Wtplot.uncertainty),markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)

# ~ ax1.set_xlabel(xlabel,fontsize=10)
# ~ plt.xticks(x,labels)
# ~ ax1.set_yscale('log')
plt.setp(ax1.xaxis.get_majorticklabels(),rotation=45,ha='right',rotation_mode='anchor')
ax1.set_ylabel(r'$W/(\rm{J\cdot m}^{-3})$',fontsize=10)
ax1.tick_params(direction='out')
ax1.tick_params(axis='x',pad=2,labelsize=8)
ax1.tick_params(axis='y',pad=2,labelsize=8)
ax1.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
ax1.xaxis.get_offset_text().set_size(8)
ax1.yaxis.get_offset_text().set_size(8)
plt.tight_layout(pad=0.5)
plt.savefig('W.pdf',transparent=True)
plt.savefig('W.png',dpi=600)
plt.close('all')
