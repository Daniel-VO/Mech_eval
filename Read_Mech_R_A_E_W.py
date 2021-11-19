"""
Created 26. May 2021 by Daniel Van Opdenbosch, Technical University of Munich

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed without any warranty or implied warranty of merchantability or fitness for a particular purpose. See the GNU general public license for more details: <http://www.gnu.org/licenses/>
"""

import numpy
import glob
import os
import fnmatch
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import signal
from scipy import stats

os.system('mv Read.log Read.alt')

data=numpy.load('data.npy')
#[filename,R,A,E,W,Ag,Wt]

nameslist=[]
for n,value in enumerate(data):
	print(value[0])
	nameslist.append(value[0])
	print(value[0])

def common_elements(list1, list2):
	return list(set(list1) & set(list2))

samples=sorted(common_elements(nameslist,nameslist))

sampledata=[]

for m,name in enumerate(samples):
	R=numpy.array([])
	A=numpy.array([])
	E=numpy.array([])
	W=numpy.array([])
	Ag=numpy.array([])
	At=numpy.array([])
	Wt=numpy.array([])
	for n,value in enumerate(data):
		if value[0]==name:
			R=numpy.append(R,float(value[1]))
			A=numpy.append(A,float(value[2]))
			E=numpy.append(E,float(value[3]))
			W=numpy.append(W,float(value[4]))
			Ag=numpy.append(Ag,float(value[5]))
			At=numpy.append(Ag,float(value[6]))
			Wt=numpy.append(Wt,float(value[7]))
	sampledata.append([name.replace('_','-'),'R:',numpy.average(R),'pm',numpy.std(R),'A:',numpy.average(A),'pm',numpy.std(A),'E:',numpy.average(E),'pm',numpy.std(E),'W:',numpy.average(W),'pm',numpy.std(W),'Ag:',numpy.average(Ag),'pm',numpy.std(Ag),'At:',numpy.average(At),'pm',numpy.std(At),'Wt:',numpy.average(Wt),'pm',numpy.std(Wt)])

with open('Read.log','a') as e:
	for s in sampledata:
		e.write(str(s)+'\n')

nameplot=[]
Rplot=[]
Rstdplot=[]
Aplot=[]
Astdplot=[]
Eplot=[]
Estdplot=[]
Econfloplot=[]
Econfupplot=[]
Wplot=[]
Wstdplot=[]
Agplot=[]
Agstdplot=[]
Atplot=[]
Atstdplot=[]
Wtplot=[]
Wtstdplot=[]
for s in samples:	#Probensaetze
	for m,name in enumerate(sampledata):
		if s in sampledata[m][0]:
			nameplot.append(r'$\rm{'+s+'}$')
			Rplot.append(sampledata[m][2])
			Rstdplot.append(sampledata[m][4])
			Aplot.append(sampledata[m][6])
			Astdplot.append(sampledata[m][8])
			Eplot.append(sampledata[m][10])
			Estdplot.append(sampledata[m][12])
			Wplot.append(sampledata[m][14])
			Wstdplot.append(sampledata[m][16])
			Agplot.append(sampledata[m][18])
			Agstdplot.append(sampledata[m][20])
			Atplot.append(sampledata[m][22])
			Atstdplot.append(sampledata[m][24])
			Wtplot.append(sampledata[m][26])
			Wtstdplot.append(sampledata[m][28])
			xlabel=r'$\rm{Probe}$'

########################################################################

plt.clf()
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
fig,ax1=plt.subplots(figsize=(7.5/2.54,5.3/2.54))

ax1.errorbar(nameplot,Rplot,marker='s',color='k',yerr=Rstdplot,markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)

# ~ ax1.set_xlabel(xlabel,fontsize=10)
plt.setp(ax1.xaxis.get_majorticklabels(),rotation=45,ha='right',rotation_mode='anchor')
ax1.set_ylabel(r'$R/\rm{Pa}$',fontsize=10)
ax1.tick_params(direction='out')
ax1.tick_params(axis='x',pad=2,labelsize=8)
ax1.tick_params(axis='y',pad=2,labelsize=8)
# ~ ax1.set_yscale('log')
ax1.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
ax1.xaxis.get_offset_text().set_size(8)
ax1.yaxis.get_offset_text().set_size(8)
plt.tight_layout(pad=0.1)
plt.savefig('R.pdf',transparent=True)
plt.savefig('R.png',dpi=600)
plt.close('all')

########################################################################

plt.clf()
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
fig,ax1=plt.subplots(figsize=(7.5/2.54,5.3/2.54))

ax1.errorbar(nameplot,Aplot,marker='s',color='k',yerr=Astdplot,markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)
ax1.errorbar(nameplot,Agplot,marker='o',color='k',yerr=Agstdplot,markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)
ax1.errorbar(nameplot,Atplot,marker='^',color='k',yerr=Agstdplot,markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)

# ~ ax1.set_xlabel(xlabel,fontsize=10)
plt.setp(ax1.xaxis.get_majorticklabels(),rotation=45,ha='right',rotation_mode='anchor')
ax1.set_ylabel(r'$A/1$',fontsize=10)
ax1.tick_params(direction='out')
ax1.tick_params(axis='x',pad=2,labelsize=8)
ax1.tick_params(axis='y',pad=2,labelsize=8)
# ~ ax1.set_yscale('log')
ax1.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
ax1.xaxis.get_offset_text().set_size(8)
ax1.yaxis.get_offset_text().set_size(8)
plt.tight_layout(pad=0.1)
plt.savefig('A.pdf',transparent=True)
plt.savefig('A.png',dpi=600)
plt.close('all')

########################################################################

plt.clf()
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
fig,ax1=plt.subplots(figsize=(7.5/2.54,5.3/2.54))

ax1.errorbar(nameplot,Eplot,marker='s',color='k',yerr=Estdplot,markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)

# ~ ax1.set_xlabel(xlabel,fontsize=10)
plt.setp(ax1.xaxis.get_majorticklabels(),rotation=45,ha='right',rotation_mode='anchor')
ax1.set_ylabel(r'$E/\rm{Pa}$',fontsize=10)
ax1.tick_params(direction='out')
ax1.tick_params(axis='x',pad=2,labelsize=8)
ax1.tick_params(axis='y',pad=2,labelsize=8)
# ~ ax1.set_yscale('log')
ax1.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
ax1.xaxis.get_offset_text().set_size(8)
ax1.yaxis.get_offset_text().set_size(8)
plt.tight_layout(pad=0.1)
plt.savefig('E.pdf',transparent=True)
plt.savefig('E.png',dpi=600)
plt.close('all')

########################################################################

plt.clf()
mpl.rc('text',usetex=True)
mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
fig,ax1=plt.subplots(figsize=(7.5/2.54,5.3/2.54))

ax1.errorbar(nameplot,Wplot,marker='s',color='k',yerr=Wstdplot,markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)
ax1.errorbar(nameplot,Wtplot,marker='^',color='k',yerr=Wtstdplot,markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0)

# ~ ax1.set_xlabel(xlabel,fontsize=10)
plt.setp(ax1.xaxis.get_majorticklabels(),rotation=45,ha='right',rotation_mode='anchor')
ax1.set_ylabel(r'$W/(\rm{J\cdot m}^{-3})$',fontsize=10)
ax1.tick_params(direction='out')
ax1.tick_params(axis='x',pad=2,labelsize=8)
ax1.tick_params(axis='y',pad=2,labelsize=8)
# ~ ax1.set_yscale('log')
ax1.ticklabel_format(style='sci',axis='y',scilimits=(-3,3))
ax1.xaxis.get_offset_text().set_size(8)
ax1.yaxis.get_offset_text().set_size(8)
plt.tight_layout(pad=0.1)
plt.savefig('W.pdf',transparent=True)
plt.savefig('W.png',dpi=600)
plt.close('all')
