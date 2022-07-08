"""
Created 08. July 2022 by Daniel Van Opdenbosch, Technical University of Munich

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed without any warranty or implied warranty of merchantability or fitness for a particular purpose. See the GNU general public license for more details: <http://www.gnu.org/licenses/>
"""

import numpy
import glob
import os
import ray
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy import signal
from scipy import stats

Punkte=100		#Punkte pro Segment zur Bestimmung von E
Dehngrenze=2e-4	#

#Groessen aus DIN 527
h=2e-3		#Probendicke
b1=5e-3		#Probenbreite
b2=10e-3	#Einspannbreite
L=58e-3		#Einspannlaenge
L0=25e-3	#Messlaenge
r=30e-3		#Radius der Fase
x=numpy.linspace(r,0,1000000)
alpha=b1/L*numpy.trapz(1/(b2-2*(r**2-(r-x)**2)**0.5),x)

def conv(t):
	return t.replace(',','.')

os.system('mv Results.log Results.alt')
sys.stdout=open('Results.log','a')

files=glob.glob('*[!_alt].txt',recursive=True)

@ray.remote
def mech(f):
	filename=os.path.splitext(f)[0].split('/')[-1]
	print(filename)

	Zeit_s,Kraft_N,Weg_mm,Spannung_MPa,Dehnung_perc=numpy.genfromtxt((conv(t) for t in open(f)),delimiter='\t',unpack=True,skip_header=1,skip_footer=0,usecols=range(5))
	Spannung=Spannung_MPa*1e6*(2*alpha+1)
	Dehnung=Dehnung_perc/1e2

	Spannung=Kraft_N/(h*b1)
	Dehnung=Weg_mm/(L0+Weg_mm[numpy.where(Kraft_N>0)][0])

	Spannung,Dehnung=Spannung[numpy.where(Spannung>0)],Dehnung[numpy.where(Spannung>0)]-Dehnung[numpy.where(Spannung>0)][0]

	R=max(Spannung)	#Punkte=int(numpy.where(Spannung==R)[0][0]/5)

	Agt0=float(Dehnung[numpy.where(Spannung==R)][0])
	steps=len(Dehnung[numpy.where(Dehnung<Agt0)])//Punkte
	theils=numpy.array([])
	for i in numpy.arange(steps):
		ind=numpy.where((Dehnung>=i*Agt0/steps)&(Dehnung<(1+i)*Agt0/steps))
		theils=numpy.append(theils,stats.theilslopes(Spannung[ind],Dehnung[ind]))
	theils=theils.reshape(steps,4)
	theil=theils[numpy.where(theils[:,0]==max(theils[:,0]))][0]
	E=theil[0]
	disp=theil[1]
	Econflo=theil[2]
	Econfup=theil[3]

	Dehnung+=disp/E
	Spannung,Dehnung=Spannung[numpy.where(Dehnung>0)],Dehnung[numpy.where(Dehnung>0)]

	Re=Spannung[numpy.where(Dehnung-Spannung/E>Dehngrenze)][0]

	indBruch=numpy.where(Spannung>=R/10)[-1][-1]
	Agt=float(Dehnung[numpy.where(Spannung==R)][0])
	At=float(Dehnung[indBruch])

	Ag=Agt-R/E
	A=At-float(numpy.median(Spannung[indBruch-5:indBruch]))/E

	Wt=numpy.trapz(Spannung[:indBruch],x=Dehnung[:indBruch])
	W=Wt-R**2/E/2

	print(filename,'R',R,'E',E,'A',A,'W',W,'Re',R,'Ag',Ag,'At',At,'Wt',Wt)

	plt.clf()
	mpl.rc('text',usetex=True)
	mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
	fig,ax1=plt.subplots(figsize=(7.5/2.54,5.3/2.54))

	ax1.plot(Dehnung,Dehnung*E,'k',linewidth=0.5,path_effects=[pe.Stroke(linewidth=1,foreground='white'),pe.Normal()])
	ax1.plot(Dehnung,Dehnung*Econflo,'k:',linewidth=0.5,path_effects=[pe.Stroke(linewidth=1,foreground='white'),pe.Normal()])
	ax1.plot(Dehnung,Dehnung*Econfup,'k:',linewidth=0.5,path_effects=[pe.Stroke(linewidth=1,foreground='white'),pe.Normal()])

	ax1.plot(Dehnung+Dehngrenze,Dehnung*E,'k--',linewidth=0.5,path_effects=[pe.Stroke(linewidth=1,foreground='white'),pe.Normal()])
	ax1.errorbar(Dehngrenze,0,marker='s',color='k',markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0,path_effects=[pe.Stroke(linewidth=2,foreground='w'),pe.Normal()],zorder=10)
	ax1.errorbar(Dehngrenze+Re/E,Re,marker='s',color='k',markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0,path_effects=[pe.Stroke(linewidth=2,foreground='w'),pe.Normal()],zorder=10)

	ax1.plot(Dehnung+Agt-R/E,Dehnung*E,'k--',linewidth=0.5,path_effects=[pe.Stroke(linewidth=1,foreground='white'),pe.Normal()])
	ax1.errorbar(Agt-R/E,0,marker='s',color='k',markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0,path_effects=[pe.Stroke(linewidth=2,foreground='w'),pe.Normal()],zorder=10)
	ax1.errorbar(Agt,R,marker='s',color='k',markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0,path_effects=[pe.Stroke(linewidth=2,foreground='w'),pe.Normal()],zorder=10)

	ax1.plot(Dehnung+At-float(numpy.median(Spannung[indBruch-5:indBruch]))/E,Dehnung*E,'k--',linewidth=0.5,path_effects=[pe.Stroke(linewidth=1,foreground='white'),pe.Normal()])
	ax1.errorbar(At-float(numpy.median(Spannung[indBruch-5:indBruch]))/E,0,marker='s',color='k',markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0,path_effects=[pe.Stroke(linewidth=2,foreground='w'),pe.Normal()],zorder=10)
	ax1.errorbar(At,float(numpy.median(Spannung[indBruch-5:indBruch])),marker='s',color='k',markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0,path_effects=[pe.Stroke(linewidth=2,foreground='w'),pe.Normal()],zorder=10)

	ax1.plot(Dehnung[numpy.where(Dehnung<=At*1.5)],Spannung[numpy.where(Dehnung<=At*1.5)],'k',linewidth=0.5,path_effects=[pe.Stroke(linewidth=1,foreground='white'),pe.Normal()])

	ax1.set_xlim([-max(Dehnung[numpy.where(Dehnung<=At*1.5)])*0.05,max(Dehnung[numpy.where(Dehnung<=At*1.5)])*1.05])
	ax1.set_ylim([-R*0.05,R*1.1])

	ax1.set_xlabel(r'$\epsilon/1$',fontsize=10)
	ax1.set_ylabel(r'$\sigma/\rm{Pa}$',fontsize=10)
	ax1.tick_params(direction='out')
	ax1.tick_params(axis='x',pad=2,labelsize=8)
	ax1.tick_params(axis='y',pad=2,labelsize=8)
	ax1.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
	ax1.xaxis.get_offset_text().set_size(8)
	ax1.yaxis.get_offset_text().set_size(8)
	plt.tight_layout(pad=0.1)
	plt.savefig(str(os.path.splitext(f)[0])+'.pdf',transparent=True)
	plt.savefig(str(os.path.splitext(f)[0])+'.png',dpi=600)
	plt.close('all')

	return filename,R,E,A,W,Re,Ag,At,Wt

data=ray.get([mech.remote(f) for f in files])

numpy.save('data.npy',data)
os.system('python3 read_Mech_R_E_A_W.py')
