"""
Created 11. January 2023 by Daniel Van Opdenbosch, Technical University of Munich

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

def Bahadur(f,wDehnung,wSpannung):
	slope,intercept,low_slope,high_slope=stats.theilslopes(numpy.log(wSpannung),x=wDehnung)
	m,sigma0=slope,numpy.exp(intercept)
	plt.plot(wDehnung,sigma0*numpy.exp(m*wDehnung))
	plt.plot(wDehnung,wSpannung)
	plt.savefig(str(os.path.splitext(f)[0])+'_Bahadur.png')
	return m,sigma0

Dehngrenze=2e-4		#

#Groessen aus DIN 527
h=2e-3		#Probendicke
b1=5e-3		#Probenbreite
b2=10e-3	#Einspannbreite
L=58e-3		#Einspannlaenge
L0=25e-3	#Messlaenge
r=(b2-b1)/2	#Radius der Fase
x=numpy.linspace(0,r,1000000)
alpha=b1/L*numpy.trapz(1/(b2-2*(r**2-(r-x)**2)**0.5),x)

os.system('mv Results.log Results.alt')
sys.stdout=open('Results.log','a')

files=glob.glob('*[!_corr].txt',recursive=True)

@ray.remote
def mech(f,Dehngrenze,L,alpha,*args):
	filename=os.path.splitext(f)[0].split('/')[-1]
	print(filename)
	if 'Weg_G_mm' in str(open(f,'r').readlines()):
		Zeit_s,Kraft_N,Weg_mm,Spannung_MPa,Dehnung_perc,Weg_F_mm,Weg_G_mm=numpy.genfromtxt((t.replace(',','.') for t in open(f)),delimiter='\t',unpack=True,skip_header=1,skip_footer=0,usecols=range(7))
		if numpy.median(Weg_G_mm)!=0:
			Weg_mm,L,alpha=Weg_G_mm,L0,0
	else:
		Zeit_s,Kraft_N,Weg_mm,Spannung_MPa,Dehnung_perc=numpy.genfromtxt((t.replace(',','.') for t in open(f)),delimiter='\t',unpack=True,skip_header=1,skip_footer=0,usecols=range(5))
	Spannung=Kraft_N/(h*b1)
	Dehnung=Weg_mm/(L*1e3+(Weg_mm[numpy.where(Kraft_N>0)][0]-Weg_mm[0]))/(2*alpha+1)

	Spannung,Dehnung=Spannung[numpy.where(Spannung>0)],Dehnung[numpy.where(Spannung>0)]-Dehnung[numpy.where(Spannung>0)][0]

	R=max(Spannung)
	Punkte=int(numpy.where(Spannung==R)[0][0]/10)	#

	Agt0=float(Dehnung[numpy.where(Spannung==R)][0])
	steps=len(Dehnung[numpy.where(Dehnung<Agt0)])//Punkte
	theils=numpy.array([])
	for i in numpy.arange(steps):
		ind=numpy.where((Dehnung>=i*Agt0/steps)&(Dehnung<(1+i)*Agt0/steps))
		theils=numpy.append(theils,stats.theilslopes(Spannung[ind],x=Dehnung[ind]))
	theils=theils.reshape(steps,4)
	theil=theils[numpy.where(theils[:,0]==max(theils[:,0]))][0]
	E=theil[0]
	disp=theil[1]
	Econflo=theil[2]
	Econfup=theil[3]

	Dehnung+=disp/E
	Spannung,Dehnung=Spannung[numpy.where(Dehnung>0)],Dehnung[numpy.where(Dehnung>0)]

	indBruch=numpy.where(Spannung>=R/10)[-1][-1]
	Agt=float(Dehnung[numpy.where(Spannung==R)][0])
	At=float(Dehnung[indBruch])

	Ag=Agt-R/E
	A=At-float(numpy.median(Spannung[indBruch-5:indBruch]))/E

	if Ag<Dehngrenze:
		Dehngrenze=Ag
		Re=R
	else:
		Re=Spannung[numpy.where(Dehnung-Spannung/E>=Dehngrenze)][0]

	Wt=numpy.trapz(Spannung[:indBruch],x=Dehnung[:indBruch])
	W=Wt-R**2/E/2

	if 'Bahadur' in [*args]:
		Bereich=numpy.where(Dehnung-Spannung/E<Dehngrenze)
		m,sigma0=Bahadur(f,numpy.log(Dehnung[Bereich]+1),Spannung[Bereich]*(Dehnung[Bereich]+1))
		print(filename,'m',m,'sigma0',sigma0)

	print(filename,'R',R,'E',E,'A',A,'W',W,'Re',R,'Ag',Ag,'At',At,'Wt',Wt)
	numpy.savetxt(filename+'_corr.txt',numpy.transpose([Dehnung,Spannung]))

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
	ax1.tick_params(axis='both',pad=2,labelsize=8)
	ax1.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
	ax1.xaxis.get_offset_text().set_size(8)
	ax1.yaxis.get_offset_text().set_size(8)
	plt.tight_layout(pad=0.1)
	plt.savefig(str(os.path.splitext(f)[0])+'.pdf',transparent=True)
	plt.savefig(str(os.path.splitext(f)[0])+'.png',dpi=600)
	plt.close('all')

	return filename,R,E,A,W,Re,Ag,At,Wt

data=ray.get([mech.remote(f,Dehngrenze,L,alpha,sys.argv) for f in files])

numpy.save('data.npy',data)
os.system('python3 read_Mech_R_E_A_W.py')
