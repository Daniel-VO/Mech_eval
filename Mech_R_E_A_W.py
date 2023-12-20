"""
Created 20. December 2023 by Daniel Van Opdenbosch, Technical University of Munich

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed without any warranty or implied warranty of merchantability or fitness for a particular purpose. See the GNU general public license for more details: <http://www.gnu.org/licenses/>
"""

import os
import ray
import sys
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy import stats

def Bahadur(f,wDehnung,wSpannung):
	slope,intercept,low_slope,high_slope=stats.theilslopes(np.log(wSpannung),x=wDehnung)
	m,sigma0=slope,np.exp(intercept)
	plt.plot(wDehnung,sigma0*np.exp(m*wDehnung))
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
x=np.linspace(0,r,1000000)
alpha=b1/L*np.trapz(1/(b2-2*(r**2-(r-x)**2)**0.5),x)

os.system('mv results.log results.alt')

files=glob.glob('**/*[!_corr].txt',recursive=True)

@ray.remote
def mech(f,Dehngrenze,L,alpha,*args):
	filename=os.path.splitext(f)[0].split('/')[-1]
	if 'Weg_F_mm' in str(open(f,'r').readlines()):
		Zeit_s,Kraft_N,Weg_mm,Spannung_MPa,Dehnung_perc,Weg_F_mm,Weg_G_mm=np.genfromtxt((t.replace(',','.') for t in open(f)),delimiter='\t',unpack=True,skip_header=1,skip_footer=0,usecols=range(7))
		if np.median(Weg_F_mm)!=0 and max(Weg_G_mm-Weg_F_mm)<1:
			Weg_mm,L,alpha=Weg_F_mm,L0,0
		elif np.median(Weg_G_mm)!=0:
			Weg_mm,L,alpha=Weg_G_mm,L0,0
	else:
		Zeit_s,Kraft_N,Weg_mm,Spannung_MPa,Dehnung_perc=np.genfromtxt((t.replace(',','.') for t in open(f)),delimiter='\t',unpack=True,skip_header=1,skip_footer=0,usecols=range(5))
	Spannung=Kraft_N/(h*b1)
	Dehnung=Weg_mm/(L*1e3+(Weg_mm[np.where(Kraft_N>0)][0]-Weg_mm[0]))/(2*alpha+1)

	Spannung,Dehnung=Spannung[np.where(Spannung>0)],Dehnung[np.where(Spannung>0)]-Dehnung[np.where(Spannung>0)][0]

	R=max(Spannung)

	Schritte=5	#
	Punkte=int(np.where(Spannung==R)[0][0]/Schritte)
	theils=np.array([])
	for i in np.arange(1,Schritte):
		ind=np.arange(i*Punkte,(i+1)*Punkte)
		theils=np.append(theils,stats.theilslopes(Spannung[ind],x=Dehnung[ind]))
	theils=theils.reshape(Schritte-1,4)
	try:
		theil=theils[np.where(theils[:,0]==max(theils[:,0]))][0]
	except:
		os.system('mv '+f+' '+filename+'.aus')
	E,disp,Econflo,Econfup=theil[0],theil[1],theil[2],theil[3]

	Dehnung+=disp/E
	Spannung,Dehnung=Spannung[np.where(Dehnung>0)],Dehnung[np.where(Dehnung>0)]

	indBruch=np.where(Spannung>=R/10)[-1][-1]
	Agt=float(Dehnung[np.where(Spannung==R)][0])
	At=float(Dehnung[indBruch])

	Ag=Agt-R/E
	A=At-float(np.median(Spannung[indBruch-5:indBruch]))/E

	if Ag<Dehngrenze:
		Dehngrenze=Ag
		Re=R
	else:
		Re=Spannung[np.where(Dehnung-Spannung/E>=Dehngrenze)][0]

	Wt=np.trapz(Spannung[:indBruch],x=Dehnung[:indBruch])
	W=Wt-R**2/E/2

	if 'Bahadur' in [*args]:
		Bereich=np.where(Dehnung-Spannung/E<Dehngrenze)
		m,sigma0=Bahadur(f,np.log(Dehnung[Bereich]+1),Spannung[Bereich]*(Dehnung[Bereich]+1))
		print(filename,'m',m,'sigma0',sigma0)

	with open('results.log','a') as logfile:
		logfile.write(str([filename,'R:',R,' E:',E,' A:',A,' W:',W,' Re:',Re,' Ag:',Ag,' At:',At,' Wt:',Wt]).replace('[','').replace(']','').replace("'","")+'\n')
	np.savetxt(filename+'_corr.txt',np.transpose([Dehnung,Spannung]))

	plt.close('all')
	mpl.rc('text',usetex=True)
	mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
	plt.figure(figsize=(7.5/2.54,5.3/2.54))

	plt.plot(Dehnung,Dehnung*E,'k',linewidth=0.5,path_effects=[pe.Stroke(linewidth=1,foreground='white'),pe.Normal()])
	plt.plot(Dehnung,Dehnung*Econflo,'k:',linewidth=0.5,path_effects=[pe.Stroke(linewidth=1,foreground='white'),pe.Normal()])
	plt.plot(Dehnung,Dehnung*Econfup,'k:',linewidth=0.5,path_effects=[pe.Stroke(linewidth=1,foreground='white'),pe.Normal()])

	plt.plot(Dehnung+Dehngrenze,Dehnung*E,'k--',linewidth=0.5,path_effects=[pe.Stroke(linewidth=1,foreground='white'),pe.Normal()])
	plt.errorbar(Dehngrenze,0,marker='s',color='k',markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0,path_effects=[pe.Stroke(linewidth=2,foreground='w'),pe.Normal()],zorder=10)
	plt.errorbar(Dehngrenze+Re/E,Re,marker='s',color='k',markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0,path_effects=[pe.Stroke(linewidth=2,foreground='w'),pe.Normal()],zorder=10)

	plt.plot(Dehnung+Agt-R/E,Dehnung*E,'k--',linewidth=0.5,path_effects=[pe.Stroke(linewidth=1,foreground='white'),pe.Normal()])
	plt.errorbar(Agt-R/E,0,marker='s',color='k',markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0,path_effects=[pe.Stroke(linewidth=2,foreground='w'),pe.Normal()],zorder=10)
	plt.errorbar(Agt,R,marker='s',color='k',markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0,path_effects=[pe.Stroke(linewidth=2,foreground='w'),pe.Normal()],zorder=10)

	plt.plot(Dehnung+At-float(np.median(Spannung[indBruch-5:indBruch]))/E,Dehnung*E,'k--',linewidth=0.5,path_effects=[pe.Stroke(linewidth=1,foreground='white'),pe.Normal()])
	plt.errorbar(At-float(np.median(Spannung[indBruch-5:indBruch]))/E,0,marker='s',color='k',markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0,path_effects=[pe.Stroke(linewidth=2,foreground='w'),pe.Normal()],zorder=10)
	plt.errorbar(At,float(np.median(Spannung[indBruch-5:indBruch])),marker='s',color='k',markersize=1,elinewidth=0.5,capthick=0.5,capsize=2,linewidth=0,path_effects=[pe.Stroke(linewidth=2,foreground='w'),pe.Normal()],zorder=10)

	plt.plot(Dehnung[np.where(Dehnung<=At*1.5)],Spannung[np.where(Dehnung<=At*1.5)],'k',linewidth=0.5,path_effects=[pe.Stroke(linewidth=1,foreground='white'),pe.Normal()])

	plt.xlim([-max(Dehnung[np.where(Dehnung<=At*1.5)])*0.05,max(Dehnung[np.where(Dehnung<=At*1.5)])*1.05])
	plt.ylim([-R*0.05,R*1.1])

	plt.xlabel(r'$\epsilon/1$',fontsize=10)
	plt.ylabel(r'$\sigma/\rm{Pa}$',fontsize=10)
	plt.tick_params(axis='both',pad=2,labelsize=8)
	plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
	plt.tight_layout(pad=0.1)
	plt.savefig(str(os.path.splitext(f)[0])+'.pdf')
	plt.savefig(str(os.path.splitext(f)[0])+'.png',dpi=300)

	return filename,R,E,A,W,Re,Ag,At,Wt

np.save('data.npy',ray.get([mech.remote(f,Dehngrenze,L,alpha,sys.argv) for f in files]))
os.system('python3 Mech_R_E_A_W_read.py')
