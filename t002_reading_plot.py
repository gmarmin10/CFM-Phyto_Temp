'''
Created on Aug 27, 2021

@author: garmin
'''

from pylab import *
from numpy import *
#from matplotlib import *
from matplotlib.pyplot import * #figure, xlabel, ylabel,colorbar,pcolormesh,show, xticks
#from final_verison2.py import *
from final_verison2 import cell


X = genfromtxt('X.csv',delimiter=',')
Y = genfromtxt('Y.csv',delimiter=',')
Z = genfromtxt('Z.csv',delimiter=',')

print(len(X))
print(len(Y))
print(len(Z))

NP=copy(Z)*nan
PC=copy(Z)*nan
NC=copy(Z)*nan

for i in arange(len(X)):
    for j in arange(len(Y)):
            T=Z[i,j]
            NP[i,j],PC[i,j],NC[i,j]=cell(T)

print(shape(NP))
print(shape(Z))

figure(1)
pcolormesh(X,Y,NP.T)
clim(0,50)
colorbar(label='(mol/mol)')
ylabel('Current Temperature')
#xlabel('Longitude')
title('N:P')
xticks((-150,-100,-50,0,50,100,150))
xticklabels = ["150\N{DEGREE SIGN}W","100\N{DEGREE SIGN}W","50\N{DEGREE SIGN}W","0","50\N{DEGREE SIGN}E","100\N{DEGREE SIGN}E","150\N{DEGREE SIGN}E"]
gca().set_xticklabels(xticklabels)
yticks((-60,-40,-20,0,20,40,60,80))
yticklabels = ["60\N{DEGREE SIGN}S","40\N{DEGREE SIGN}S","20\N{DEGREE SIGN}S","0","20\N{DEGREE SIGN}N","40\N{DEGREE SIGN}N","60\N{DEGREE SIGN}N","80\N{DEGREE SIGN}N"]
gca().set_yticklabels(yticklabels)

PC[isnan(Z)==1] = nan

figure(2)
pcolormesh(X,Y,PC.T)
clim(0,0.030)
colorbar(label='(mol/mol)')
#ylabel('Latitude')
#xlabel('Longitude')
title("P:C")
xticks((-150,-100,-50,0,50,100,150))
xticklabels = ["150\N{DEGREE SIGN}W","100\N{DEGREE SIGN}W","50\N{DEGREE SIGN}W","0","50\N{DEGREE SIGN}E","100\N{DEGREE SIGN}E","150\N{DEGREE SIGN}E"]
gca().set_xticklabels(xticklabels)
yticks((-60,-40,-20,0,20,40,60,80))
yticklabels = ["60\N{DEGREE SIGN}S","40\N{DEGREE SIGN}S","20\N{DEGREE SIGN}S","0","20\N{DEGREE SIGN}N","40\N{DEGREE SIGN}N","60\N{DEGREE SIGN}N","80\N{DEGREE SIGN}N"]
gca().set_yticklabels(yticklabels)

figure(3)
pcolormesh(X,Y,NC.T)
clim(0,0.35)
colorbar(label='(mol/mol)')
#ylabel('Latitude')
#xlabel('Longitude')
title('N:C')
xticks((-150,-100,-50,0,50,100,150))
xticklabels = ["150\N{DEGREE SIGN}W","100\N{DEGREE SIGN}W","50\N{DEGREE SIGN}W","0","50\N{DEGREE SIGN}E","100\N{DEGREE SIGN}E","150\N{DEGREE SIGN}E"]
gca().set_xticklabels(xticklabels)
yticks((-60,-40,-20,0,20,40,60,80))
yticklabels = ["60\N{DEGREE SIGN}S","40\N{DEGREE SIGN}S","20\N{DEGREE SIGN}S","0","20\N{DEGREE SIGN}N","40\N{DEGREE SIGN}N","60\N{DEGREE SIGN}N","80\N{DEGREE SIGN}N"]
gca().set_yticklabels(yticklabels)

NP2=copy(Z)*nan
PC2=copy(Z)*nan
NC2=copy(Z)*nan

for i in arange(len(X)):
    for j in arange(len(Y)):
            T=Z[i,j]+4
            NP2[i,j],PC2[i,j],NC2[i,j]=cell(T)
figure(4)
pcolormesh(X,Y,NP2.T)
clim(0,50)
colorbar(label='(mol/mol)')
ylabel('4 \u2103 increase')
#xlabel('Longitude')
#title('N:P')
xticks((-150,-100,-50,0,50,100,150))
xticklabels = ["150\N{DEGREE SIGN}W","100\N{DEGREE SIGN}W","50\N{DEGREE SIGN}W","0","50\N{DEGREE SIGN}E","100\N{DEGREE SIGN}E","150\N{DEGREE SIGN}E"]
gca().set_xticklabels(xticklabels)
yticks((-60,-40,-20,0,20,40,60,80))
yticklabels = ["60\N{DEGREE SIGN}S","40\N{DEGREE SIGN}S","20\N{DEGREE SIGN}S","0","20\N{DEGREE SIGN}N","40\N{DEGREE SIGN}N","60\N{DEGREE SIGN}N","80\N{DEGREE SIGN}N"]
gca().set_yticklabels(yticklabels)

PC2[isnan(Z)==1] = nan

figure(5)
pcolormesh(X,Y,PC2.T)
clim(0,0.030)
colorbar(label='(mol/mol)')
#ylabel('Latitude')
#xlabel('Longitude')
#title('P:C')
xticks((-150,-100,-50,0,50,100,150))
xticklabels = ["150\N{DEGREE SIGN}W","100\N{DEGREE SIGN}W","50\N{DEGREE SIGN}W","0","50\N{DEGREE SIGN}E","100\N{DEGREE SIGN}E","150\N{DEGREE SIGN}E"]
gca().set_xticklabels(xticklabels)
yticks((-60,-40,-20,0,20,40,60,80))
yticklabels = ["60\N{DEGREE SIGN}S","40\N{DEGREE SIGN}S","20\N{DEGREE SIGN}S","0","20\N{DEGREE SIGN}N","40\N{DEGREE SIGN}N","60\N{DEGREE SIGN}N","80\N{DEGREE SIGN}N"]
gca().set_yticklabels(yticklabels)

figure(6)
pcolormesh(X,Y,NC2.T)
clim(0,0.35)
colorbar(label='(mol/mol)')
#ylabel('Latitude')
#xlabel('Longitude')
#title('N:C')
xticks((-150,-100,-50,0,50,100,150))
xticklabels = ["150\N{DEGREE SIGN}W","100\N{DEGREE SIGN}W","50\N{DEGREE SIGN}W","0","50\N{DEGREE SIGN}E","100\N{DEGREE SIGN}E","150\N{DEGREE SIGN}E"]
gca().set_xticklabels(xticklabels)
yticks((-60,-40,-20,0,20,40,60,80))
yticklabels = ["60\N{DEGREE SIGN}S","40\N{DEGREE SIGN}S","20\N{DEGREE SIGN}S","0","20\N{DEGREE SIGN}N","40\N{DEGREE SIGN}N","60\N{DEGREE SIGN}N","80\N{DEGREE SIGN}N"]
gca().set_yticklabels(yticklabels)



NP3=copy(Z)*nan
PC3=copy(Z)*nan
NC3=copy(Z)*nan

for i in arange(len(X)):
    for j in arange(len(Y)):
            NP3[i,j]=NP2[i,j]-NP[i,j]
            PC3[i,j]=PC2[i,j]-PC[i,j]
            NC3[i,j]=NC2[i,j]-NC[i,j]


figure(7)
pcolormesh(X,Y,NP3.T)
clim(0,6)
colorbar(label='(mol/mol)')
ylabel('Change')
#xlabel('Longitude')
#title('N:P')
xticks((-150,-100,-50,0,50,100,150))
xticklabels = ["150\N{DEGREE SIGN}W","100\N{DEGREE SIGN}W","50\N{DEGREE SIGN}W","0","50\N{DEGREE SIGN}E","100\N{DEGREE SIGN}E","150\N{DEGREE SIGN}E"]
gca().set_xticklabels(xticklabels)
yticks((-60,-40,-20,0,20,40,60,80))
yticklabels = ["60\N{DEGREE SIGN}S","40\N{DEGREE SIGN}S","20\N{DEGREE SIGN}S","0","20\N{DEGREE SIGN}N","40\N{DEGREE SIGN}N","60\N{DEGREE SIGN}N","80\N{DEGREE SIGN}N"]
gca().set_yticklabels(yticklabels)

PC3[isnan(Z)==1] = nan

figure(8)
pcolormesh(X,Y,PC3.T)
set_cmap('viridis_r')
clim(-0.010,0)
colorbar(label='(mol/mol)')                      #ticks=range(),label='digit values'
#ylabel('Latitude')
#xlabel('Longitude')
#title('P:C')
xticks((-150,-100,-50,0,50,100,150))
xticklabels = ["150\N{DEGREE SIGN}W","100\N{DEGREE SIGN}W","50\N{DEGREE SIGN}W","0","50\N{DEGREE SIGN}E","100\N{DEGREE SIGN}E","150\N{DEGREE SIGN}E"]
gca().set_xticklabels(xticklabels)
yticks((-60,-40,-20,0,20,40,60,80))
yticklabels = ["60\N{DEGREE SIGN}S","40\N{DEGREE SIGN}S","20\N{DEGREE SIGN}S","0","20\N{DEGREE SIGN}N","40\N{DEGREE SIGN}N","60\N{DEGREE SIGN}N","80\N{DEGREE SIGN}N"]
gca().set_yticklabels(yticklabels)

figure(9)
pcolormesh(X,Y,NC3.T)
clim(-0.12,0)
set_cmap('viridis_r')
colorbar(label='(mol/mol)')
#ylabel('Latitude')
#xlabel('Longitude')
#title('N:C')
xticks((-150,-100,-50,0,50,100,150))
xticklabels = ["150\N{DEGREE SIGN}W","100\N{DEGREE SIGN}W","50\N{DEGREE SIGN}W","0","50\N{DEGREE SIGN}E","100\N{DEGREE SIGN}E","150\N{DEGREE SIGN}E"]
gca().set_xticklabels(xticklabels)
yticks((-60,-40,-20,0,20,40,60,80))
yticklabels = ["60\N{DEGREE SIGN}S","40\N{DEGREE SIGN}S","20\N{DEGREE SIGN}S","0","20\N{DEGREE SIGN}N","40\N{DEGREE SIGN}N","60\N{DEGREE SIGN}N","80\N{DEGREE SIGN}N"]
gca().set_yticklabels(yticklabels)

show()
