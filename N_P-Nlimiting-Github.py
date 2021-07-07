'''
Created on March 9, 2021

final for the N-limiting

@author: Keisuke/garmin
'''

from pylab import *
from af001_energy_calculation import *
from matplotlib import pyplot
from matplotlib import *
import matplotlib.patches as mpat

global What_is_limiting

What_is_limiting=1                      #0: P-limiting  1:N-limiting

if What_is_limiting==0:
    #for P limiting data
    from Healey85_data_08_Chl_from_Nlimited_added import healey85C
elif What_is_limiting==1:
    #for N limiting data
    from Healey85_data_09_Nlimited_case import healey85C

#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#Function beging here
#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
def kkI():   

  
    I=70                                #irradiance (umol photons m-2 s-1)

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Parameters
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    if What_is_limiting==0:
    #For P-limiting case
        Pin=0.002                       #(mol/m3) Phosphorus concentration in the incoming medium (Healey 1985)
        Nin=0.2                         #(mol/m3) Nitrate concentration in the incoming medium (Healey 1985)
        Qc=1.00*10**(-12)/12            #(molC/cell) biomass C per cell (196-18)(from Healey 1985)
    
    elif What_is_limiting==1:
    #For N-limiting case
        Pin=0.02                        #(mol/m3) Phosphorus concentration in the incoming medium (Healey 1985)
        Nin=0.05                        #(mol/m3) Nitrate concentration in the incoming medium (Healey 1985)
        Qc=10**(-12)/12                 #(molC/cell) biomass C per cell (196-18)(from Healey 1985)

    E3=evalue()
    E=E3.E

    C=6

    Dd=0.25                             #growth rate 
    D=Dd/(3600*24)
    Mchl=893.49                         #(g / mol chlorophyll) mollar mass of chlorophyll
    
    #==============================
    #parameter sets
    #==============================
    m=3.791E-19                         #(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
    Pmax=0.003205
    OT=0.008634
    Ynphoto_chl=2.5                     #Chlamydomonas reinhardtii
    #Ynphoto_chl=3.561                  #((molN cell-1)/(molC chl cell-1)) the stoichiometric ratio for cell photosynthetic enzyme (Rubisco etc.) nitrogen to chlorophyll (193-25)
    Cnbiosynth=4.347E-10                #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
    Nconst_protein=4.453E-15            #(molN cell-1) Constant protein pool in nitrogen (193-25)
    Nstore_max=0                        #(molN cell-1) Constant protein pool in nitrogen (193-25)
    Cnrna_variable=6213*2.5             #(s) Constant for Variable part of RNA (193-26)
    Ypthylakoid_chl=0.02816             #((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
    Pconst_other=1.3*5.445E-17          #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
    Qp_max=25.26/(3.097e16)     
    Cessential=1.518E-15                #(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%
    Ea=70000                            #(J/mol) activation energy for phytoplankton (Li, 1980) cited in (Geider, 1997)
    R=8.3                               #(J/mol*k) universal gas constant 
    Tref=293                            #Reference temperature (K) from Healey experiment 20 celsius
    #------------------------------
    #Photosynthesis
    #------------------------------
    Tref=293 
    K = 273                                 #Reference temperature (K) 
    Tmax= 35 + K
    Tmin= K
    #Following steps to match 1600 array for the N:P output
    Tstep=3.8889
    Ttstep=3.88889
    Tt=arange(Tmin,Tmax+Ttstep,Ttstep)
    U = arange(size(Tt))
    
    
    A=Ea/R
    Arr=exp(-A*((1/Tt)-(1/Tref)))           #arrehenius equation (Geider, 1997) function of temperature 
    
    Tc = Tt - K

    Pchl=Arr*Pmax*(1-exp(-OT*I))            #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
    Pchl=Pchl/2                             #12:12 dark:light cycle leading to half photosynthesis
    
    Cnbiosynth = Cnbiosynth/Arr
    Cnrna_variable = Cnrna_variable/Arr
    
    ls=D*Qc                                 #(molC s-1) Biomass synthesis rate (193-25)
    #------------------------------
    
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Key parameters
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    
    Daverage=0.6*86400
    #================================
    #Constant parameters
    #================================

    Cnproteinsynth=24000                    #(s) Constant for protein synthesis nitrogen calculation (193-25)
 

    Cndna_variable=1.60e-15                 #(molN cell-1 s) Constant for variable part of DNA (193-26) 
        
    #================================
    
    #================================
    #Y intercept parameters
    #================================

    Nconst_dnarna=2.78e-15/5/2.8            #(molN cell-1) Constant part of nitrogen in DNA and RNA (193-26)
    Nconst_other=2.78e-16                   #(molN cell-1) Other kinds of nitrogen assuming constant (193-33)

    #================================
    
    Molar_mass_DNA_AT_average=307.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_DNA_CG_average=307.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    
    Molar_mass_RNA_AT_average=316.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_RNA_CG_average=323.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    
    #================================
    #E coli
    #================================
    CG_Ecoli=0.506                          #(dimensionless) from [http://www.ncbi.nlm.nih.gov/genome/167 (accessed 06/18/2016)]
    AT_Ecoli=1-CG_Ecoli                     #(dimensionless) 
    
    Molar_mass_DNA_Ecoli=Molar_mass_DNA_AT_average*CG_Ecoli+Molar_mass_DNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of DNA unit
    Molar_mass_RNA_Ecoli=Molar_mass_RNA_AT_average*CG_Ecoli+Molar_mass_RNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of RNA unit
    
    RNA_DNA_mass_ratio=20/7.6               #(ug/ug) Bremer and Dennis 1996
    RNA_DNA_mass_ratio=17.844/6.5239        #(ug/ug) from values ad D=0 "07 Bremer and Dennis 1996 data plot.xlsx"
    
    RNA_DNA_molar_ratio=RNA_DNA_mass_ratio/Molar_mass_RNA_Ecoli*Molar_mass_DNA_Ecoli               #(mol mol-1)
    #================================
    #Stoichiometric parameters
    #================================
    YcyanoC_N=2                                         #(molC molN) C/N molar ratio of cyanophycin
    YpgC_P=40                                           #(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)
    
    CG=0.563                                            #GC%    [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
    YnucacidP_N=1/(3.5*(1-CG)+4*CG)                     #(molP molN-1) P/N molar ratio of RNA (193-26) values (193-28) excel file "08 N to P ratio in DNA and RNA.xlsx"
    
    YdnaC_N=3.5*(1-CG)+2.5*CG                           #(molC molN-1) C/N molar ratio of dna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
    YrnaC_N=3.25*(1-CG)+2.5*CG                          #(molC molN-1) C/N molar ratio of rna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)

    DNAmb=2.1269                                        #(Mb) Megabase pair of synechococcus DNA in mega (million) base pairs [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
    Avogadro=6.022*10**23                               #(molecules mol-1) Avogadro constant
    Pdna_const=DNAmb*2*10**6/Avogadro                   #(molP cell-1) Constant part of DNA in phosphorus 
    Prna_const=Pdna_const*RNA_DNA_molar_ratio           #(molP cell-1) Constant part of RNA in phosphorus
    #* Make sure to multiply by 2 as they are base PAIRs"
    Ndna_const=Pdna_const/YnucacidP_N                   #(molN cell-1) Constant part of DNA in nitrogen
    Nrna_const=Ndna_const*RNA_DNA_molar_ratio           #(molN cell-1) Constatn part of RNA in phosphorus
    
    Ynrnaconst_dnarnaconst=1/2                          #(dimensionless) N molar ratio of RNA (constant part) to DNA + RNA (constant part) (193-33) reffering to around p.112 of Biology of Prokyariotes
    Yndnaconst_dnarnaconst=1-Ynrnaconst_dnarnaconst     #(dimensionless) N molar ratio of DNA (constant part) to DNA + RNA (constant part)  (193-33) refering to around p.112 of Biology of Prokyariotes
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Calculation
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    Chl=((1+E)*ls+m)/Pchl                               #(molC chl cell-1) Chlrophyll concentration (193-25) 
    Nphoto=Chl*Ynphoto_chl                              #(molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nproteinsynth=(Nphoto+Nconst_protein)*D*Cnproteinsynth/(1-D*Cnproteinsynth)   #(molN cell-1) protein synthesis related protein in N (193-25)
    Nbiosynth=D*Cnbiosynth                              #(molN cell-1) various part of biosynthesis related protein in N (193-37)
    Nprotein=Nphoto+Nconst_protein+Nbiosynth            #(molN cell-1) All the proteins in N (193-26)
    Nrna_variable=Nprotein*D*Cnrna_variable             #(molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
    Ndna_variable=Ndna_const*Dd/1.2*(18.3-7.6)/7.6      #(molN cell-1) variable part of nitrogen in DNA (193-26) Increasing ratio based on Bremmer 1996
    Ndna_variable=0*Dd                                  #(molN cell-1) While Bremer and Dennis shows increasing trend, Parrott 1980 shows decreasing trend. 
   
    Nchl=Chl*4/55                                       #(molN cell-1) Chlorophyll nitrogen (actually almost negligiable)
    Pthylakoid=Chl*Ypthylakoid_chl                      #(molP cel-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
    Prna_variable=Nrna_variable*YnucacidP_N             #(molP cell-1) variable part of phosphorus in RNA (193-26)
    Pdna_variable=Ndna_variable*YnucacidP_N             #(molP cell-1) variable part of phosphorus in DNA (193-26)
    
    #=================================
    #Total calculation
    #=================================
    Qn_max=Nprotein+Nrna_variable+Nrna_const+Ndna_variable+Ndna_const+Nchl+Nstore_max       #(molN cell-1)  nitrogen content in the cell (193-26)                                 #(molN cell-1) total phosphorus content in the cell (193-26)
                                                                                            #(molP cell-1) total phosphorus content in the cell (193-26)
    Qn_min=Nprotein+Nrna_variable+Nrna_const+Ndna_variable+Ndna_const+Nchl                  #(molN cell-1) total nitrogen in the cell without storage
    Qp_min=Pconst_other+Pthylakoid+Prna_variable+Prna_const+Pdna_variable+Pdna_const        #(molP cell-1) total phosphorus in the cell without storage
    
    #=================================
    #Vector preparation
    #=================================
    Nstore=zeros(size(Tt))
    X=zeros(size(Tt))
    Qn_test=zeros(size(Tt))
    Qp_test=copy(X)
    Qp=copy(X)
    Qn=copy(X)
    Pstore=copy(X)
    #=================================
    #Population calculation
    #=================================
    Xn_max=Nin/Qn_min
    Xp_max=Pin/Qp_min
    for i in U:
        if Xn_max[i]>Xp_max[i]:
            X[i]=Xp_max[i]
            Qp[i]=Qp_min[i]
            Qn_test[i]=Nin/X[i]
            if Qn_test[i]<Qn_max[i]:
                Qn[i]=Qn_test[i]
                Nstore[i]=Qn_test[i]-Nprotein[i]-Nrna_variable[i]-Nrna_const-Ndna_variable[i]-Ndna_const-Nchl[i]  #(molN cell-1) Nitrogen storage in the cell
            else:
                Qn[i]=Qn_max[i]
                Nstore[i]=Nstore_max
        else:
            X[i]=Xn_max[i]
            Qn[i]=Qn_min[i]
            Qp_test[i]=Pin/X[i]
            if Qp_test[i]<Qp_max:
                Qp[i]=Qp_test[i]
            else:
                Qp[i]=Qp_max
            Pstore[i]=Qp[i]-Pconst_other-Pthylakoid[i]-Prna_variable[i]-Prna_const-Pdna_variable-Pdna_const       #(molP cell-1) Stored phosphorus in the cell
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 1
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    BiomassC=12*X*Qc                            #(mg C L-1) Biomass concentration
    NtoCplot=Qn/Qc*14*10**6/(12*10**3)          #(ug N / mg C) biomass N to C ratio (164-20)
    PtoCplot=Qp/Qc*30.97*10**6/(12*10**3)       #(ug P / mg C) biomass P to C ratio (164-20)
    NtoPplot=Qn/Qp*14*10**6/(30.97*10**6)       #(ug N /ug P) biomass N to P ratio (164-20)
    ChltoC0=Chl/Qc                              #(mol C chl mol C -1) Chlorophyll to carbon ratio
    Mchl=893.49                                 #(g / mol chlorophyll) mollar mass of chlorophyll
    ChltoCplot=ChltoC0/12/1000*Mchl/55*10**6    #(ug chlorophyll a mg C-1) (see 157-36 for conversion)
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 2  (calculation of dna, rna, etc)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 3 (unit adjustment)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Nunit=1/Qc                                  #((ug N / mgC)/(molN cell-1) unit conversion term (164-20)
    Punit=1/Qc                                  #((ug P / mgC)/(molP cell-1) unit conversion term (164-20)
    Numbertoarray=ones(size(Tt))                #(dimensionless) Number to array converter
    
    #=======================================
    #Calculation of carbon usage (195-16)
    #=======================================
    CNprotein=4.49                              #(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
    
    #-------------------
    #C protein
    #-------------------
    Cphoto=Nphoto*CNprotein                     #(molC cell-1) carbon in photosystem protein (195-16)
    Cbiosynth=Nbiosynth*CNprotein               #(molC cell-1) carbon in biosynthesis protein (195-16)
    Cconst_protein=Nconst_protein*CNprotein     #(molC cell-1) carbon in other protein assumed constant (195-16)
 
    #----------------------
    #C chlorophyll
    #----------------------
    Cchl=Chl                                    #(molC cell-1) carbon in chlorophyll (195-16)
    
    #----------------------
    #C DNA RNA
    #----------------------
    Crna_const=Nrna_const*YrnaC_N               #(molC cell-1) carbon in variable part of RNA (195-16)
    Crna_variable=Nrna_variable*YrnaC_N         #(molC cell-1) carbon in variable part of RNA (195-16)
    
    Cdna_const=Ndna_const*YdnaC_N               #(molC cell-1) carbon in constant part of DNA (195-16)
    Cdna_variable=Ndna_variable*YdnaC_N         #(molC cell-1) carbon in variable part of DNA (195-16)
    

    Cnstore=Nstore*YcyanoC_N                    #(molC cell-1) carbon in nitrogen storage (cyanophycin)
    CthylakoidPG=Pthylakoid*YpgC_P              #(molC cell-1) carbon in PG (phosphatidyl glycerol) in thylakoid membranes
    
    #---------------------------------------------------
    #C other: Here revised to include Nstore reduction
    #---------------------------------------------------

    #print(Qc)
    Cother_without_Nstore=Qc-Cphoto-Cbiosynth-Cconst_protein-Cchl\
            -Crna_const-Crna_variable-Cdna_const-Cdna_variable\
            -Cessential-CthylakoidPG
    
    Cother_with_full_Nstore=Qc-Cphoto-Cbiosynth-Cconst_protein-Cchl\
            -Crna_const-Crna_variable-Cdna_const-Cdna_variable\
            -Cessential-Cnstore-CthylakoidPG
            
    Cother=Cother_with_full_Nstore            
    
    Nstore_reduce=logical_and(Cother_without_Nstore>0, Cother_with_full_Nstore<0)
    
    Cnstore[Nstore_reduce]=Cother_without_Nstore[Nstore_reduce]
    Nstore0=copy(Nstore)
    Nstore[Nstore_reduce]=Cother_without_Nstore[Nstore_reduce]/YcyanoC_N
    Cother[Nstore_reduce]=0
    Qn[Nstore_reduce]=Qn[Nstore_reduce]+Nstore[Nstore_reduce]-Nstore0[Nstore_reduce]
    
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 1
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    BiomassC=12*X*Qc                            #(mg C L-1) Biomass concentration
    NtoCplot=Qn/Qc*14*10**6/(12*10**3)          #(ug N / mg C) biomass N to C ratio (164-20)
    PtoCplot=Qp/Qc*30.97*10**6/(12*10**3)       #(ug P / mg C) biomass P to C ratio (164-20)
    NtoPplot=Qn/Qp                              #(ug N /ug P) biomass N to P ratio (164-20)
    ChltoC0=Chl/Qc                              #(mol C chl mol C -1) Chlorophyll to carbon ratio
    Mchl=893.49                                 #(g / mol chlorophyll) mollar mass of chlorophyll
    ChltoCplot=ChltoC0/12/1000*Mchl/55*10**6    #(ug chlorophyll a mg C-1) (see 157-36 for conversion)
    
    #=======================================
    #C for plot in
    #=======================================
    percentorratio=100                          #100: percent, 1:ratio
    Cphoto_plot=Cphoto/Qc*percentorratio           
    Cbiosynth_plot=Cbiosynth/Qc*percentorratio
    Cconst_protein_plot=Cconst_protein/Qc*percentorratio*Numbertoarray
    Cchl_plot=Cchl/Qc*percentorratio
    Crna_const_plot=Crna_const/Qc*percentorratio*Numbertoarray
    Crna_variable_plot=Crna_variable/Qc*percentorratio
    Cdna_const_plot=Cdna_const/Qc*percentorratio*Numbertoarray
    Cdna_variable_plot=Cdna_variable/Qc*percentorratio

    Cother_plot=Cother/Qc*percentorratio
    Cessential_plot=Cessential/Qc*percentorratio*Numbertoarray
    Cnstore_plot=Cnstore/Qc*percentorratio
    CthylakoidPG_plot=CthylakoidPG/Qc*percentorratio
    
    #=======================================
    #N for plot ***(Mainly from 193-33)***
    #=======================================
    Nphoto_plot=Nphoto*Nunit                                    #(ug N/ mgC) Photosynthesis related protein nitrogen (193-25)(193-33)
    Nbiosynth_plot=Nbiosynth*Nunit                              #(ug N/ mgC) biosynthesis related protein in N (193-37)
    Nconst_protein_plot=Nconst_protein*Nunit*Numbertoarray      #(ug N/ mgC) constant protein pool in nitrogen (193-25)(193-33)
    Nchl_plot=Nchl*Nunit                                        #(ug N/ mgC) Chlorophyll nitrogen (actually almost negligiable) (193-33)
    
    Nrna_variable_plot=Nrna_variable*Nunit                      #(ug N/ mgC) Nitrogen in Variable part of nucleic acid (193-37)
    Ndna_variable_plot=Ndna_variable*Nunit                      #(ug N/ mgC) Nitrogen in Variable part of nucleic acid (193-37)
    
    Ndna_const_plot=Ndna_const*Nunit*Numbertoarray              #(ug N/ mgC) Nitrogen in constant part of DNA
    Nrna_const_plot=Nrna_const*Nunit*Numbertoarray              #(ug N/ mgC) Nitrogen in constant part of RNA
    
    Nstore_plot=Nstore*Nunit                                    #(ug N/ mgC) Nitrogen in storage
    #=======================================
    #P for plot ***(Mainly from 193-33)***
    #=======================================
    Prna_variable_plot=Prna_variable*Punit                      #(ug P/mgC) Phosphorus in variable part of RNA (193-37)
    Pdna_variable_plot=Pdna_variable*Punit                      #(ug P/mgC) Phosphorus in variable part of DNA (193-37)
    
    Pthylakoid_plot=Pthylakoid*Punit                            #(ug P/ mgC) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)(193-33)
    Pconst_other_plot=Pconst_other*Punit*Numbertoarray          #(ug P/ mgC) PHosphorus in other parts (ex. phospholipid in outer membrane, ATP, ADP, etc. (assuming constant) (193-33)
    
    Pdna_const_plot=Pdna_const*Punit*Numbertoarray              #(ug P/ mgC) Phosphorus in constant part of DNA
    Prna_const_plot=Prna_const*Punit*Numbertoarray              #(ug P/ mgC) Phosphorus in constant part of RNA
    
    Pstore_plot=Pstore*Punit                                    #(ug P/ mgC) Phosphorus in phosphorus storage

    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Applying chl max
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Chlmax=0.048    
    
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#Plot free parameters
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    from decimal import Decimal
    def sci(Numb):
        Numb='%.2E' % Decimal(Numb)
        return Numb
    
    print("m",sci(m/Qc*86400))
    print("Pmax",sci(Pmax*86400))
    print("Apho",sci(OT))
    print("Ynphoto_chl",sci(Ynphoto_chl*CNprotein))
    print("Cother_protein",sci(Nconst_protein/Qc*CNprotein))
    print("Ythylakoid_chl_P",sci(Ypthylakoid_chl))
    print("Pconst_other",sci(Pconst_other/Qc))
    print("Nstore_max",sci(Nstore_max/Qc))
    print("Cessential",sci(Cessential/Qc))
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #5.Plot
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    rcParams.update({'font.size': 25,
                     'lines.markersize':12,
                     'lines.markeredgewidth':1})
    rcParams.update({'xtick.major.pad': 15})
    rcParams.update({'xtick.major.pad': 15})
    rcParams.update({'font.serif': 'Times New Roman'})
    rcParams.update({'figure.autolayout': True})
    rcParams['figure.figsize']=8,6.5
    rcParams.update({'figure.facecolor':'w'})
    rcParams.update({'lines.linewidth':3})   
    rcParams.update({'patch.edgecolor':'none'})
    
    rcParams.update({'axes.linewidth':1.5})
    rcParams.update({'xtick.major.width':1})
    rcParams.update({'ytick.major.width':1})
    rcParams.update({'mathtext.default': 'regular' })
    
     
    #==================================
    #Plot control * 1=on other=off
    #==================================
    Plot_C=1
    Plot_NC=1
    Plot_PC=1
    Plot_Chl=1
    Plot_NC_stack=1
    Plot_PC_stack=1
    Plot_C_stack=1
    Plot_C_stack_small=1
    #============================================
    #What is limiting: output folder controll
    #============================================
    if What_is_limiting==0:
        Whatislimiting="P-limiting"
    elif What_is_limiting==1:
        Whatislimiting="N-limiting"
    #============================================
         
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Stack plot part
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Color_DNA_const='orange'
    Color_RNA_const='yellow'
    Color_protein_const='blue'
    Color_photo='yellow'
    Color_RNA_variable='red'
    Color_DNA_variable='blue'
    Color_protein_biosynthesis='pink'
    Color_chl='#66FF66'
    Color_other='#006400'
    Color_other='#008000'
     
    Color_P_const='#CCFFFF'
    Color_DNA='black'
    Color_RNA='red'
    Color_Nstore='purple'
    Color_Pstore='#D0CECE'
    Color_Pstore='#88CCEE'  #cyan
    Color_Cessential='brown'
    Color_thylakoid='#FFD966'
    
    
    #needs to be changed colors for more color blind friendly figures
    Color_Photo='#CC6677'   #rose
    Color_Bio='#44AA99'     #teal
    Color_Other='#4B0082'   #indigo 
    Color_Nstore='#999933'  #olive
    Color_other='#DDCC77'   #sand
    
     
    StackPlotColorsN=(Color_Nstore,Color_protein_const,Color_photo,Color_protein_biosynthesis,Color_DNA,Color_RNA,Color_chl)
    StackPlotColorsP=(Color_P_const,Color_thylakoid,Color_DNA,Color_RNA,Color_Pstore)
    
    StackPlotColorsN=(Color_Other,Color_Photo,Color_Bio,Color_Nstore,Color_other)
    StackPlotColorsP=(Color_Other,Color_Photo,Color_Bio,Color_Pstore)
    StackPlotColorsC=(Color_Other,Color_Photo,Color_Bio,Color_Nstore,Color_other)
    StackPlotColorsCsmalls=(Color_DNA,Color_RNA,Color_chl,Color_thylakoid) 
     
    FigurenumberN_C=10

    FigurenumberC=20
    if Plot_C_stack==1:
        pyplot.figure(1,figsize=(7.7,6.5))
        pyplot.stackplot(Tc,Cessential_plot+Cconst_protein_plot+Cdna_const_plot+Cdna_variable_plot,\
                  Cphoto_plot+CthylakoidPG_plot+Cchl_plot,Crna_const_plot+Crna_variable_plot+Cbiosynth_plot,Cnstore_plot,Cother_plot,colors=StackPlotColorsC,alpha=0.75)
        print(Cbiosynth_plot) 
        pyplot.xlabel('Temperature (\u2103)', fontsize=25)                       
        pyplot.ylabel('C allocation ($\%$)', fontsize=25)  
        pyplot.xticks(arange(10,22,step=2),fontsize=20)
        pyplot.xlim(10,20)
        pyplot.ylim(ymax=percentorratio+1e-20)
        pyplot.yticks(fontsize=20)
        pyplot.title('Light: '+str(I)+' $\mu$mol m$^{-2}$ s$^{-1}$',y=1.02,fontsize=25)
        eleg=mpat.Patch(color='#4B0082', label="Other",alpha=0.75)
        pleg=mpat.Patch(color='#CC6677', label="Photosynthesis",alpha=0.75)
        bleg=mpat.Patch(color='#44AA99',label='Biosynthesis',alpha=0.75)
        sleg=mpat.Patch(color='#DDCC77',label='Carbon Storage',alpha=0.75)
        pyplot.legend(handles=[eleg,pleg,bleg,sleg],loc='upper center',bbox_to_anchor=(0.5,-0.25), ncol=2,fontsize='x-small',frameon=False)


        pyplot.figure(2,figsize=(7.7,6.5))
        pyplot.stackplot(Tc,Nconst_protein_plot+Ndna_const_plot+Ndna_variable_plot,Nphoto_plot+Nchl_plot,Nbiosynth_plot+Nrna_variable_plot+Nrna_const_plot,Nstore_plot,colors=StackPlotColorsN,alpha=0.75)
        pyplot.xlabel('Temperature (\u2103)', fontsize=25)                       
        pyplot.ylabel('N:C (mol/mol)', fontsize=25)  
        pyplot.xticks(arange(10,22,step=2),fontsize=20)
        pyplot.yticks(arange(0.00,0.35,step=0.05),fontsize=20)
        pyplot.xlim(10,20)
        pyplot.ylim(top=0.30)
        eleg=mpat.Patch(color='#4B0082', label="Other",alpha=0.75)
        pleg=mpat.Patch(color='#CC6677', label="Photosynthesis",alpha=0.75)
        bleg=mpat.Patch(color='#44AA99',label='Biosynthesis',alpha=0.75)
        psleg=mpat.Patch(color='#88CCEE',label='Phosphorus Storage',alpha=0.75)
        pyplot.legend(handles=[psleg,eleg],loc='upper center',bbox_to_anchor=(0.65,-0.25), ncol=3,fontsize='x-small',frameon=False)
        
        pyplot.figure(3,figsize=(7.7,6.5))
        pyplot.stackplot(Tc,Pconst_other_plot+Pdna_const_plot+Pdna_variable_plot,Pthylakoid_plot,Prna_variable_plot+Prna_const_plot,Pstore_plot,colors=StackPlotColorsP,edgecolor='none',alpha=0.75)
        pyplot.xlabel('Temperature (\u2103)',fontsize=25)                       
        pyplot.ylabel('P:C (mol/mol)',fontsize=25)  
        pyplot.xticks(arange(10,22,step=2),fontsize=20)
        pyplot.yticks(fontsize=20)
        pyplot.xlim(10,20)
        pyplot.ylim(top=0.01)
        eleg=mpat.Patch(color='#4B0082', label="Other",alpha=0.75)
        pleg=mpat.Patch(color='#CC6677', label="Photosynthesis",alpha=0.75)
        bleg=mpat.Patch(color='#44AA99',label='Biosynthesis',alpha=0.75)
        pyplot.legend(handles=[pleg,bleg],loc='upper center',bbox_to_anchor=(0.25,-0.25), ncol=2,fontsize='x-small',frameon=False)

  ######Entering the Thrane Data
    thrane=genfromtxt('thrane_d.csv',delimiter=',') 
      
    Thranet=thrane[:,0]
    ThraneNP=thrane[:,1]

 #####plotting   
    pyplot.figure(4)
    pyplot.plot(Tc,NtoPplot, color='#44AA99', linewidth=2,zorder=-1)
    pyplot.xlim(10,20)
    pyplot.ylim(bottom=0,top=30)
    pyplot.xlabel('Temperature (\u2103)', fontsize=25)
    pyplot.xticks(fontsize=20)
    pyplot.ylabel('N:P (mol/mol)', fontsize=25)
    pyplot.yticks(fontsize=20)
    modelleg=pyplot.plot(Tc,NtoPplot,'-',color='#44AA99', label="Model",zorder=-1)
    pyplot.legend(loc='lower right',fontsize='x-small',frameon=False)
    return

#AAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    
healey85data=healey85C()
kkI()

pyplot.show()
    
    