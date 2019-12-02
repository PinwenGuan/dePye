import os
import numpy as np
from readst import *
from sympy import *
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from pymatgen.analysis.eos import *
import pymatgen as mg
import time

"""Read experimental data in expt file"""

def read_expt(expt='expt'):
    alias=''
    with open(expt,'r') as fexpt:
        cexpt0=fexpt.readlines()
        for i in cexpt0:
            if '#alias:' in i:
                alias=i.strip('#alias:').strip('\n').split(',')            
        cexpt=[i.replace('\n',' ').replace('\t',' ').strip(' ') for i in cexpt0 if not '#' in i]
    for i in range(cexpt.count('')):
        cexpt.remove('')
    idata=[]
    for i in range(len(cexpt)):
        cexpt[i]=cexpt[i].strip(' ').split(' ')
        for j in range(cexpt[i].count('')):
            cexpt[i].remove('')
        if len(cexpt[i])==2:
            cexpt[i].append('')
        if len(cexpt[i])>3:
            for j in range(3,len(cexpt[i])):
                cexpt[i][2]=cexpt[i][2]+' '+cexpt[i][j]
        if cexpt[i][2]=='' or cexpt[i][2]=='?':
            cexpt[i][2]='0'
        if cexpt[i][0]=='T':
            idata.append(i)
    idata.append(len(cexpt))

    data=[]
    for i in range(len(idata)-1):
        ref=[]
        for j in range(idata[i]+1,idata[i+1]):
            ref.append(cexpt[j][2])
        reff=[]
        [reff.append(k) for k in ref if not k in reff]
        for j in reff:
            t=[];v=[]
            for k in range(idata[i]+1,idata[i+1]):
                if cexpt[k][2]==j:
                    t.append(float(cexpt[k][0]))
                    v.append(float(cexpt[k][1]))
            data.append({'property':cexpt[idata[i]][1], 'T':t, 'value':v, 'ref':str(j)})    

    ref_num=[int(data[i]['ref']) for i in range(len(data))]    
    if alias=='':
        alias=[str(i+1) for i in range(max(ref_num))]

    return data,alias

"""Obtain highest common factor"""

def hcf(a):
    amin=min(a)
    for i in range(1,amin+1):
        s=1
        for j in a:
            if j%i==0:
                s=s
            else:
                s=0
        if s==1:
            hcf=i
    return hcf  

"""trim E-V data"""

def trim_ve(v,e):
    emin=min(e)
    for i in range(len(e)):
        if e[i]==emin:
            imin=i
    v=v[max(0,imin-3):min(len(e),imin+4)]
    e=e[max(0,imin-3):min(len(e),imin+4)]
    return v,e

"""EOS fitting based on pymatgen. Return equilibrium volume, energy and bulk modulus in angstrom-eV unit"""
"""     "murnaghan": Murnaghan,
        "birch": Birch,
        "birch_murnaghan": BirchMurnaghan,
        "pourier_tarantola": PourierTarantola,
        "vinet": Vinet,
        "deltafactor": DeltaFactor,
        "numerical_eos": NumericalEOS         """
#VEB0_PYMATGEN
def VEB0(v,e,eos='vinet'):
    Eos = EOS(eos_name=eos)
    eosfit=Eos.fit(v,e)
    V0=eosfit.v0
    E0=eosfit.e0
    B0=eosfit.b0
    BP0=eosfit.b1
    vfit=[(1.3-i/100)*v[0]+(-0.3+i/100)*v[-1] for i in range(161)]
    efit=eosfit.func(vfit)
    return [V0,E0,B0,BP0],[vfit,efit]

"""get volume from pressure based on the vinet equation"""

def p2v(p,b0,bp0):
    from scipy.optimize import fsolve
    def pv(x):
        y=3*b0*(1-x)/x**2*np.exp(1.5*(bp0-1)*(1-x))-p
        return y
    v=fsolve(pv,0.5)
    w=[]
    for i in range(len(v)):
        if v[i]>0:
            w.append(v[i])
    return w[0]

def d(z):
    f=(z)**3/(np.exp(z)-1)
    return f

"""Use Simpson's rule for integral in the Debye function, since the simpler method does not work"""
"""n is the number of intervals in 1"""

def D(x,n=100):
    N=math.floor(n*x)
    f=[0 for i in range(N)]
    f[0]=(4*d(0.5*x/N)+d(x/N))/6
    for i in range(1,min(N,int(700*N/x-1))):
        f[i]=(d(i*x/N)+4*d((i+0.5)*x/N)+d((i+1)*x/N))/6
    s=sum(f)*x/N*3/x**3    
    return s

def d2(z):
    f=(z)**4*np.exp(z)/(np.exp(z)-1)**2
    return f

def D2(x,n=100):
    N=math.floor(n*x)
    f=[0 for i in range(N)]
    f[0]=(4*d(0.5*x/N)+d(x/N))/6
    for i in range(1,min(N,int(700*N/x-1))):
        f[i]=(d2(i*x/N)+4*d2((i+0.5)*x/N)+d2((i+1)*x/N))/6
    s=sum(f)*x/N*3/x**3
    return s

"""calculation based on the Debye model from a single set of E-V data"""
"""N is number of atoms corresponding to V and E, nat is number of atoms in output. Only F is normalized in per atom""" 

def debye_single(V,E,N,M,s=0.617,Tmax=1000,nat=1,eos_only='False',eos='vinet'):
    hbar=1.054571800e-34;kB=1.38064852e-23;q=1.60217662e-19;NA=6.022140e23
    A=(6*np.pi**2)**(1/3)*hbar/kB

    veby=VEB0(V,E,eos=eos)
    [V0,E0,B0,BP0]=veby[0]
    B0=B0*q*1e30    # to Pa
    EVfit=veby[1]

    y=(1+BP0)/2-2/3        # high-temperature

    T_debye0=s*A*(V0*1e-30/N)**(1/6)*(B0/M)**(1/2)
    T_debye=[0 for i in range(len(V))]
    for i in range(len(V)):
        T_debye[i]=T_debye0*(V0/V[i])**y

    if eos_only=='True':
        return [V0,E0,B0/1e9,BP0,y,T_debye0],EVfit    

    start = time.time()
    T=list(range(0,Tmax+3,10))
    F=[[(E[i]-E0)/N+9/8*kB/q*T_debye[i] for i in range(len(V))] for j in range(len(T))]
    for i in range(len(V)):
        for j in range(1,len(T)):
            F[j][i]=(E[i]-E0)/N+9/8*kB/q*T_debye[i]-kB/q*T[j]*(D(T_debye[i]/T[j])-3*math.log(1-np.exp(-T_debye[i]/T[j])))

    Cv=[0 for j in range(len(T))]
    for j in range(1,len(T)):
        Cv[j]=3*nat*kB*NA*D2(T_debye0/T[j])        # J/mol/K
    end = time.time()
    print('Time for F(V,T) and Cv(T):'+str(end-start))

    start = time.time()
    G=[0 for j in range(0,len(T))]
    VT=[0 for j in range(0,len(T))]
    BT=[0 for j in range(0,len(T))]
    BP=[0 for j in range(0,len(T))]
    v=[i/N for i in V]
    for j in range(0,len(T)):
        veb0=VEB0(v,F[j],eos=eos)[0]
        G[j]=veb0[1]*q*NA/1e3*nat        # kJ/mol 
        VT[j]=veb0[0]*1e-24*NA*nat         # cm^3/mol 
        BT[j]=veb0[2]*q*1e30/1e9        # GPa
        BP[j]=veb0[3]
    
    S=[0 for j in range(0,len(T)-1)]
    H=[0 for j in range(0,len(T)-1)];H[0]=G[0]
    TEC=[0 for j in range(0,len(T)-1)]
    for j in range(1,len(T)-1):
        S[j]=-(G[j+1]*1e3-G[j-1]*1e3)/(T[j+1]-T[j-1])        # J/mol/K
        H[j]=G[j]+T[j]*S[j]/1e3                    # kJ/mol 
        TEC[j]=(VT[j+1]-VT[j-1])/(T[j+1]-T[j-1])/VT[j]*1e6    # 1e-6/K

    Cp=[0 for j in range(0,len(T)-2)]
    for j in range(1,len(T)-2):
        Cp[j]=(H[j+1]*1e3-H[j-1]*1e3)/(T[j+1]-T[j-1])    # J/mol/K

    end = time.time()
    print('Time for the other quantities:'+str(end-start))

    return [V0,E0,B0/1e9,BP0,y,T_debye0],EVfit,[T,G,VT,BT,S,H,TEC,Cp,Cv,BP],F

"""calculate, plot and save thermodynamic properties based on the Debye model"""

def debye(ev='e-v.dat',trimve='F',expt='expt',struc='POSCAR',out='Debye-thermo',s=0.617,poisson='poisson',Tmax=1000,Tstep_write=10,show='T',write='T',nat=0,prop='all',eos_only='False',eos='vinet'):
    hbar=1.054571800e-34;kB=1.38064852e-23;q=1.60217662e-19;NA=6.022140e23
    A=(6*np.pi**2)**(1/3)*hbar/kB

    print('\n')
    print('                  _            ______                                    ')
    print('                 | |          |  ____ \                   1.0            ')
    print('                 | |          | |    \ \                                 ')
    print('             ____| |   ____   | |     ) |_       _    ____               ')
    print('           / ____  | / ____ \ | |____/ /| |     | | / ____ \             ')
    print('          / /    \ |/ /____\ \|  _____/ | |     | |/ /____\ \            ')
    print('         | (     | |  ________| |       | |     | |  ________)           ')
    print('          \ \____/ |\ \______ | |        \ \___/ / \ \______             ')
    print('           \_______| \______/ |_|         \___  /   \______/             ')
    print('                                             / /                         ')
    print('                                          __/ /                          ')
    print('                                         |___/                           ')
    print('\n')
    print('       ______________________________________________________________________________________________________________           ')
    print('      |                                                                                                              |          ')
    print('      | The Debye model was first proposed in 1912 by Peter Joseph William Debye (March 24, 1884 - November 2, 1966) |          ')
    print('      |                                                                                                              |          ')
    print('       --------------------------------------------------------------------------------------------------------------        ')
    print('\n')
    print('Please use the Debye model with caution: it may fail when there is large variation of stength of different bonds in the solid')
    print('\n')
    print('                                                        - Pinwen Guan')
    print('\n')

    print('Use '+eos+' equation')
    if os.path.exists(poisson):
        with open(poisson,'r') as fpoisson:
            cpoisson=[i.strip('\n').split(' ') for i in fpoisson.readlines() if i!='\n' and i!=' \n' and not '#' in i]
        for j in range(cpoisson[0].count('')):
            cpoisson[0].remove('')
        cpoisson[0]=[float(k) for k in cpoisson[0]]
        poisson_ratio=cpoisson[0][0]
        s=3**(5/6)*(4*2**0.5*((1+poisson_ratio)/(1-2*poisson_ratio))**1.5+((1+poisson_ratio)/(1-poisson_ratio))**1.5)**(-1/3)
        print('Scaling factor from Poisson ratio: '+str(s))

    with open(ev,'r') as fev:
        cev=[i.strip('\n').replace('\t',' ').split(' ') for i in fev.readlines() if i!='\n' and i!=' \n' and not '#' in i]
    cname=[i for i in cev if 'v' in i or 'V' in i]
    cev=[i for i in cev if not 'v' in i and not 'V' in i]
    for i in range(len(cev)):
        for j in range(cev[i].count('')):
            cev[i].remove('')
        cev[i]=[float(k) for k in cev[i]] # A, eV

    if len(cname) != 0:
        cna=cname[-1]
        for j in range(cna.count('')):
            cna.remove('')
        cna=[str(k) for k in cna]
    else:
        cna=['' for i in range(len(cev[-1]))]        

    V=[i[0] for i in cev]
    E=[[] for i in range(len(cev[0])-1)]
    for i in range(len(cev[0])-1):
        E[i]=[j[i+1] for j in cev]
    pp=readst(struc)[2]
    an=[pp[i][0] for i in range(len(pp))]
    dict = {}
    for key in an:
        dict[key] = dict.get(key, 0) + 1

    ann = []
    [ann.append(i) for i in an if not i in ann]
    st=[dict[i] for i in ann]
    N=sum(st)
    if nat==0:
        sthcf=hcf(st)
        nat=N/sthcf
    M=1
    for i in range(len(ann)):
        amass=mg.Element(ann[i]).atomic_mass
        M=M*amass**st[i]

    M=M**(1/N)*1.66053904e-27
    [V0,E0,B0,BP0,y,T_debye0]=[[1 for i in range(len(cev[0])-1)] for j in range(6)]
    EVfit=[[] for i in range(len(cev[0])-1)]
    [T,G,VT,BT,S,H,TEC,Cp,Cv,BP]=[[[] for i in range(len(cev[0])-1)] for j in range(10)]
    F=[[] for i in range(len(cev[0])-1)]
    for i in range(len(cev[0])-1):
        ds=debye_single(V,E[i],N,M,s,Tmax,nat=nat,eos_only=eos_only,eos=eos)
        if trimve=='T':
            Vtrim,E[i]=trim_ve(V,E[i])
            ds=debye_single(Vtrim,E[i],N,M,s,Tmax,nat=nat,eos_only=eos_only,eos=eos)
        [V0[i],E0[i],B0[i],BP0[i],y[i],T_debye0[i]]=ds[0]
        EVfit[i]=ds[1]
        if eos_only=='False':
            [T[i],G[i],VT[i],BT[i],S[i],H[i],TEC[i],Cp[i],Cv[i],BP[i]]=ds[2]
            F[i]=ds[3]
        
    if os.path.exists(expt):
        data,alias=read_expt(expt)
        
    """plot"""
    start = time.time()
    fts=15 # font size
    if show=='T':
        colors=['b','r','k','m','g','y','c']
        markers=['*','o','s','D','P','^','v']
        properties0=[G,VT,BT,S,H,TEC,Cp,Cv,BP]
        plotnum='abcdefghijklmn'
        if prop!='all':
            properties=[properties0[i] for i in prop]
        else:
            properties=properties0
        if eos_only=='True':
            properties=[]
        nrow=math.ceil((len(properties)+1)/3)
        plt.rcParams["font.size"] = 15
        plt.rcParams["font.weight"] = "bold"
        plt.rcParams['lines.linewidth'] = 1.5
        plt.rcParams['axes.linewidth'] = 1.5
        plt.rcParams['xtick.major.width'] = 1.5
        plt.rcParams['ytick.major.width'] = 1.5
        plt.subplot(nrow,3,1)
        for i in range(len(cev[0])-1):
            plt.plot(EVfit[i][0],EVfit[i][1],color=colors[i],label=cna[i+1])
            plt.scatter(V,E[i],color=colors[i],marker='o')
        plt.xlabel("V ($\AA^3$)",fontsize=fts,weight='bold')
        plt.ylabel("E (eV)",fontsize=fts,weight='bold')
        plt.tick_params(axis='both',which='major',labelsize=12)
        plt.legend(frameon=False,prop={'size': 12})
        ax=plt.gca()
        plt.text(0.5, 0.9,'('+plotnum[0]+')', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
        plt.subplot(nrow,3,2)
        for i in range(len(cev[0])-1):
            for j in range(0,len(F[i]),5):
                fvfit=VEB0(V,np.array(F[i][j])*N,eos=eos)[1]
                plt.plot(fvfit[0],fvfit[1],color=colors[i])
            plt.plot(np.array(VT[i])*N/(1e-24*NA*nat),np.array(G[i])*N/(q*NA/1e3*nat),color='r')
        plt.xlabel("V ($\AA^3$)",fontsize=fts,weight='bold')
        plt.ylabel("F (eV)",fontsize=fts,weight='bold')
        plt.tick_params(axis='both',which='major',labelsize=12)
        plt.legend(frameon=False,prop={'size': 12})
        ax=plt.gca()
        plt.text(0.5, 0.9,'('+plotnum[1]+')', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
        k=0
        property_names0=['G (kJ/mol)','V ($cm^3$/mol)','B (GPa)','S (J/mol/K)','H (kJ/mol)','TEC ($10^-6$/K)','Cp (J/mol/K)','Cv (J/mol/K)','BP']
        property_symbols0=['G','V','B','S','H','TEC','Cp','Cv','BP']
        if prop!='all':
            property_names=[property_names0[i] for i in prop]
            property_symbols=[property_symbols0[i] for i in prop]
        else:
            property_names=property_names0
            property_symbols=property_symbols0
        for i in properties:    
            plt.subplot(nrow,3,3+k)
            for j in range(len(cev[0])-1):
                plt.plot(T[j][0:len(T[j])-2],i[j][0:len(T[j])-2],color=colors[j],label=cna[j+1])
            if os.path.exists(expt):
                ref_set=[]
                ref_set=[data[j]['ref'] for j in range(len(data)) if not data[j]['ref'] in ref_set]
                ref_id={}
                for j in range(len(ref_set)):
                    ref_id[ref_set[j]]=j
                for j in range(len(data)):
                    if data[j]['property']==property_symbols[k]:
                        plt.scatter(data[j]['T'],data[j]['value'],color='orangered',marker=markers[int(data[j]['ref'])],label='Expt., '+alias[int(data[j]['ref'])-1])
            plt.xlabel("T (K)",fontsize=fts,weight='bold')
            plt.ylabel(property_names[k],fontsize=fts,weight='bold')
            plt.xlim(0,T[0][len(T[0])-3])
            plt.tick_params(axis='both',which='major',labelsize=12)
            plt.legend(frameon=False,prop={'size': 12})
            ax=plt.gca()
            plt.text(0.5, 0.9,'('+plotnum[2+k]+')', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
            k=k+1
        plt.rcParams["axes.labelweight"] = "bold"
        plt.subplots_adjust(left=0.085,bottom=0.1,top=0.95,right=0.99,hspace=0.3,wspace=0.3)
        fig = plt.gcf()
        plt.show()
        fig.savefig('thermo.png',dpi=600,transparent=True)

    if write=='T':
        rf=5
        for i in [G,VT,BT,S,H,TEC,Cp,Cv,BP]:
            for k in range(len(cev[0])-1):
                i[k]=[round(i[k][j],rf) for j in range(0,len(i[k]))]

        with open(out,'w') as f:
            for i in range(len(cev[0])-1):            
                f.write("EOS fitting results for the dataset %g\n"%(i+1))
                f.write('\n')
                f.write("%-15s %-15s %-15s %-15s %-15s %-15s\n"%('V0 (angstrom^3)','E0 (eV)','B0 (GPa)','BP0','y','T_debye0'))
                f.write("%-15g %-.13g %-15g %-15g %-15g %-15g\n"%(V0[i],E0[i],B0[i],BP0[i],y[i],T_debye0[i]))
                f.write('\n')
                f.write("%-15s %-15s %-15s %-15s %-15s %-15s\n"%('V0 (cm^3/mol)','E0 (kJ/mol)','B0 (GPa)','BP0','y','T_debye0'))
                f.write("%-15g %-.13g %-15g %-15g %-15g %-15g\n"%(V0[i]/N*1e-24*NA*nat,E0[i]/N*q*NA/1e3*nat,B0[i],BP0[i],y[i],T_debye0[i]))
                f.write('\n')
                f.write('Thermodynamic properties based on the Debye model (groundstate static energy at 0 K is substracted)\n')
                f.write('\n')
                f.write("%-8s %-15s %-15s %-10s %-15s %-15s %-15s %-15s %-15s %-8s\n"%('T (K)','G (kJ/mol)','V (cm^3/mol)','B (GPa)','S (J/mol/K)','H (kJ/mol)','TEC (1e-6/K)','Cp (J/mol/K)','Cv (J/mol/K)','BP'))
                for j in range(0,len(T[i])-2,Tstep_write):
                    f.write("%-8g %-15g %-15g %-10g %-15g %-15g %-15g %-15g %-15g %-8g\n"%(T[i][j],G[i][j],VT[i][j],BT[i][j],S[i][j],H[i][j],TEC[i][j],Cp[i][j],Cv[i][j],BP[i][j]))

        end = time.time()
        print('Time for plotting and writing:'+str(end-start))

    return [V0,E0,B0,BP0,y,T_debye0],T,G,VT,BT,S,H,TEC,Cp,Cv,BP    

