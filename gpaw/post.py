import numpy as np
import os

def vol(ftraj='atoms.traj'):
    from ase.io import read
    a=read(ftraj,':')
    v=a[-1].get_volume()
    with open('vol','w') as f:
        f.write(str(v))
    return v

# do eos fitting and Debye model for one functional and save results in files

def eos_gpaw(ev='e-v.dat',struc='rel.in',eos_only='False'):
    from gthermo_open import debye
    try:
        with open(ev, 'r') as f:
            c=f.readlines()
        cc=[]
        for i in c:
            if not '#' in i and i!='\n' and i!=' \n':
                cc.append(i.strip('\n'))

        v=[1 for i in cc]
        e=[0 for i in cc]
        for i in range(len(cc)):
            cc[i]=cc[i].split(' ')
            for j in range(cc[i].count('')):
                cc[i].remove('')
            v[i]=float(cc[i][0])
            e[i]=float(cc[i][1])

        emin=min(e)
        for i in range(len(e)):
            if e[i]==emin:
                imin=i

        v=v[max(0,imin-2):min(len(v),imin+3)]
        e=e[max(0,imin-2):min(len(e),imin+3)]
        with open(ev,'w') as ff:
            ff.write('\n')
            for i in range(len(v)):
                ff.write(str(v[i])+'   '+str(e[i])+'\n')

        D=debye(expt='noread',ev=ev,Tmax=1000,Tstep_write=10,struc=struc,show='F',write='F',eos_only=eos_only)
        name=['T','G','VT','BT','S','H','TEC','Cp','Cv','BP']
        name2=['V0','E0','B0','BP0','y','T_debye0']
        if eos_only=='False':    
            for j in range(1,len(D)):
                with open('beef/'+name[j-1],'a') as f:
                    f.write(str(D[j][0])[1:-1]+' ')
                    f.write('\n')
        for j in range(0,len(D[0])):
            with open('beef/'+name2[j],'a') as f3:
                f3.write(str(D[0][j][0])+' ')
                f3.write('\n')
        with open('beef/'+'cunt','r') as fcunt:
            ccunt=fcunt.readlines()
        ncunt=len(ccunt)+1
        with open('beef/'+'cunt','a') as f2:
            f2.write(str(ncunt)+'\n')
    except:
        print('EOS fails!')
        name=['T','G','VT','BT','S','H','TEC','Cp','Cv','BP']
        name2=['V0','E0','B0','BP0','y','T_debye0']
        if eos_only=='False':
            for j in range(len(name)):
                with open('beef/'+name[j],'a') as f:
                    f.write('\n')
        for j in range(len(name2)):
            with open('beef/'+name2[j],'a') as f3:
                f3.write('\n')
        with open('beef/'+'cunt','a') as f2:
            f2.write('\n')

# calculate probability distribution function (PDF) from property ensembles by BEEF  

def pdf(dir='beef/'):
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    name=['G','VT','BT','S','H','TEC','Cp','Cv','BP'] 
    for j in range(len(name)):
        if os.path.exists(dir+'/'+name[j]):
            a=np.loadtxt(dir+'/'+name[j],delimiter=',')
            if os.path.exists(dir+'/Cp'):
                acp=np.loadtxt(dir+'/Cp',delimiter=',')
            else:
                acp=a
            aa=np.append(a[:,0],a[:,-1])
            maxx=np.ceil(max(aa)/10)*10
            minn=np.floor(min(aa)/10)*10
            if j==3 or j==4 or j==6 or j==7:
                minn=0
            ystep=min((maxx-minn)/100,1) 
            bin=[minn+ystep*i for i in range(int((maxx-minn)/ystep+1))]
            y=[i+ystep/2 for i in bin[0:-1]]
            y=np.array(y)
            np.savetxt(dir+'/ybin_'+name[j],y)
            n=[[] for i in range(acp.shape[1])]
            for i in range(acp.shape[1]):
                with open(dir+'/count','a') as ff:
                    ff.write(str(i)+' ')
                n[i], bins, patches = plt.hist(a[:,i], bin, density=True,histtype='step')
            np.savetxt(dir+'/n_'+name[j],n)

# plot histograms of properties at a given temperature 

def plothist(dir='beef/',T=300,nat=1):
    import time
    import os
    import matplotlib
    import matplotlib.pyplot as plt
    name=['G','VT','BT','S','H','TEC','Cp','Cv','BP0','y','T_debye0']
    acp=np.loadtxt(dir+'/Cp',delimiter=',')
    t=np.loadtxt(dir+'/T',delimiter=',')
    x=t[0,:][0:-2]
    start = time.time()

    #filter outliners
    cp=np.loadtxt(dir+'/'+'Cp',delimiter=',')
    cp=cp[:,-1]
    nfil=[]
    for i in range(len(cp)):
        if abs(cp[i])>500 or cp[i]<=0:
            nfil.append(i)
        print(nfil)

    for j in range(len(name)):
        if j==2 or j==5 or j>=8:
            naty=1
        else:
            naty=nat
        if os.path.exists(dir+'/'+name[j]):
            if j>=8:
                a=naty*np.loadtxt(dir+'/'+name[j],delimiter=',')
            else:
                ilook=int(T/(x[1]-x[0]))
                a=naty*np.loadtxt(dir+'/'+name[j],delimiter=',')[:,ilook]
            a=[a[m] for m in range(len(a)) if m not in nfil]
            plt.subplot(3,4,1+j)
            plt.hist(a, 100, density=True,histtype='step')
            if j==1:
                plt.title('T='+str(T)+' K',loc='right')
            plt.xlabel(name[j])
            np.savetxt(dir+name[j]+'_'+str(T),a)

    end = time.time()
    print('Time for histograms:'+str(end-start))
    plt.subplots_adjust(left=0.075,bottom=0.1,top=0.9,right=0.99,hspace=0.3,wspace=0.35)
    plt.show()

# plot PDF of properties

def plotpdf(dir='beef/',expt='expt',nat=1,prop='all'):
    import matplotlib.pyplot as plt
    import os
    import math
    from matplotlib import ticker,cm
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    from gthermo_open import debye,read_expt

    cna=['PBE','RPBE','optPBE-vdW','PBEsol','PW91'] # functionals used to generate data in e-v-multi.dat if existing

    ev='e-v-multi.dat'
    if os.path.exists(ev):
        de=debye(show='F',write='F',trimve='F',ev=ev,struc='POSCAR')

    cmap=cm.get_cmap('hot_r')

    if os.path.exists(expt):
        data=read_expt(expt)[0]
    
    lines=[]
    legends=[]
    legendcount=[]
    name0=['G','VT','BT','S','H','TEC','Cp','Cv']
    property_symbols0=['G','V','B','S','H','TEC','Cp','Cv']
    if prop!='all':
        name=[name0[i] for i in prop]
        property_symbols=[property_symbols0[i] for i in prop]
    else:
        name=name0
        property_symbols=property_symbols0
    ncol=min(3,len(name))
    nrow=math.ceil(len(name)/ncol)
    markers=['*','o','s','D','P','^','v']
    n=[[] for j in range(len(name))]
    y=[[] for j in range(len(name))]
    t=np.loadtxt(dir+'/T',delimiter=',')
    x=t[0,:][0:-2]
    for j in range(len(name)):
        print(name[j])
        if os.path.exists(dir+'/n_'+name[j]):
            if j==2 or j==5:
                naty=1
            else:
                naty=nat
            y[j]=naty*np.loadtxt(dir+'/ybin_'+name[j])
            X,Y=np.meshgrid(x,y[j])
            n[j]=np.loadtxt(dir+'/n_'+name[j])
            n[j]=n[j].transpose()
            plt.subplot(nrow,ncol,1+j)
            plt.contourf(X,Y,n[j], cmap=cmap)
            if os.path.exists(ev):
                for k in range(len(de[j+2])):
                    pcalc,=plt.plot(de[1][k][0:len(de[1][k])-2],np.array(de[j+2][k][0:len(de[1][k])-2]),label=cna[k])
                    if j==0:
                        lines.append(pcalc)
                        legends.append(cna[k])
            if os.path.exists(expt):
                for k in range(len(data)):
                    if data[k]['property']==property_symbols[j]:
                        pexpt=plt.scatter(data[k]['T'],data[k]['value'],color='blue',marker=markers[int(data[k]['ref'])],label='Expt.['+data[k]['ref']+']')
                        if not data[k]['ref'] in legendcount:
                            lines.append(pexpt)
                            legends.append('Expt.['+data[k]['ref']+']')
                            legendcount.append(data[k]['ref'])
            plt.axis([min(x),max(x),min(y[j]),max(y[j])])
            plt.xlabel("T (K)")
            plt.ylabel(name[j])
    plt.figlegend(lines,legends,'lower right')
    plt.subplots_adjust(left=0.075,bottom=0.1,top=0.9,right=0.99,hspace=0.3,wspace=0.35)
    fig = plt.gcf()
    plt.show()

# plot property ensembles

def plot(dir='beef/',expt='expt',nat=1):
    import matplotlib.pyplot as plt
    import os
    from gthermo_open import read_expt

    if os.path.exists(expt):
        data=read_expt(expt)

    name=['G','VT','BT','S','H','TEC','Cp','Cv']
    property_symbols=['G','V','B','S','H','TEC','Cp','Cv']
    markers=['*','o','s','D','P','^','v']
    y=[[] for j in range(len(name))]
    t=np.loadtxt(dir+'/T',delimiter=',')
    x=t[0,:][0:-2]
    xlen=len(x)
    for j in range(len(name)):
        if j==2 or j==5:
            naty=1
        else:
            naty=nat
        y[j]=naty*np.loadtxt(dir+'/'+name[j],delimiter=',')
        y[j]=y[j][:,0:xlen]
        plt.subplot(241+j)
        for k in range(y[j].shape[0]):
            plt.plot(x,y[j][k],'b')

        if os.path.exists(expt):
            for k in range(len(data)):
                if data[k]['property']==property_symbols[j]:
                    plt.scatter(data[k]['T'],data[k]['value'],color='orangered',marker=markers[int(data[k]['ref'])],label='Expt.['+data[k]['ref']+']')
                    plt.legend()
        plt.axis([min(x),max(x),y[j].min(),y[j].max()])
        plt.xlabel("T (K)")
        plt.ylabel(name[j])
    plt.subplots_adjust(left=0.075,bottom=0.1,top=0.9,right=0.99,hspace=0.3,wspace=0.35)
    plt.show()

# calculate, plot and save statistical quantities of property ensembles

def stat(dir='beef/',expt='expt',nat=1):
    import matplotlib.pyplot as plt
    import os
    from gthermo_open import read_expt
    import numpy as np
    from scipy import stats

    if os.path.exists(expt):
        data=read_expt(expt)

    cna=['mu','sigma','COV','S','K','JB']
    #stat=[mu,sigma,skew,kurt,jb]
    name=['G','VT','BT','S','H','TEC','Cp','Cv']
    property_symbols=['G','V','B','S','H','TEC','Cp','Cv']
    markers=['*','o','s','D','P','^','v']
    y=[[] for j in range(len(name))]
    t=np.loadtxt(dir+'/T',delimiter=',')
    x=t[0,:][0:-2]
    xlen=len(x)

    #filter outliners
    cp=np.loadtxt(dir+'/'+'Cp',delimiter=',')
    cp=cp[:,-1]
    nfil=[]
    for i in range(len(cp)):
        if abs(cp[i])>500 or cp[i]<=0:
            nfil.append(i)
    print(nfil)
    for j in range(len(name)):
        if j==2 or j==5:
            naty=1
        else:
            naty=nat
        y[j]=naty*np.loadtxt(dir+'/'+name[j],delimiter=',')
        y[j]=y[j][:,0:xlen]
        mu=[0 for j in range(y[j].shape[1])]
        sigma=[0 for j in range(y[j].shape[1])]
        cov=[0 for j in range(y[j].shape[1])]
        skew=[0 for j in range(y[j].shape[1])]
        kurt=[0 for j in range(y[j].shape[1])]
        jb=[0 for j in range(y[j].shape[1])]
        plt.subplot(241+j)
        for k in range(y[j].shape[1]):
            xstat=list(y[j][:,k])
            xstat=[xstat[m] for m in range(len(xstat)) if m not in nfil]
            jb[k]=stats.jarque_bera(xstat)[0]
            skew[k]=stats.skew(xstat)
            kurt[k]=stats.kurtosis(xstat)
            cov[k]=stats.variation(xstat)
            mu[k]=np.average(xstat)        
            sigma[k]=np.std(xstat)
        plt.xlabel("T (K)")
        plt.ylabel(name[j])
        with open(dir+'stat_'+name[j],'w') as f:
            f.write("%-8s %-15s %-15s %-15s %-15s %-15s\n"%('#T',cna[0],cna[1],cna[2],cna[3],cna[4]))
            for i in range(len(x)):
                f.write("%-8g %-15g %-15g %-15g %-15g %-15g\n"%(x[i],mu[i],sigma[i],cov[i],skew[i],kurt[i])) 
    plt.legend()
    plt.subplots_adjust(left=0.075,bottom=0.1,top=0.9,right=0.99,hspace=0.3,wspace=0.35)
    plt.show()
