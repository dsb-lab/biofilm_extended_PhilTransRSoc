from __future__ import division
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import math
import subprocess as subp

def compile_only(program_name="name.c", basicmf="basic_mf_2017"): 
    #done in 17th august 2017, to compile only
    
    Mf=open(basicmf,"r")
    lines=Mf.readlines()
    for lnum,l in enumerate(lines):
        if l.startswith("PROG"):
            print l
            l="PROG = %s\n"%program_name.strip(".c")
            lines[lnum]=l
    Mf.close()
    print lines
    Mf=open("Makefile","w")
    for l in lines:
        Mf.write(l)
    Mf.close()

    subp.call("make clean", shell=True) #removes *.o files (RMC)
    subp.call("rm %s" % program_name.strip(".c"), shell=True)
    #inittime=time.time()
    subp.check_call("make", shell=True)
    #print "compilation done in", time.time()-inittime,"s"
    
def fill_zeros_with_last(arr):
    prev = np.arange(len(arr))
    prev[arr == 0] = 0
    prev = np.maximum.accumulate(prev)
    return arr[prev]

def identify_peaks_histeresis(data,timev=None,direction='+1',threshold=None,tolerance_p=0.05):
    
    mean=data.mean()
    
    tolerance=mean*tolerance_p
    
    if threshold is None:
        threshold=mean
    
    a = np.zeros_like(data)
    if timev is None:
        timev=np.arange(len(data))
    a[ data < threshold-tolerance] = -1
    a[ data > threshold+tolerance] = +1
    a=fill_zeros_with_last(a) #now a is a vector where -1 is "down" and +1 is "up"
    diff=np.diff(a)
    if direction=='+1':
        beginnings=np.where(diff==2)[0] #from -1 to 1 (so that 1-(-1)=2)
        ends=np.where(diff==-2)[0] #from 1 to -1 (so that -1-1=-2)
    elif direction=='-1':
        beginnings=np.where(diff==-2)[0]
        ends=np.where(diff==2)[0]
    else:
        print "direction of pulse not understood"
        return None
    if len(beginnings)==0: #no pulse
        return None
    elif len(ends)==0: #the pulse doesn't end...
        return None
    else:
       
        timeb=timev[beginnings] #time points at which pulses begin
        timee=timev[ends] #time points at which pulses ends
        if timeb[0]>timee[0]:#means the data started in a pulse
            if len(timeb)==1: #if there is no other beginning, then it cannot be
                return None
            else:
                timee=timee[1:]
                
        if timee[-1]<timeb[-1]:#data ends with pulse
            timeb=timeb[0:-1]
            
    return [timeb,timee]
def parse_results_intonumpy_dtype(output_dir=None,basename=None,varlist=None,dtype=None,ndim=1):
    results=[]
    for vnum,variable in enumerate(varlist): #changed into c,g RMC
        if type(dtype)==str:
        	dtype_=dtype
        else:
            dtype_=dtype[vnum]
            
        print "Getting %s" % variable
        var_basename = "{b}_{v}".format(b=basename, v=variable)
        fnin="{o}/{vb}.pat".format(o=output_dir,vb=var_basename)
        fd = open(fnin, 'rb')
        L, = np.fromfile(fd,'I',1)
        npat, = np.fromfile(fd,'I',1) #number of timesteps
        if ndim==2:
            
            m=np.zeros((npat,L,L))
            reshape=[L,L]
            len_=L*L
        else:
            m=np.zeros((npat,L))
            reshape=[L]
            len_=L
            
        print '#',L, npat
        for i in range(0,npat):
            #if True:
            try:
                t_seq, = np.fromfile(fd,'d',1);
                p = np.fromfile(fd,dtype_,len_).reshape(reshape)
                m[i]=p
            except:
                print 'except,',i,
                m=m[:i]
        results.append(m)
    return results
def return_edges(kym):
    cdenlim=0.01
    mask=np.array(np.ma.greater(kym,cdenlim))
    rows=len(kym)
    if kym[0][0]>cdenlim or kym[0][-1]:
        wall=True
        print 'wall!'
    else:
        wall=False
    diffs=np.where(np.diff(mask))
    if wall is False:
        if diffs[0][2]==diffs[0][1]:
            nbiofilms=2
            cols=4
        else:
            nbiofilms=1
            cols=2
        edges=np.array(diffs[1]).reshape(rows,cols)
        #now the cell that is biofilm and is at the edge is actually one more from the left
        edges[:,0]+=1 #
        if nbiofilms==2:
            edges[:,2]+=1 
    else:
        
        if diffs[0][1]==diffs[0][0]:
            nbiofilms=2
            cols=4
            edges_=np.array(diffs[1]).reshape(rows,2)
            
        else:
            nbiofilms=1
            cols=2
            edges_=np.array(diffs[1])
        edges=np.zeros((rows,cols),dtype='int')
        if nbiofilms==1:
            edges[:,1]=edges_
        else:
            edges[:,1:3]=edges_
            edges[:,3]=np.ones(rows)*len(kym[0])
    
    
    return edges




def plot_all(outdir='',basename='', st=1, color_gr1='', color_gr2='',per_int=5, threshold_Ek=9, wsize_tht=1, wsize_gr=1, plotkym=True,plotEk=False, plotThT=False,varlist=None):
    rsc=parse_results_intonumpy_dtype(output_dir=outdir,basename=basename,varlist=['c'],dtype='d')
    #make C kymograph
    plt.figure(figsize=(3,3))
    kym=rsc[0]
    plt.imshow(kym[::-1],interpolation='None')
    plt.grid()

    #get edges
    edges=return_edges(kym)
    stpoints=np.arange(len(rsc[0])) #number of printed time points
    time_=stpoints*st #time vector in simulation time units

    #plot edges on kymograph
    for i in range(len(edges[0])): #edges of one or two biofilms
        plt.scatter(edges[:,i][::-1],stpoints,s=0.5,color='white')
    #plt.xlim([250,350])
    #plt.ylim([0,60])
    plt.show()


    #parse and plot Ek and ThT (or any other variable)

    if varlist is None:
        varlist=['g','gi','Rg','s','nk','Ek','Ik','tht','v']
    varlist_peaks=['Ek','tht']
    ncols=3
    nrows=int(math.ceil(len(varlist)/ncols))
    rs=parse_results_intonumpy_dtype(output_dir=outdir,basename=basename,varlist=varlist,dtype='d')
    print 'nrows', nrows
    if plotkym is True:
        fig,axes=plt.subplots(nrows,ncols,figsize=(12,nrows*3))
    
        for vnum,var in enumerate(varlist):
            try:
                ax=axes[vnum//ncols][vnum%ncols]
            except TypeError:
                ax=axes[vnum%ncols]

            kym=rs[vnum]
            im=ax.imshow(kym[::-1],interpolation='None',aspect='auto',cmap='YlGnBu')
            plt.colorbar(im,ax=ax)
            ax.set_title(var)
            #ax.grid()
            for i in range(len(edges[0])): #edges of one or two biofilms
                ax.scatter(edges[:,i][::-1],stpoints,s=0.5,color='k')
            #plt.xlim([250,350])
            #plt.ylim([0,60])
            npat=len(kym)
            L=len(kym[0])
            ax.set_yticks(np.linspace(0,npat,11))
            ax.set_yticklabels(map(str,st*np.linspace(0,npat,11)[::-1]))
            ax.set_ylabel("simulation time")
            ax.set_xticks(range(L)[::100])
            ax.set_xticklabels(map(str,np.arange(0,L)[::100]))
            ax.set_xlabel("lattice position")

        plt.tight_layout()
        plt.show()


    #now get peaks at a given distance from the periphery (per_int), and at the center


    trajectories_overlay=dict()
    #first define the arrays of coordinates at center and periphery
    if len(edges[0])==2:
        centers=[np.ones(len(stpoints),dtype='int')*int(edges[0][0]+(edges[0][1]-edges[0][0])/2)]
        peripheries=[edges[:,1]-per_int]
        nb=1
    else:
        centers=[np.ones(len(stpoints),dtype='int')*int(edges[0][0]+(edges[0][1]-edges[0][0])/2),np.ones(len(stpoints),dtype='int')*int(edges[0][2]+(edges[0][3]-edges[0][2])/2)]
        peripheries=[edges[:,1]-per_int, edges[:,2]+per_int]
        nb=2


    peak_times_d=dict()
    wsize_h=int(wsize_tht/2.)
    print 'wsize is', wsize_tht, 'wsize_h is', wsize_h
    for vnum,var in enumerate(varlist):
        if var in varlist_peaks:
            peak_times_d[var]=[[] for i in range(2*nb)]
            peak_times=peak_times_d[var]
            tt_plot=[] #timetrace to plot with growth rate (peripheries)
            kl=0
            for lnum,lines in enumerate([centers,peripheries]):
                for line in lines:#in case there are two centers or peripheries
                    kym=rs[vnum]
                    timetrace=kym[stpoints,line]
                    timetrace_=timetrace.copy()                    

                    if var=='Ek':
                        threshold=threshold_Ek
                        ttrace=timetrace

                    elif var=='tht':
                        #detrend (actually only necessary at center but it will not change much in periphery anyway)

                        for i in range(len(timetrace_)):
                            if i > wsize_h and i<len(timetrace-wsize_h):
                                timetrace_[i]=timetrace[i]-np.mean(timetrace[i-wsize_h:i+wsize_h])
                            elif i>wsize_h:
                                timetrace_[i]=timetrace[i]-np.mean(timetrace[i-wsize_h:])
                            else:
                                timetrace_[i]=timetrace[i]-np.mean(timetrace[0:i+wsize_h])
                        ttrace=timetrace_
                        threshold=0    
                    if lnum==1:
                        tt_plot.append(ttrace)

                    beg_end=identify_peaks_histeresis(ttrace,threshold=threshold,tolerance_p=0.1) #peak identification based on normalise
                    if beg_end is not None:
                        beg,end=beg_end
                        #print beg, end
                        centerst=beg+(end-beg)/2
                        peak_times[kl]=centerst

                    #plot kymograph and timetrace, with peaks
                    if (var=='Ek' and plotEk) or (var=='tht' and plotThT):
                        fig,axes=plt.subplots(1,2,figsize=(10,3))
                        ax=axes[0]
                        ax.imshow(kym,origin='lower',aspect='auto')
                        #plot edges on kymograph
                        for i in range(len(edges[0])): #edges of one or two biofilms
                            ax.scatter(edges[:,i],stpoints,s=0.5,color='white')
                        ax.scatter(line,stpoints,s=0.5)
                        ax.set_title(var)
                        ax1=axes[1]
                        ax1.plot(time_,timetrace)
                        if var=='tht':
                            ax1.plot(time_,timetrace_)

                        ax1.set_title(var)
                        ax1.set_xlabel('time')
                        ax1.set_ylabel(var)

                        if beg_end is not None:
                            for y in centerst:
                                ax.axhline(y=y)
                                ax1.axvline(x=y*st)

                        plt.tight_layout()
                        plt.show()
                    kl+=1


            #now plot growthrate and timetrace 
            fig2,axes2=plt.subplots(1,2,figsize=(10,4))
            fig2.suptitle('with variable %s'%varlist[vnum])

            ax=axes2[0]
            #plot growth
            width=edges[:,1]-edges[:,0]
            ax.scatter(time_,width,color=color_gr1)
            ax.set_xlabel('time')
            ax.set_ylabel('width')


            if nb==1:
                for peak in peak_times[1]:#periphery of biofilm 1
                    ax.axvline(x=peak*st,color=color_gr1,linestyle='--')
            elif nb==2:
                for peak in peak_times[2]: #periphery of biofilm 1
                    ax.axvline(x=peak*st,color=color_gr1,linestyle='--')
                for peak in peak_times[3]: #periphery of biofilm 2
                    ax.axvline(x=peak*st, color=color_gr2,linestyle='--')

            #plot growth rate
            ax=axes2[1]
            dif=np.diff(edges[:,1])
            gr=np.zeros(len(dif))
            wsize=wsize_gr
            wsize_h=int(wsize/2)
            for i in range(len(dif)):
                if i > wsize_h and i<len(dif)-wsize_h:
                    gr[i]=np.mean(dif[i-wsize_h:i+wsize_h])
                elif i>wsize_h:
                    gr[i]=np.mean(dif[i-wsize_h:])
                else:
                    gr[i]=np.mean(dif[0:i+wsize_h])
            ax.plot(time_[:-1], gr,marker='o',color=color_gr1,ms=2)
            ax.set_xlabel('time')
            ax.set_ylabel('growth rate')
            print 'gr', wsize
            #plot timetrace on twin axis
            ax_=ax.twinx()
            ax_.plot(time_,tt_plot[0],color='orange')
            ax_.set_ylabel('%s (periphery)'%var)


            if len(edges[0])==4:#2 biofilm
                ax=axes2[0]
                width=edges[:,3]-edges[:,2]
                ax.scatter(time_,width,color=color_gr2)
                ax=axes2[1]
                wsize=wsize_gr

                wsize_h=int(wsize/2)
                dif=np.diff(edges[:,2])
                gr=np.zeros(len(dif))
                for i in range(len(dif)):
                    if i > wsize_h and i<len(dif)-wsize_h:
                        gr[i]=np.mean(dif[i-wsize_h:i+wsize_h])
                    elif i>wsize_h:
                        gr[i]=np.mean(dif[i-wsize_h:])
                    else:
                        gr[i]=np.mean(dif[0:i+wsize_h])
                ax.plot(time_[:-1], gr,color=color_gr2,marker='o',ms=2)
                ax_=ax.twinx()
                ax_.plot(time_,tt_plot[1],color='brown')


            if nb==1:
                colors_=['r','orange'] #center, periphery
            else:
                colors_=['r','g','orange','cyan'] #center, center, periphery, periphery
            #for p,peakset in enumerate(peak_times):
            #    for peak in peakset:
            #        ax.axvline(x=peak*st,label=key,color=colors_[p])
            fig2.tight_layout()
            fig2.show()
    return [edges, peak_times_d]
