from __future__ import division
import os
import subprocess as subp
import numpy as np
import sys
import time

#functions to be used for the fitting purposes
#I modified the get_timetraces_and_width to return None in the values if the cell density went over the field size

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


def get_periods_and_phasedif(peaks1,peaks2):
    sign1=np.asarray(peaks1,dtype='float') #! careful!! if 'float' is not specified, then it is integers and only integer values of phase difference are obtained!!This was done initially, now repeated like this
    sign2=np.asarray(peaks2,dtype='float')

    if len(sign1)>len(sign2):
        print "***careful! sign1 is not the shortest! "
        temp=sign1
        sign1=sign2
        sign2=temp
    inst_per=sign1[1:]-sign1[:-1]

    period1_ = np.hstack([inst_per, [inst_per[-1]]]) #last period cannot be calculated
    
    inst_per2=sign2[1:]-sign2[:-1]
    period2_ = np.hstack([inst_per2, [inst_per2[-1]]]) #last period cannot be calculated
    
    peak_shifts = (sign1[:,None]-sign2[None,:])
    peak_shifts = peak_shifts[range(len(sign1)),np.abs(peak_shifts).argmin(axis=1)]
    peak_shifts *= 2*np.pi/period1_

    #mdif=peak_shifts.mean()

    peak_shifts = abs(peak_shifts)

    mdif=peak_shifts.mean()
    return [mdif,period1_,period2_]

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
