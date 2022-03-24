#this is what we had until now
    

from __future__ import division
import os
import subprocess as subp
import numpy as np
import sys
import time
import functions_fitting_2018_03_16 as ffit
import itertools


i=int(sys.argv[1])

np.random.seed(i)

program_name="electricG_v1.c"
parlist_names_ordered=["Krg", "nrg", "alpha_Rg", "delta_Rg", "alpha_rRg", "a_fl", "a_dif", "gamma_dif", "gamma_fl", "std_ic","Gi01", "Gi02", "Kgt", "dperturbK", "tperturbK", "K_perturb", "tperturbG", "GE_perturb", "V0_tht", "g_tht", "dl", "alpha_t", "gamma_t",  "stopflow1", "stopflow2", "dx", "tt", "dt", "st", "zero_tol", "gv", "Radius2", "Radius1", "freeze", "mk", "Ikmax", "gl", "gk", "fl_Gext", "V0", "Dg", "F", "Pgrow", "alpha_gt", "freeze2", "tprint", "Dp", "Sth", "K_media", "a0", "GE", "center1", "delta_g", "center2", "seed_", "GS0", "b", "read_input", "S0", "Vk", "Vl", "Dke", "unfreeze", "gamma_s", "ug", "fl_Ekext", "Ls", "dist_pK", "alpha_sv"]






def norm(x):
    return (x-x.min())/(x.max()-x.min())

   
pars_setup={
    "Ls":200,
    "Radius1":np.random.uniform(15,25),
    "center1":0,
    "Radius2":0,
    "center2":200,
    "freeze":0,
    "unfreeze":0.1,
    "freeze2":50,
    "dx":1,
    "st":0.05,
    "tt":24,
    "zero_tol":0.05,
    "dt":0.000005,
    "read_input":-1.0,
    "tprint":50,
    "seed_":i,
}
pars_flow={
    "fl_Gext":5,
    "fl_Ekext":5,
    "K_media":8,
    "GE":30,
    "stopflow1":50,
    "stopflow2":50,
    "tperturbG":50,
    "tperturbK":30,
    "dperturbK":0.05,
    "K_perturb":30,
    "dist_pK":0,
    "GE_perturb":0.5,
}


pars_list=[('Dg', 36000), ('alpha_gt', 24), ('Kgt', 0.75), ('gv', 1), ('V0', -150), 
           ('delta_g', 4.8), ('g_dg', 0.5), ('Vth_dg', -150), 
           ('S0', 1.12), ('gamma_s', 2.8), ('GS0', 0.2), ('ug', 2),('alpha_sv',0.0),
           ('gl', 18), ('gk', 70), ('dl', 4), ('Vk', 100), ('Vl', -156), 
           ('F', 0.05), ('Dp', 0.12), ('Ikmax', 300), ('Dke', 72000), 
           ('a0', 91), ('Sth', 0.03), ('b', 34), ('mk', 1), 
           ('g_tht', 0.3), ('alpha_t', 20.), ('gamma_t', 10.), ('V0_tht', -165), 
           ('Pgrow', 0.5), ('Gi01', 3), ('Gi02', 0), ('std_ic', 0.001), 
           ('gamma_fl', 0.085), ('gamma_dif', 0.085), ('a_fl', 0.012), ('a_dif', 0.012), 
         ('alpha_Rg', 4.5), ('delta_Rg', 24), ('alpha_rRg', 31), ('Krg', 2.25), ('nrg', 2)]




pars=dict(pars_list)
pars.update(pars_setup)
pars.update(pars_flow)

basal_pars=pars.copy()


conditions=['WT30','WT15','TrkA4','GltA']
#conditions=['TrkA3']
#conditions=['TrkA4'] #added 28th may

cond=conditions[(i-1)%len(conditions)]
#cond=conditions[0] #28th may for one single condition
if True:
    pars=basal_pars.copy()
    if cond=='WT30':
        pars['GE']=30
    elif cond=='WT15':
        pars['GE']=15
    elif cond=='TrkA4':
        pars['a0']=pars['a0']*2 #1.5
        pars['gk']=pars['gk']*0.5 #0.75
    elif cond=='GltA':
        pars['delta_g']=pars['delta_g']*1.5
        pars['alpha_gt']=pars['alpha_gt']*1.3

    np.random.seed(i)

    outf_log=(open('log_%d_%s.log'%(i,cond),'w'))


    infile="parset_%i_%s.txt"%(i,cond)
    inf=open(infile,"w")
    for pnum in range(len(parlist_names_ordered)):
        parname=parlist_names_ordered[pnum]
        parval=pars[parname]
        print >> inf, "%s %g"%(parname,parval)
    inf.close()

    output_dir="out_%d_%s"%(i,cond)   
    outfile_name="out_name"
    basename=outfile_name
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    

    ic_dir="ic"
    ic_name="ic"
    #ic_dir=""
    #ic_name=""

    toprintinputdir="ic_print"
    if not os.path.isdir(toprintinputdir):
        os.mkdir(toprintinputdir)
    toprintinputdir_name="ic_print_name"




    inittime=time.time()
    run=False
    try:
        subp.check_call('./%s %s %s %s %s %s %s %s' %(program_name.strip('.c'), output_dir, outfile_name, ic_dir, ic_name, toprintinputdir, toprintinputdir_name, infile), shell=True)
        run=True
        print >> outf_log, 'run OK'

    except(subp.CalledProcessError) as error:
        wenttonan=(error.returncode==5) #if it went to Nan, it exited with code 5
        if wenttonan:
            print >> outf_log, 'went to nan'
        else:
            print >> outf_log, 'other things went wrong', error
        #print 'done in', time.time()-inittime,'s'

    #curval=pars[partochange]
    if run is True:
        varlist=["c","Ek"]
        rs=ffit.parse_results_intonumpy_dtype(output_dir=output_dir,basename=outfile_name,varlist=varlist,dtype='d')
        kym=rs[1]
           
        stpoints=np.arange(len(kym))
        edges=ffit.return_edges(rs[0])
        Ek_timetrace=kym[stpoints,edges[:,1]-5]
        np.savetxt('Ek_%i_%s.txt'%(i,cond),Ek_timetrace)
        np.savetxt('edge_%i_%s.txt'%(i,cond),edges[:,1])
    

    outf_log.flush()
 

    outf_log.close()






