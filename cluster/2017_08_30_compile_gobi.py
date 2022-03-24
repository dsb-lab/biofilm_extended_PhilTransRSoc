import subprocess as subp
import time
def compile_only(program_name="name.c", basicmf="basic_mf_2017_gobi"): 
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
    inittime=time.time()
    subp.check_call("make", shell=True)
    print "compilation done in", time.time()-inittime,"s"

compile_only('electricG_v1.c')

