#!/usr/bin/env python
import sys,os
reduce_home=os.path.abspath(os.path.dirname(sys.argv[0])+"/..")
sys.path.append(reduce_home)
import bmxreduce as br
dm = br.datamanager()
dm.setTags(new=True)
dm.setPassages()


##
## let's wait with more sophisticated argument parsing for later
command=sys.argv[1] if len(sys.argv)>1 else "list"

def get_cmd(tags):
    return "bin/reduce_do.py -x pas -n %s -t %s"%(tags[0],",".join(tags))

def write_condor_job(names, tags, fn):
    f=open(fn,'w')
    f.write("""
Universe        = vanilla
Executable      = //direct/astro+u/bmx/bmxreduce/bin/reduce_do.py
request_memory = 17000M
request_cpus = 1
Priority        = 4
GetEnv          = True
Initialdir      = /direct/astro+u/bmx/bmxreduce
Input           = /dev/null 
""")
    for n,t in zip(names,tags):
        f.write ('Arguments       = "-n %s -t %s" \n'%(n,t))
        f.write ('Output          =  /astro/u/bmx/jobfiles/BMX_%s.log \n'%n)
        f.write ('Error           = /astro/u/bmx/jobfiles/BMX_%s.err \n'%n)
        f.write ('Queue\n')
    f.close()



print ("\n--------------\n")
if command=='list':
    print ("Passages: ")
    for p in dm.passages():
        print ("         %s : %i tags."%(p[0], len(p)))
    print ("Fragments: ")
    for f in dm.fragments():
        print ("         %s : %i tags."%(f[0], len(f)))
elif command=="reduce":
    try:
        pas=sys.argv[2]
    except:
        print ('specify passage')
        sys.exit(1)
    for p in dm.passages():
        if pas==p[0]:
            print ("Passage ",pas,"reduced before to stage:",dm.reduce_stage('pas',pas))
            command=get_cmd(p)
            print ("Executing: ",command)
            os.system(command)
            break
elif command=="cron":
    paslist=[]
    for pas in dm.passages():
        pstage=dm.reduce_stage('pas',pas[0])
        print ("Passage ",pas[0],"reduced before to stage:",pstage)
        if (pstage>=0) and (pstage<2):
            paslist.append(pas)

    names=[p[0] for p in paslist]
    tags=[",".join(p) for p in paslist]
    if len(names)>0:
        write_condor_job(names,tags,'/astro/u/bmx/jobfiles/cron.job')
        os.system("condor_submit /astro/u/bmx/jobfiles/cron.job")
    else:
        print ("Nothing to do.")
else:
    print ("Command %s not understood."%command)

    

