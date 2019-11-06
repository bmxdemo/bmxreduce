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
    return "bin/reduce_do.py pas %s %s"%(tags[0],",".join(tags))


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
        if (pstage<1):
            paslist.append(pas)

    names=[p[0] for p in paslist]
    tags=[",".join(p) for p in paslist]
    f = br.farmit ('bin/reduce_do.py', args={'n':names,'t':tags},
                         names=['BMX_'+n for n in names],
                  reqs={'N':1,'X':0,'priority':'low','mode':'bynode'})
    f.writejobfiles()
    f.runjobs(maxjobs=10)
            
else:
    print ("Command %s not understood."%command)

    

