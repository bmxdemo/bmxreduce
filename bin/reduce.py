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
    cmdlist=[]
    for pas in dm.passages():
        pstage=dm.reduce_stage('pas',pas[0])
        print ("Passage ",pas[0],"reduced before to stage:",pstage)
        if (pstage<1):
            command=get_cmd(pas)
            cmdlist.append(command)



    sys.exit(1)
    f = br.farmit.farmit('bin/reduce_batch.py', args={'t':cmdlist},
                  names=['BMX_'+tag for tag in dm.tags],
                  reqs={'N':2,'X':0,'priority':'low','mode':'bycore1'})
    f.writejobfiles()
    #f.runjobs(maxjobs=500)
            
else:
    print ("Command %s not understood."%command)

    

