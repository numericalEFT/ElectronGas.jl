#!/usr/bin/python
import random
import os
import sys

##### Modify parameters here  ###############
# Cluster="PBS"
Cluster = "local"

############################################
inlist = open("./inlist", "r")
rootdir = os.getcwd()
execute = "gapfunc_MS.jl"
execute1 = "eigen.jl"
# assert len(sys.argv)==2, "Number of jobs is needed."
# Number=int(sys.argv[1])

for index, eachline in enumerate(inlist):
    para = eachline.split()
    if len(para) == 0:
        print("All submitted!")
        break
    rs = float(para[0])
    beta0 = float(para[1])
    num = int(para[2])
    minchan = int(para[3])
    maxchan = int(para[4])
    dim = int(para[5])
    sigmatype = int(para[-1])

    print("Creating {0} jobs...".format(maxchan-minchan+1))
    fname = "rs{0}_sigtype{1}".format(rs, sigmatype)
    homedir = os.getcwd() + "/"+fname

    if(os.path.exists(homedir) != True):
        os.system("mkdir "+homedir)
    # else:
        # print(homedir+" alreadly exists!")
        # break

    os.system("cp {0} {1}".format(execute, homedir))
    os.system("cp {0} {1}".format(execute1, homedir))

    outfilepath = homedir+"/outfile"
    if(os.path.exists(outfilepath) != True):
        os.system("mkdir "+outfilepath)
    jobfilepath = homedir+"/jobfile"
    if(os.path.exists(jobfilepath) != True):
        os.system("mkdir "+jobfilepath)

    for pid in range(minchan, maxchan+1):

        ### terminal output goes here #############
        outfile = "_out"+str(pid)
        ### job file to submit to cluster  ########
        jobfile = "_job"+str(pid)+".sh"

        pname = "para"+str(pid)
        with open(homedir+"/"+pname, "w") as file:
            parameters = ' '.join(para[:3])+' '+str(dim)+' '+str(pid)
            file.write(parameters+"\n\n")
            file.write("#rs, Beta0, num_beta, dim, channel")

        if Cluster == "local":
            os.chdir(homedir)
            os.system("julia "+execute+' '+str(sigmatype) + ' ' + pname + " &")
            os.chdir("..")
        elif Cluster == "PBS":
            with open(jobfilepath+"/"+jobfile, "w") as fjob:
                fjob.write("#!/bin/sh\n"+"#PBS -N "+jobfile+"\n")
                fjob.write("#PBS -o "+homedir+"/Output\n")
                fjob.write("#PBS -e "+homedir+"/Error\n")
                fjob.write("#PBS -l walltime=2000:00:00\n")
                fjob.write("echo $PBS_JOBID >>"+homedir+"/id_job.log\n")
                fjob.write("cd "+homedir+"\n")
                fjob.write("julia "+execute+' ' + str(sigmatype) + ' ' + pname)

            os.chdir(homedir)
            os.system("qsub "+jobfilepath + "/"+jobfile)
            os.system("rm "+jobfilepath + "/"+jobfile)
            os.chdir("..")
        else:
            print("I don't know what is {0}".format(Cluster))
            break

print("Jobs manage daemon is ended")
sys.exit(0)
