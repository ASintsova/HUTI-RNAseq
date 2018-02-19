import os
import logging

def generatePBSScript(script_name, script):

    pbs_preamble = """
    ####  PBS preamble

    #PBS -N {}

    #PBS -M annasint@umich.edu
    #PBS -m ea
    #PBS -j oe

    #PBS -l nodes=1:ppn=4,pmem=4gb,walltime=12:00:00
    #PBS -V

    #PBS -A hmobley_fluxod
    #PBS -l qos=flux
    #PBS -q fluxod

    ####  End PBS preamble
    ####  Commands follow this line


    if [ -s "$PBS_NODEFILE" ] ; then
        echo "Running on"
        cat $PBS_NODEFILE
    fi

    if [ -d "$PBS_O_WORKDIR" ] ; then
        cd $PBS_O_WORKDIR
        echo "Running from $PBS_O_WORKDIR"
    fi

    """.format(os.path.basename(script_name).split(".")[0])

    pbs = open(script_name, "w")
    pbs.write(pbs_preamble + "\n" + script + "\n")
    pbs.close()

### Test
#import datetime
#script_name = datetime.datetime.now().strftime("%Y-%m-%d")+ '_trimmomatic.pbs'
#script = "TEST TEST TEST"
#generatePBSScript(script_name, script)
