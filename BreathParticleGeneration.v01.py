#!/home/wayne/bin/pyenv/bin/python3

# perry
###!/home/wayne/bin/pyenv/bin/python3

# seawulf
# module load python/3.8.6
###!/gpfs/software/python/3.8.6/bin/python3

# zagros
###!/mnt/lustre/wayne/bin/python3/bin/python3

import math
import random
import sys

debug = False

def BreathParticleGeneration():
    if debug:
        import pdb;
        pdb.set_trace()

    pymajor = 3
    if (sys.version_info[0] <= 2):
        pymajor = 2

# test code
#    for i, t in enumerate(range(0, 525, 25)):
#        locPrint(pymajor, '{}\t{}\t{}'.format(i, t/100, math.sin(2*math.pi*t/500)))
    ymin = 0.98
    ymax = 1.02
    ydel = ymax - ymin
    zmin = 1.67
    zmax = 1.68
    zdel = zmax - zmin

# One breath cycle: 2.5 exhale
#                   2.5 inhale
#                   5.0 TOTAL

# Version 2: output 5 particles every 5ms
# Each time step is 5ms.  There are currently 20000 time steps.
# I expect more iterations, so I'me going to go to 30000 time steps.
# That will be 150s.

# 
# Verions 1: output 5 particles every .25s
# every 0.25s is 50 iterations
# 5s is 1000 iterations
# currently at 15700 iterations 2021/03/05
# going to go to 20000 iterations which is 0s to 100s inclusive

    locPrint(pymajor, '# location  velocity  "start time"  diameter  density')
    locPrint(pymajor, '#  x y z     u v w')
    for ti in range(0, 10025, 25):
        tm = ti % 500
        if tm >= 0 and tm <= 225:

            for i in range(0,5):
                rydel = random.uniform(0.0, ydel)
                y = ymin + rydel
                rzdel = random.uniform(0.0, zdel)
                z = zmin + rzdel

                locPrint(pymajor, '{}\t{}\t{}\t0.0\t0.0\t0.0\t{}\t10.1e-6\t977.0'.format(
                    0.005, y, z, ti/100.0))

def locPrint(pymajor, mess):
    if pymajor == 2:
        print mess
    else:
        print(mess)
    

            
if __name__ == '__main__':
    BreathParticleGeneration()
