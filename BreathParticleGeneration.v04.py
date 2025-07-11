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

debug = False

def BreathParticleGeneration():
    if debug:
        import pdb;
        pdb.set_trace()

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

# 
# Verions 1: output 5 particles every .25s
# every 0.25s is 50 iterations
# 5s is 1000 iterations
# currently at 15700 iterations 2021/03/05
# going to go to 20000 iterations which is 0s to 100s inclusive

# Version 2: output 5 particles every 5ms
# Each time step is 5ms.  There are currently 20000 time steps.
# I expect more iterations, so I'me going to go to 30000 time steps.
# That will be 150s.
# For the exhale step, there will be 2.5 / 0.005 = 500 particle released cycles
# release from 0 to 360 time steps for 1.8s
# loop from 0, 0.005, 150.005
# multiply by 1000 and use integer arithmatic: 0, 5, 150005
# The cycle is 5s or 5000 which us used to mod time to compare when
# particles need to be generated.

    pId = 0
    print('# location  velocity  "start time"  diameter  density')
    print('#  x y z     u v w')
    for ti in range(0, 200005, 5):
        tm = ti % 5000
        if tm >= 0 and tm < 2500:

            for i in range(0,5):
                rydel = random.uniform(0.0, ydel)
                y = ymin + rydel
                rzdel = random.uniform(0.0, zdel)
                z = zmin + rzdel

                print('{:.3f}\t{:.15f}\t{:.15f}\t0.0\t0.0\t0.0\t{:7.3f}\t10.0e-6\t977.0\t{}'.format(
                    0.015, y, z, ti/1000.0, pId))
                pId += 1
    

            
if __name__ == '__main__':
    BreathParticleGeneration()
