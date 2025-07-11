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
# partdebug = True



def Output(ti, pId, coord):
    print('{:.15f}\t{:.15f}\t{:.15f}\t0.0\t0.0\t0.0\t{:7.3f}\t10.0e-6\t977.0\t{}'.format(
           coord[0], coord[1], coord[2], ti/1000.0, pId))
    return pId + 1


def OutputStar(ti, pId, coord):
    print('v, {}, {:.15f}, {:.15f}, {:.15f}'.format(pId+1, coord[0], coord[1], coord[2]))
    return pId + 1

whichOutput = Output
# whichOutput = OutputStar



def NoseParticleGeneration():
    if debug:
        import pdb;
        pdb.set_trace()

    '''
    Generate tracking particles for the mouth and nose.

    mouth:
    One breath cycle: 2.5 exhale
                      2.5 inhale
                      5.0 TOTAL

    Version 2: output 5 particles every 5ms
    Each time step is 5ms.  There are currently 20000 time steps.
    I expect more iterations, so I'me going to go to 30000 time steps.
    That will be 150s.
    For the exhale step, there will be 2.5 / 0.005 = 500 particle released cycles
    release from 0 to 360 time steps for 1.8s
    loop from 0, 0.005, 150.005
    multiply by 1000 and use integer arithmatic: 0, 5, 150005
    The cycle is 5s or 5000 which us used to mod time to compare when
    particles need to be generated.

    For the with nose simulation, I changed the mesh so that the mouth is
    centered around y = 0. The previous simulation centered mount around y = 1.
    '''

    ymin = 0.98 - 1.00
    ymax = 1.02 - 1.00
    ydel = ymax - ymin
    zmin = 1.67
    zmax = 1.68
    zdel = zmax - zmin

    '''
    nostrils:
    Generations circle, radius: 3.75 mm

    rotation | translation maxtrix
    left nostril   [ [ 0.948, 0.000, -0.319, |  0.007   ],
                     [ 0.000, 1.000,  0.000, | -0.00875 ],
                     [ 0.319, 0.000,  0.948, |  1.7122  ] ]

    right nostril  [ [ 0.948, 0.000, -0.319, |  0.007   ],
                     [ 0.000, 1.000,  0.000, |  0.00875 ],
                     [ 0.319, 0.000,  0.948, |  1.7122  ] ]
    '''

    # Nose parameters
    radius = 0.00375

    leftRotTrans  = [ [ 0.948, 0.000, -0.319,  0.007   ],
                      [ 0.000, 1.000,  0.000, -0.00875 ],
                      [ 0.319, 0.000,  0.948,  1.7122  ] ]

    rightRotTrans = [ [ 0.948, 0.000, -0.319,  0.007   ],
                      [ 0.000, 1.000,  0.000,  0.00875 ],
                      [ 0.319, 0.000,  0.948,  1.7122  ] ]


    pId = 0
    if whichOutput == Output:
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
                pId = whichOutput(ti, pId, [ 0.015, y, z ])

            # Generate nostril particles
            pId = GenNostril(ti, radius, pId, leftRotTrans)
            pId = GenNostril(ti, radius, pId, rightRotTrans)
            

def GenNostril(ti, radius, pId, rotTrans):
    # Two particles per nostril per exhale cycle
    x = 0.0
    y = 0.0
    
    for i in range(0,2):
        while True:
            # These calculations gave the best particle distribution throughout
            # the model volume.
            x = radius * random.uniform(-1.0, 1.0)
            y = radius * random.uniform(-1.0, 1.0)
            if x*x + y*y < radius*radius:
                break

        coord = RotTrans(rotTrans, [ x, y, 0.0 ])
        pId = whichOutput(ti, pId, coord)
            
    return pId


def RotTrans(rotTrans, coord):
    return [ rotTrans[0][0] * coord[0] + rotTrans[0][1] * coord[1] + rotTrans[0][2] * coord[2] + rotTrans[0][3],
             rotTrans[1][0] * coord[0] + rotTrans[1][1] * coord[1] + rotTrans[1][2] * coord[2] + rotTrans[1][3],
             rotTrans[2][0] * coord[0] + rotTrans[2][1] * coord[1] + rotTrans[2][2] * coord[2] + rotTrans[2][3] ]

            
if __name__ == '__main__':
    NoseParticleGeneration()
