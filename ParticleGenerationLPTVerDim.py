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
# debug = True

version = '20210706'

def ParticleGenerationLPTVerDim():
    if debug:
        import pdb;
        pdb.set_trace()

    print('Running ParticleGenerationLPTVerDim, v{}'.format(version))

# test code
#    for i, t in enumerate(range(0, 525, 25)):
#        locPrint(pymajor, '{}\t{}\t{}'.format(i, t/100, math.sin(2*math.pi*t/500)))
    x    = 6.667e-04
    ymin = 0.195
    ymax = 0.262
    ydel = ymax - ymin
    zmin = 0.0267
    zmax = 0.0667
    zdel = zmax - zmin
    dia = 6.82E-05
    den = 8800.0
    timeOffset = 3.0

# Time step: 0.0001s
# Total release time: 0.07s
# 70 particles per time step

    pId = 0
    with open('ParticleInitialLPTVerDim.dat', 'w') as fo:
        fo.write('# location  velocity  "start time"  diameter  density pId\n')
        fo.write('#  x y z     u v w\tgenerate by ParticleGenerationLPTVerDim, v{}\n'
              .format(version))
        for time in range(0, 700, 1):
            for piter in range(70):
                rydel = random.uniform(0.0, ydel)
                y = ymin + rydel
                rzdel = random.uniform(0.0, zdel)
                z = zmin + rzdel

                fo.write('{:.15f}\t{:.15f}\t{:.15f}\t0.0\t0.0\t0.0\t{:.4f}\t{:.15f}\t{:.15f}\t{}\n'.format(
                    x, y, z, time/10000.0 + timeOffset, dia, den, pId))
                pId += 1
    

            
if __name__ == '__main__':
    ParticleGenerationLPTVerDim()
