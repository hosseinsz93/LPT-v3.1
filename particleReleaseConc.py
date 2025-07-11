#!/mnt/lustre/wayne/bin/python3/bin/python3

# perry
###!/home/wayne/bin/pyenv/bin/python3

# seawulf
# module load python/3.8.6
###!/gpfs/software/python/3.8.6/bin/python3

# zagros
###!/mnt/lustre/wayne/bin/python3/bin/python3

# find the 10 min and max coordinates for x, y and z.


# Finds min/max x, y and z values for Particle999999_0.dat files.


import argparse
# import jsonpickle
import re

DEBUG = False
DEBUG = True

X =     0
Y =     1
Z =     2
U =     3
V =     4
W =     5
STIME = 6
ETIME = 7
DIA =   8
DEN =   9
OUT =  10
TIME = 11
PID =  12
BID =  13

def particleReleaseConc():
    if DEBUG:
        import pdb
        pdb.set_trace()

    args = myparser()

    origin = [ -7.7418, 3.6865, 1.0000 ]
    rad = 0.16

    print('iter\tcnt')
    
    for fn in args.filename:
        cnt = 0
        m = re.search('.+([0-9]{6}).+$', fileIter)
        
        with open(fn, 'r') as f:
            for rec in f:
                # First character must be '+-[0-9].'
                if not re.search('^[0-9+-.]', rec):
                    continue

                fields = rec.split()
                fields[OUT] = int(fields[OUT])
                
                # Ignore part stat != 0
                if fields[OUT] != 0:
                    continue

                fields[X] = float(fields[X])
                fields[Y] = float(fields[Y])
                fields[Z] = float(fields[Z])

                if (fields[X]*fields[X] + fields[Y]*fields[Y]
                    + fields[Z]*fields[Z] <= rad*rad):
                    cnt += 1

        print('{}\t{}\n'.format(m.group(1), cnt))
        
        
def myparser():
    mainparser = argparse.ArgumentParser()

    mainparser.add_argument('filename', nargs='+',
                            help='List of file names to be concated')

    arg = mainparser.parse_args()
    return (arg)


if __name__ == '__main__':
    particleReleaseConc()


# Variables =     X,      Y,      Z,      U,      V,      W,      STIME,  ETIME,  DIA,    DEN,    OUT,    TIME,   PID,    BID
# Zone I = 1,     J = 5,  DATAPACKING = POINT
# 1.500000e-02    9.830118e-01    1.677469e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    0       1
# 1.500000e-02    9.818900e-01    1.676593e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    1       1
# 1.500000e-02    1.011689e+00    1.679744e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    2       1
# 1.500000e-02    9.987945e-01    1.670595e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    3       1
# 1.500000e-02    9.913230e-01    1.676918e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    4       1
