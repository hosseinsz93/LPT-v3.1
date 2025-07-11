#!/home/wayne/bin/pyenv/bin/python3

# perry
###!/home/wayne/bin/pyenv/bin/python3

# seawulf
# module load python/3.8.6
###!/gpfs/software/python/3.8.6/bin/python3

# zagros
###!/mnt/lustre/wayne/bin/python3/bin/python3

import argparse
import math
import re

DEBUG = False
# DEBUG = True

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

def particleMapReduce():
    if DEBUG:
        import pdb
        pdb.set_trace()

    args = myparser()

    # Setup -p list. Check that they are all integers.  Open outputfiles
    # for each pId.
    err = False
    try:
        pIds = args.part.split(',')
        for idx, p in enumerate(pIds):
            pIds[idx] = int(p)

    except ValueError as ve:
        exp = re.search("'(.+)'", str(ve))
        print('All part parameters must be integers, \'{}\' is not.'.format(exp.group(1)))
        err = True

    if err:
        return 1

    
    # Open output files.
    listfd = { }
    for p in pIds:
        listfd[p] = open('Particle.p{:06d}_0.dat'.format(p), 'w')
        if args.x:
            listfd[p].write('x\ty\tz\tvmag\n')
    
    for fn in args.filename:
        print('Reading file {}'.format(fn), flush=True)
              
        with open(fn, 'r') as f:
            for rec in f:
                # First character must be '+-[0-9].'
                if not re.search('^[0-9+-.]', rec):
                    continue

                fields = rec.split()
                out = int(fields[OUT])
                
                # Ignore part stat != 0
                if out != 0:
                    continue

                pId = int(fields[PID])
                if not pId in listfd:
                    continue

                if not args.x:
                    listfd[pId].write(rec)
                else:
                    u = float(fields[U])
                    v = float(fields[V])
                    w = float(fields[W])
                    vmag = math.sqrt(u*u + v*v + w*w)
                    listfd[pId].write(f'{fields[X]}\t{fields[Y]}\t{fields[Z]}\t{vmag:.6e}\n')

    # Close output files.
    for fd in listfd.values():
        fd.close()

        
def outputRecs(listfd, rec):
    for fd in listfd:
        fd.write(rec)

    
def myparser():
    mainparser = argparse.ArgumentParser()

    mainparser.add_argument('-p', '--part', help='The particles to extact, ' +
                            'a comma separated list with no spaces',
                            required=True)
    mainparser.add_argument('-s', '--skip',
                    help='Number of records to skip between each output record',
                            type=int, default=-1)
    mainparser.add_argument('-f', '--fheader',
                            help='1: all files have header, ' +
                                 '2: only first file has header',
                            type=int, default=1, choices=[1, 2])
    mainparser.add_argument('filename', nargs='+',
                            help='List of file names to be concated')
    # mainparser.add_argument('-v', '--vmag', action='store_true', help='Output format needed to particle traces: x, y, z, vmag')
    mainparser.add_argument('-x', '--x', action='store_true', help='Output format needed to particle traces: x, y, z, vmag')

    arg = mainparser.parse_args()
    return (arg)


if __name__ == '__main__':
    particleMapReduce()


# Variables =     X,      Y,      Z,      U,      V,      W,      STIME,  ETIME,  DIA,    DEN,    OUT,    TIME,   PID,    BID
# Zone I = 1,     J = 5,  DATAPACKING = POINT
# 1.500000e-02    9.830118e-01    1.677469e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    0       1
# 1.500000e-02    9.818900e-01    1.676593e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    1       1
# 1.500000e-02    1.011689e+00    1.679744e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    2       1
# 1.500000e-02    9.987945e-01    1.670595e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    3       1
# 1.500000e-02    9.913230e-01    1.676918e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    4       1
