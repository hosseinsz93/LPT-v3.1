#!/home/wayne/bin/pyenv/bin/python3

# perry
###!/home/wayne/bin/pyenv/bin/python3

# seawulf
# module load python/3.8.6
###!/gpfs/software/python/3.8.6/bin/python3

# zagros
###!/mnt/lustre/wayne/bin/python3/bin/python3


import argparse
import re


DEBUG = False
# DEBUG = True


def particleBlanking():
    if DEBUG:
        import pdb
        pdb.set_trace()

    args = myparser()

    reForY = re.compile(r'[^\s]+\s+([^\s]+)\s')
    
    for fni in args.filename:
        print(f'Reading file {fni}', flush=True)

        # Rename filename to filename.save which becomes the input file
        # then open the output file with filename.
        # Particle043850_0.dat --> Particle043850_0_side_blanked.dat
        # fni = fno + '.save'

        recs = { }
        recCnt = -1
        with open(fni, 'r') as fi:
            for rec in fi:
                recCnt += 1

                if recCnt == 0:
                    header1 = rec

                elif recCnt == 1:
                    header2 = rec

                else:
                    m = reForY.match(rec)
                    y = float(m.group(1))
                    if y > 1.025 or y < 0.975:
                        continue
                    recs[recCnt-2] = rec

        fno = fni[:fni.rfind('.dat')] + '_side_blanked.dat'
        with open(fno, 'w') as fo:
            fo.write(header1)
            m = re.match(
         r'(Zone\s+I\s+=\s+1,\s+J\s+=\s+)[0-9]+(,\s+DATAPACKING\s+=\s+POINT\n)',
                header2)
            fo.write(f'{m.group(1)}{len(recs)}{m.group(2)}')
            for rec in sorted(recs.keys()):
                fo.write(recs[rec])
            
            debug = 1


def myparser():
    mainparser = argparse.ArgumentParser()

    mainparser.add_argument('filename', nargs='+',
                            help='List of file names to be concated')

    arg = mainparser.parse_args()
    return (arg)


if __name__ == '__main__':
    particleBlanking()


'''
Variables =     X,      Y,      Z,      U,      V,      W,      STIME,  ETIME,  DIA,    DEN,    OUT,    TIME,\
   PID,    BID
Zone I = 1,     J = 105226,     DATAPACKING = POINT
2.176427e+00    7.284295e-01    1.772439e+00    3.282556e-03    -5.393515e-04   -1.149734e-03   0.000000e+00 \
   2.192500e+02    1.000000e-07    9.770000e+02    0       2.192500e+02    0       1
'''
