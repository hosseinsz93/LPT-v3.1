#!/home/wayne/bin/pyenv/bin/python3

# perry
###!/home/wayne/bin/pyenv/bin/python3

# seawulf
# module load python/3.8.6
###!/gpfs/software/python/3.8.6/bin/python3

# zagros
###!/mnt/lustre/wayne/bin/python3/bin/python3

import argparse

DEBUG = False
# DEBUG = True


def partTraj():
    if DEBUG:
        import pdb
        pdb.set_trace()

    # The first record in the first file will be copied as the header.
    first = False

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
    listp = []
    listfd = [ ]
    for p in pIds:
        listp.append(p)
        listfd.append(open('Part.Traj.p{:06d}.dat'.format(p), 'w'))
    
    for fn in args.filename:
        print('Reading file {}'.format(fn), flush=True)
              
        with open(fn, 'r') as f:
            # Handle header record
            if not first or args.fheader == 1:
                buff = f.readline()
                if not buff:
                    break
                if not first:
                    first = True
                    outputBuffs(listfd, buff)

            while True:
                buff = f.readline()
                if not buff:
                    break

                # Try find to see if I can speed up the search.
                # 1 is pid, 12 is stat
                i1 = buff.find('\t')
                i2 = buff.find('\t', i1+1)
                field1 = int(buff[i1+1:i2])
                try:
                    i3 = listp.index(field1)
                except ValueError as ve:
                    continue

                i4 = buff.rfind('\t')
                field12 = int(buff[i4+1:len(buff)-1])
                # debugging code
                if field1 ==  1 or \
                   field1 ==  2 or \
                   field1 ==  4 or \
                   field1 ==  6 or \
                   field1 ==  7 or \
                   field1 == 10 or \
                   field1 == 20:
                    listfd[i3].write(buff)
                elif field12 == 0:
                    listfd[i3].write(buff)

                # regular code
                # if field12 == 0:
                #     listfd[i3].write(buff)


    # Close output files.
    for fd in listfd:
        fd.close()

        
def outputBuffs(listfd, buff):
    for fd in listfd:
        fd.write(buff)

    
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

    arg = mainparser.parse_args()
    return (arg)


if __name__ == '__main__':
    partTraj()


# time    pId     coor.x  coor.y  coor.z  vel.x   vel.y   vel.z   stime   etime   dia     den     partStat
# 3.175031e+01    0       1.002787e+00    9.073364e-01    1.759363e+00    1.488854e-02    -7.370685e-04   5.955675e-04    0.000000e+00    3.175031e+01    1.000000e-07    9.770000e+02    0
