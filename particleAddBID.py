#!/home/wayne/bin/pyenv/bin/python3

# perry
###!/home/wayne/bin/pyenv/bin/python3

# seawulf
# module load python/3.8.6
###!/gpfs/software/python/3.8.6/bin/python3

# zagros
###!/mnt/lustre/wayne/bin/python3/bin/python3


# Add breath id to existing Particle999999_0.dat files.


import argparse
import os
import re

DEBUG = False
# DEBUG = True


def particleAddBID():
    if DEBUG:
        import pdb
        pdb.set_trace()

    # The first record in the first file will be copied as the header.
    first = False

    args = myparser()
    
    for fno in args.filename:
        print('Reading file {}'.format(fno), flush=True)

        # Rename filename to filename.save which becomes the input file
        # then open the output file with filename.
        fni = fno + '.save'
        os.rename(fno, fni)
        os.system('chmod 400 {}\n'.format(fni))
              
        with open(fni, 'r') as fi:
            with open(fno, 'w') as fo:

                recCnt = 0

                for buf in fi:
                    recCnt += 1

                    if recCnt == 1:
                        fo.write(buf[0:len(buf)-1] + ',\tBID\n')

                    elif recCnt == 2:
                        fo.write(buf)

                    else:
                        se = re.search('^[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t', buf)
                        bId = float(se.group(1)) / 5.0 + 1
                        fo.write(buf[0:len(buf)-1] + '\t{}\n'.format(int(bId)))


def myparser():
    mainparser = argparse.ArgumentParser()

    mainparser.add_argument('filename', nargs='+',
                            help='List of file names to be concated')

    arg = mainparser.parse_args()
    return (arg)


if __name__ == '__main__':
    particleAddBID()
