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

DEBUG = False
DEBUG = True


def linkFilesWithDequName():
    if DEBUG:
        import pdb
        pdb.set_trace()

    args = myparser()

    cnt = 0
    
    for fno in args.filename:
        print('linking file {}'.format(fno), flush=True)

        dot = fno.rfind('.')
        os.symlink(fno, f'ffmpeg{cnt:05d}{fno[dot:]}')
        cnt += 1


def myparser():
    mainparser = argparse.ArgumentParser()

    mainparser.add_argument('filename', nargs='+',
                            help='List of file names to be concated')

    arg = mainparser.parse_args()
    return (arg)


if __name__ == '__main__':
    linkFilesWithDequName()
