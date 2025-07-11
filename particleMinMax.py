#!/home/wayne/bin/pyenv/bin/python3

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

def particleMinMax():
    if DEBUG:
        import pdb
        pdb.set_trace()

    args = myparser()

    xmin = MinMaxList("x", "min")
    xmax = MinMaxList("x", "max")
    ymin = MinMaxList("y", "min")
    ymax = MinMaxList("y", "max")
    zmin = MinMaxList("z", "min")
    zmax = MinMaxList("z", "max")

    for fn in args.filename:
        print('Reading file {}'.format(fn), flush=True)
              
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

                fields[PID] = int(fields[PID])
                fields[X] = float(fields[X])
                fields[Y] = float(fields[Y])
                fields[Z] = float(fields[Z])

                xmin.add(fields[PID], fields[X])
                # xmin.check()
                xmax.add(fields[PID], fields[X])
                # xmax.check()

                ymin.add(fields[PID], fields[Y])
                # ymin.check()
                ymax.add(fields[PID], fields[Y])
                # ymax.check()

                zmin.add(fields[PID], fields[Z])
                # zmin.check()
                zmax.add(fields[PID], fields[Z])
                # zmax.check()

    xmin.print()
    xmax.print()
    ymin.print()
    ymax.print()
    zmin.print()
    zmax.print()

'''
    with open('Part.Traj.Min.Max.json', 'w') as f:
        json = jsonpickle.encode(xmin)
        f.write('{}\n'.format(json))
        json = jsonpickle.encode(xmax)
        f.write('{}\n'.format(json))
        json = jsonpickle.encode(ymin)
        f.write('{}\n'.format(json))
        json = jsonpickle.encode(ymax)
        f.write('{}\n'.format(json))
        json = jsonpickle.encode(zmin)
        f.write('{}\n'.format(json))
        json = jsonpickle.encode(zmax)
        f.write('{}\n'.format(json))
'''


class MinMaxList(object):
    def __init__(self, coord, type):
        self.coord = coord
        self.type = type
        self.size = 0
        self.maxSize = 10
        self.cnt = 0
        self.coor = []
        self.pId = []

    def add(self, pId, x):
        # self.cnt += 1
        
        if self.size == 0:
            self.pId.append(pId)
            self.coor.append(x)
            self.size += 1
            return;

        if self.type == 'min':
            insLoc = self.placeMin(x)
            if insLoc == -1:
                return
        
        if self.type == 'max':
            insLoc = self.placeMax(x)
            if insLoc == -1:
                return

        # See if same pId exists "above" this entry.  Above means negative
        # index direction.  If found, continue without adding to list.
        for i in range(insLoc-1, -1, -1):
            if pId == self.pId[i]:
                return

        # Insert into lists
        self.pId.insert(insLoc, pId)
        self.coor.insert(insLoc, x)
        self.size += 1
        
        # Look from max index to insertion index to see if pId is in the
        # range.  If so, remove from list.
        for i in range(self.size-1, insLoc, -1):
            if (self.pId[i] == pId):
                self.pId.pop(i)
                self.coor.pop(i)
                self.size -= 1
                break

        # Make sure lists don't exceed maximum size.
        while self.size > self.maxSize:
            self.pId.pop()
            self.coor.pop()
            self.size -= 1

            
    def placeMin(self, x):
        if x < self.coor[0]:
            return 0

        if x > self.coor[self.size-1]:
            if self.size == self.maxSize:
                return -1
            else:
                return self.size

        # Find location where it should be inserted.
        for i in range(0, self.size):
            if x < self.coor[i]:
                return i
            
        return -1

        
    def placeMax(self, x):
        if x > self.coor[0]:
            return 0

        if x < self.coor[self.size-1]:
            if self.size == self.maxSize:
                return -1
            else:
                return self.size

        # Find location where it should be inserted.
        for i in range(0, self.size):
            if x > self.coor[i]:
                return i
            
        return -1

        
    def print(self):
        print("\n{} {}".format(self.coord, self.type))
        for i in range (self.size):
            print("  {:5d}\t{:.15e}".format(self.pId[i], self.coor[i]))

    # See if table is in order.
    def check(self):
        if self.type == 'min':
            for i in range(0, len(self.pId)-1):
                if self.coor[i] <= self.coor[i+1]:
                    continue
                break
            else:
                return
        else:
            for i in range(0, len(self.pId)-1):
                if self.coor[i] >= self.coor[i+1]:
                    continue
                break
            else:
                return

        print('\narray out of order on cnt {}\n'.format(self.cnt))
        self.print()
        
        
def myparser():
    mainparser = argparse.ArgumentParser()

    mainparser.add_argument('filename', nargs='+',
                            help='List of file names to be concated')

    arg = mainparser.parse_args()
    return (arg)


if __name__ == '__main__':
    particleMinMax()


# Variables =     X,      Y,      Z,      U,      V,      W,      STIME,  ETIME,  DIA,    DEN,    OUT,    TIME,   PID,    BID
# Zone I = 1,     J = 5,  DATAPACKING = POINT
# 1.500000e-02    9.830118e-01    1.677469e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    0       1
# 1.500000e-02    9.818900e-01    1.676593e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    1       1
# 1.500000e-02    1.011689e+00    1.679744e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    2       1
# 1.500000e-02    9.987945e-01    1.670595e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    3       1
# 1.500000e-02    9.913230e-01    1.676918e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    1.000000e-05    9.770000e+02    0       0.000000e+00    4       1
