#!/home/wayne/bin/pyenv/bin/python3

# perry
###!/home/wayne/bin/pyenv/bin/python3

# seawulf
# module load python/3.8.6
###!/gpfs/software/python/3.8.6/bin/python3

# zagros
###!/mnt/lustre/wayne/bin/python3/bin/python3

# find the 10 min and max coordinates for x, y and z.

import argparse
import jsonpickle
import re

DEBUG = False
# DEBUG = True


def partTrajMinMax():
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

                # Ignore part stat != 0
                if int(rec[rec.rfind('\t')+1:len(rec)-1]) != 0:
                    continue

                # Field 1: time = rec[i0:i1-1]
                i0 = 0
                i1 = rec.find('\t', i0)

                # Field 2: pId = rec[i0:i1-1]
                i0 = i1+1
                i1 = rec.find('\t', i0)
                pId = int(rec[i0:i1])

                # Field 3: x = rec[i0:i1-1]
                i0 = i1+1
                i1 = rec.find('\t', i0)
                x = float(rec[i0:i1])
                xmin.add(pId, x)
                # xmin.check()
                xmax.add(pId, x)
                # xmax.check()

                # Field 4: y = rec[i0:i1-1]
                i0 = i1+1
                i1 = rec.find('\t', i0)
                y = float(rec[i0:i1])
                ymin.add(pId, y)
                # ymin.check()
                ymax.add(pId, y)
                # ymax.check()

                # Field 5: z = rec[i0:i1-1]
                i0 = i1+1
                i1 = rec.find('\t', i0)
                z = float(rec[i0:i1])
                zmin.add(pId, z)
                # zmin.check()
                zmax.add(pId, z)
                # zmax.check()

    xmin.print()
    xmax.print()
    ymin.print()
    ymax.print()
    zmin.print()
    zmax.print()

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
    partTrajMinMax()


# time    pId     coor.x  coor.y  coor.z  vel.x   vel.y   vel.z   stime   etime   dia     den     partStat
# 3.175031e+01    0       1.002787e+00    9.073364e-01    1.759363e+00    1.488854e-02    -7.370685e-04   5.955675e-04    0.000000e+00    3.175031e+01    1.000000e-07    9.770000e+02    0
