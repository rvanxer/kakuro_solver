"""written: Wed, 2019-11-06 14:53:38, author: Rajat Verma

Kakuro puzzle solver: This program first uses multiple simple tricks to shrink 
the possible "pool" of values each cell can attain. If it cannot
get the result then, it uses a pseudo-brute force approach to try as many 
feasible combinations as possible. Correct solution is not always guaranteed.

Usage: either instantiate a `Problem` by providing the path of the CSV file
of the puzzle (see examples for the format) or uncomment last line to try"""

import logging
from time import time
from itertools import groupby, combinations, product
from operator import itemgetter
import csv
import numpy as np
import matplotlib.pyplot as plt

# *****************************************************************************
class Cell(object):
    def __init__(self, i, j, digits=(0)):
        self.i, self.j = i, j
        self.loc = (i, j)
        self.val = 0
        self.final = False
        self.pool = set(digits)
        self.vects = {}

    def __repr__(self):
        return '<%d,%d>' % (self.i+1, self.j+1)

    def deterministic_solve(self):
        """Solve self and its neighbors recursively if its pool size is 1"""
        if len(self.pool) == 1 and not self.final:
            self.val = list(self.pool)[0]
            self.final = True
            for v in self.vects.values():
                v.total -= self.val
                for c in v.cells():
                    c.pool.discard(self.val)

# *****************************************************************************
class Vector(object):
    def __init__(self, id, cells, type, total):
        self.id = id # not a unique identifier but a serial num
        self.type = type # type is either 'row' or 'col'
        self.origCells = cells
        self.origTotal = total # original total
        self.total = total # effective total
        for c in self.cells():
            c.vects[self.type] = self

    def __repr__(self):
        return '%s%d<%d,%d>' % (self.type[0].upper(), self.id,
            self.origCells[0].i + 1, self.origCells[0].j + 1)

    def cells(self): # effective cells, i.e., those that are not final
        return [c for c in self.origCells if not c.final]

    def length(self): # effective length
        return len(self.cells())

    def feasible_combs(self):
        """Get the list of all feasible combinations of self"""
        return [f for f in product(*(c.pool for c in self.cells())) if
                sum(f) == self.total and len(set(f)) == len(f)]

    def filter_poolmap(self, poolmap):
        """Filter the pools of self k cells using the poolmap of the problem"""
        if self.length() > 0:
            feasiblePool = poolmap[self.length()][self.total]
            for c in self.cells():
                c.pool = c.pool.intersection(feasiblePool)

    def solve_tup1(self):
        """For size-1 vector, solve it deterministically"""
        if self.length() == 1:
            self.cells()[0].pool = {self.total}
            self.cells()[0].deterministic_solve()

    def complement_tup2(self):
        """Filter pools of size-2 vectors so that they complement each other"""
        if self.length() == 2:
            for i, c in enumerate(self.cells()):
                complement = {self.total - j for j in self.cells()[i-1].pool}
                c.pool = c.pool.intersection(complement)

    def min_max_check(self):
        """Check for each cell if its pool values satisfy the vector sum 
        within bounds (very loose check)"""
        for cell in self.cells():
            minVals = [min(c.pool) for c in self.cells() if c != cell]
            maxVals = [max(c.pool) for c in self.cells() if c != cell]
            for val in list(cell.pool):
                remainder = self.total - val
                if remainder < sum(minVals) or remainder > sum(maxVals):
                    cell.pool.discard(val)

    def paired_uncertainty(self):
        """If two of self k cells have size-2 pool, then remove these values
        from pools of other self k cells"""
        if self.length() > 2:
            pools = [c.pool for c in self.cells()]
            counts = [pools.count(p) for p in pools]
            idx = [i for i, x in enumerate(counts) if x==2] # ids for 2-pools
            if len(idx) == 2:
                p = pools[idx[0]] # pool values to be discarded
                if len(p) == 2:
                    for i, c in enumerate(self.cells()):
                        if i not in idx:
                            for val in p:
                                c.pool.discard(val)

# *****************************************************************************
class Problem(object):
    def __init__(self, file, digits=list(range(1,10)), nTryVect=-1, 
                 logFile='logger.log'):
        self.time = {'start': time()}
        self.file = file
        self.digits = digits
        self.logFile = logFile
        self.solved = False
        # pre-processing
        self.configure_logging()
        self.set_poolmap()
        self.import_data()
        # main solving
        self.solve(firstTime=True)
        if not self.solved:
            self.try_solve(nTryVect)
            
    def __repr__(self):
        return 'Problem<"%s">' % self.file.split('/')[-1]
            
    def configure_logging(self):
        logging.basicConfig( filename=self.logFile, level=logging.INFO,
            filemode='w', format='%(levelname)s: %(message)s' )

    def set_poolmap(self):
        """Given a set of allowable digits, store in a dict all the possibile 
        unordered sets of values a vector of size n & total s can achieve"""
        A = {} # dict of allowable pools
        for n in self.digits:
            A[n] = {}
            for comb in combinations(self.digits, n):
                s = sum(comb)
                if len(set(comb)) == n:
                    if s not in A[n].keys():
                        A[n][s] = set(comb)
                    else:
                        A[n][s] = A[n][s].union(set(comb))
        self.poolmap = A

    def import_data(self):
        """Read the problem CSV into a matrix, cell list & vector list"""
        self.mat = []
        self.cells = []
        self.rows, self.cols = [], []
        with open(self.file) as f:
            for i, line in enumerate(csv.reader(f)):
                row = []
                for j, val in enumerate(line):
                    if val == '':
                        c = Cell(i, j, self.digits)
                        row.append(c)
                        self.cells.append(c)
                    elif val == 'x':
                        row.append(None)
                    elif '\\' in val:
                        left, right = val.split('\\')
                        left = int(left) if left != '' else 0
                        right = int(right) if right != '' else 0
                        row.append((left, right, val))
                self.mat.append(row)
        for i, row in enumerate(self.mat):
            colIdx = [c.j for c in row if isinstance(c, Cell)]
            for k, g in groupby(enumerate(colIdx), lambda x: x[0]-x[1]):
                J = list(map(itemgetter(1), g)) # column IDs
                cells = [self.mat[i][j] for j in J]
                rowSum = self.mat[i][J[0]-1][1]
                rowObj = Vector(len(self.rows), cells, 'row', rowSum)
                self.rows.append(rowObj)
        for j, col in enumerate(map(list, zip(*self.mat))):
            rowIdx = [c.i for c in col if isinstance(c, Cell)]
            for k, g in groupby(enumerate(rowIdx), lambda x: x[0]-x[1]):
                I = list(map(itemgetter(1), g)) # row IDs
                cells = [self.mat[i][j] for i in I]
                colSum = self.mat[I[0]-1][j][0]
                colObj = Vector(len(self.cols), cells, 'col', colSum)
                self.cols.append(colObj)
        self.vects = self.rows + self.cols
        self.time['import_data'] = time() - self.time['start']

    def save_state(self):
        """Save current problem state (set of values of cells and vects)"""
        self.state = {'vects': [(v.cells, v.total) for v in self.vects],
            'cells': [(c.final, c.val, c.pool.copy()) for c in self.cells]}

    def reset_state(self):
        """Restore the problem state to the state after the first solve()"""
        for i in range(len(self.cells)):
            final, val, pool = self.state['cells'][i]
            self.cells[i].final = final
            self.cells[i].val = val
            self.cells[i].pool = pool
        for i in range(len(self.vects)):
            cells, total = self.state['vects'][i]
            self.vects[i].cells = cells
            self.vects[i].total = total

    def print_state(self):
        """Print the current problem state as a matrix (heatmap chart)"""
        A = self.mat
        result = np.zeros((len(A), len(A)), np.int8)
        fig, ax = plt.subplots(figsize=(8,8))
        params = {'ha': 'center', 'va': 'center'}
        for i, row in enumerate(A):
            for j, elem in enumerate(row):
                if elem is None:
                    result[i, j] = -1
                elif isinstance(elem, Cell):
                    result[i, j] = elem.val
                    if elem.final:
                        ax.text(j, i, result[i,j], size=16, **params)
                    else:
                        ax.text(j, i, size=10, **params, s=''.join(
                                sorted([str(x) for x in elem.pool])))
                elif isinstance(elem, tuple):
                    result[i, j] = -1
                    ax.text(j, i, elem[2], size=12, color='grey', **params)
        ax.imshow(np.sign(result), cmap='Pastel1')
        ax.xaxis.set_ticks_position('top')
        positions, labels = range(len(A)), range(1, len(A)+1)
        plt.xticks(positions, labels, size=12)
        plt.yticks(positions, labels, size=12)
        plt.show()
        plt.close()
        self.time['total'] = time() - self.time['start']
        resultLog = '%s (%d cells)' % (self, len(self.cells)) + \
            '\nTotal time elapsed: %.4fs' % self.time['total']
        print(resultLog)
        logging.info(resultLog + '\n' + '*'*50)

    def check_solution(self):
        """Check if the current state is the correct (unique) solution"""
        for v in self.vects:
            if v.total != sum([c.val for c in v.cells()]):
                return False
        for c in self.cells:
            if c.val == 0:
                return False
        return True

    def solve(self, numIter=3, firstTime=False):
        """Main solution step for shrinking the pools and solving as many
        cells as possible"""
        for v in self.vects:
            v.filter_poolmap(self.poolmap)
        for _ in range(numIter):
            for v in self.vects:  v.complement_tup2()
            for v in self.vects:  v.min_max_check()
            for v in self.vects:  v.paired_uncertainty()
            for v in self.vects:  v.solve_tup1()
            for c in self.cells:  c.deterministic_solve()
        if firstTime:
            self.time['first_solve'] = time() - self.time['start']
            self.save_state()
        self.solved = self.check_solution()
        if self.solved:
            self.time['solve'] = time() - self.time['start']
            self.print_state()

    def try_solve(self, nTryVect):
        """Try solving the problem by trying all feasible combinations on 
        at most a given number of high-combination vectors"""
        fcSizes = [len(v.feasible_combs()) for v in self.vects]
        priorV = sorted(enumerate(fcSizes), key=lambda x: x[1], reverse=True)
        for vi, _ in priorV[:nTryVect]:
            v = self.vects[vi]
            fc = v.feasible_combs()
            print('Trying %s (#fc=%d)' % (v, len(fc)))
            for f in fc:
                self.reset_state()
                for i, c in enumerate(v.cells()):
                    c.pool = {f[i]}
                try:
                    self.solve()
                    if self.solved:
                        return
                except ValueError:
                    logging.warning('f:%s: value error on %s' % (str(f), v))
            self.reset_state()
        self.time['try_solve'] = time() - self.time['start']
        print('Solution not found :( *_*')
        print('Time elapsed: %.4f' % self.time['try_solve'])

# Example problem (uncomment and run file to see) *****************************
#p = Problem('Problems/hard2.csv')
