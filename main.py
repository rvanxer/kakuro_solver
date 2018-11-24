import csv
import sympy
import itertools
import time
import numpy as np

def importData(fileName):
	mat = []
	with open(fileName) as file:
		reader = csv.reader(file)
		for row in reader:
			row2 = []
			for elem in row:
				if elem == '':
					elem = 0
				if elem == 'x':
					elem = -1
				try:
					elem = int(elem)
				except ValueError:
					elem = elem
				if isinstance(elem,str) and '\\' in elem:
					elem = elem.split('\\')
					for i,e in enumerate(elem):
						if e == '':
							elem[i] = -1
						else:
							elem[i] = int(e)
				row2.append(elem)
			mat.append(row2)
	return mat

def print_mat(mat):
	"""Displays the current state of the main problem matrix, i.e. trial cell values"""
	# print('***********************')
	for r in mat: # r = row
		R = []
		for c in r: # c = cell
			if c.editable:
				R.append('{:2}'.format(int(c.trialVal)))
			else:
				R.append('{:2}'.format(' .'))
		R = ' '.join(R)
		print(R)
	print('*****************************')

# ****************** CELL class ******************
class Cell(object):
	def __init__(self,mat,i,j):
		# self.data = mat
		self.i = i
		self.j = j
		self.defVal = mat[i][j] if mat[i][j] != -1 else None
		self.val = 0 # actual symbolic/numerical value
		self.blank = True # if self.val == 0
		self.trialVal = 0
		self.editable = True if self.defVal == 0 else False
		self.sumBlock = True if isinstance(mat[i][j],list) else False

	def __repr__(self):
		return 'cell(%d,%d)' % (self.i+1,self.j+1)

	def make_row(self,mat):
		"""Returns the constituent row of an editable cell of a given matrix"""
		i, j = self.i, self.j
		if self.editable:
			start, end, sr = -1, -1, 0 # sr = row sum
			while True:
				j -= 1 # move 1 cell left
				# if cell is at left end of matrix OR left neighbour (j--) is uneditable
				if j == -1 or not A[i][j].editable:
					start = j+1
					sr = mat[i][j][1]
					break
			while True:
				j += 1 # move 1 cell right
				# if cell is at right end of matrix OR right neighbour is uneditable
				if j == len(mat[0]) or not A[i][j].editable: 
					end = j-1
					break
			row = Vector([(i,start),(i,end)])
			row.Sum = sr
			return row
		return None

	def make_col(self,mat):
		"""Returns the editable column of an editable cell of a given matrix"""
		i, j = self.i, self.j
		if self.editable:
			start, end, sc = -1, -1, 0 # sc = column sum
			while True:
				i -= 1 # move 1 cell up
				if i == -1 or not A[i][j].editable:
					start = i+1
					sc = mat[i][j][0]
					break
			while True:
				i += 1 # move 1 cell down
				if i == len(mat) or not A[i][j].editable:
					end = i-1
					break
			col = Vector([(start,j),(end,j)])
			col.Sum = sc
			if col not in Cols:
				Cols.append(col)
			return col
		return None

	def neighbour(self,direction):
		try:
			if 		direction == 'left': 	X = A[self.i][self.j-1]
			elif 	direction == 'right': 	X = A[self.i][self.j+1]
			elif 	direction == 'up':		X = A[self.i-1][self.j]
			elif 	direction == 'down':	X = A[self.i+1][self.j]
			return X
		except Exception as e: 
			# print(e)
			return None

	def populate(self,nx):
		"""Main function for filling a cell with symbolic expressions (1st step of algo)"""
		directions = ['right','down','left','up']
		if self.editable:
			if self.blank:
				if all([not self.calc_expr(d) for d in directions]):
					x = sympy.Symbol('x%d' % (nx+1)) # create a new variable
					self.val = x # store this first sym var as a new variable in the cell val
					Vars[x] = {'bounds':[]}
					nx += 1 # nx = no. of assigned vars from sym var list
					self.blank = False
			for direction in directions:
				child = self.neighbour(direction)
				if child:
					child.calc_expr(direction)
		return nx

	def calc_expr(self,direction):
		"""Calculate the cell value (symbolic expression) in terms of row/col member variables
		If value cannot be calculated right now, return False, else True"""
		if self.editable and self.blank:
			if direction in ['left','right']:
				rowCells = self.row.cells[:]
				rowCells.remove(self)				
				if rowCells and any([x.val == 0 for x in rowCells]): # if any other row member is blank
					return False
				elif self.blank:
					self.val = self.row.Sum - sum([x.val for x in rowCells]) # applying the equation rule
					self.blank = False
					return True
			elif direction in ['up','down']:
				colCells = self.col.cells[:]
				colCells.remove(self)
				if colCells and any([x.val == 0 for x in colCells]):
					return False
				elif self.blank:
					self.val = self.col.Sum - sum([x.val for x in colCells])
					self.blank = False
					return True
		return False

	def get_bounds(self):
		"""Returns the limits of possible values of the vars in cell value"""
		var = self.val.free_symbols # returns the sym vars contained in cell value
		if len(var) == 1: # if the cell contains only one sym var
			lim1 = sympy.solve(self.val - digits[0],var,simplify=False,rational=False)[0]
			lim2 = sympy.solve(digits[-1] - self.val,var,simplify=False,rational=False)[0]
			bounds = [min([lim1,lim2]),max([lim1,lim2])]
			# print('%s: |%s|: %s =' % (self,self.val,list(var)[0]),bounds)
			self.var = list(var)[0]
			self.bounds = bounds
			return self.var, bounds
		return False, False

	def get_coeffs(self):
		"""Create a list of var coefficients from cell value expression"""
		if self.editable:
			coeffs = [0]*(nx+1) # blank coeff array (default 0)
			vars_ = list(self.val.free_symbols)
			k0 = int(self.val.evalf(subs={var: 0 for var in vars_}))
			coeffs[0] = k0
			for var in vars_:
				i = list(Vars.keys()).index(var) # index of current var in the varList
				k = int(self.val.coeff(var)) # coeff of current var in cell value
				coeffs[i+1] = k # i+1 because 1st coeff is constant
			coeffs = np.array(coeffs)
			self.coeffs = coeffs

	def try_sol(self,trialVals):
		"""Put values obtained from 'trial' dict or 'vals' list into cell value expression"""
		if self.editable:
			# vars_ = self.val.free_symbols # set of sym vars in the cell value expr
		# Method 1: evaluating with normal sympy
			# subs = {var: trial[var] for var in vars_}
			# val = int(self.val.evalf(subs=subs)) # num value after putting trial var values
		# Method 2: evaluating with numpy & lambdify
			# f = sympy.lambdify(tuple(vars_),self.val,'numpy')
			# val = f(*tuple([trial[var] for var in vars_]))
		# Method 3: using coeffs array
			K = self.coeffs[1:]
			k0 = self.coeffs[0] # the constant
			# val = sum([a*b for a,b in zip(K,trialVals)]+[k0])
			val = int(np.dot(K,trialVals)) + k0
			self.trialVal = val

	def check_sol(self):
		"""Chekcing if the current cell trial val is the right solution.
		Using design constraints (inequalities & inequations) to discard faulty results"""
		if self.editable:
			rowVals = [cell.trialVal for cell in self.row.cells] # trial vals of row members
			colVals = [cell.trialVal for cell in self.col.cells]
			if any([val==0 for val in rowVals+colVals]): # if any row/col member is blank
				return True
			if any([val not in digits for val in rowVals+colVals]): # any value not in [1,9]
				return False
			if len(rowVals) != len(set(rowVals)): # duplicates found in row
				return False
			if len(colVals) != len(set(colVals)):
				return False
			if sum(rowVals) != self.row.Sum or sum(colVals) != self.col.Sum:
				return False
			# return self.trialVal
			return True

# ******************** VECTOR class *****************
class Vector(object):
	def __init__(self,coords):
		self.xy = xy = coords
		self.isRow = xy[0][0] == xy[1][0]
		self.isCol = xy[0][1] == xy[1][1]
		self.Sum = 0
		self.get_cells()

	def __repr__(self):
		xy = self.xy
		if self.isRow:
			return 'row(%d,%d:%d)' % (xy[0][0]+1,xy[0][1]+1,xy[1][1]+1)
		elif self.isCol:
			return 'col(%d:%d,%d)' % (xy[0][0]+1,xy[1][0]+1,xy[1][1]+1)
		elif self.isRow and self.isCol:
			return 'cell(%d,%d)' % (xy[0][0]+1,xy[1][1]+1)
		else:
			return '--'

	def __eq__(self,other):
		return self.xy == other.xy

	def get_cells(self):
		self.cells = []
		if self.isRow:
			i = self.xy[0][0]
			jStart = self.xy[0][1]
			jEnd = self.xy[1][1]
			for j in range(jStart,jEnd+1):
				self.cells.append(A[i][j])
			for cell in self.cells:
				cell.row = self
		elif self.isCol:
			j = self.xy[1][1]
			iStart = self.xy[0][0]
			iEnd = self.xy[1][0]
			for i in range(iStart,iEnd+1):
				self.cells.append(A[i][j])
			for cell in self.cells:
				cell.col = self

	def one_var_ineq(self):
		"""Reduce possible values of vector cells by using the concept that 
		vector-wise pair values cannot be equal (e.g 'a'!='8-a')"""
		pairs = itertools.combinations(self.cells,2)
		for pair in pairs:
			expr1, expr2 = pair[0].val, pair[1].val
			vars1, vars2 = expr1.free_symbols, expr2.free_symbols
			if vars1 == vars2 and len(vars1) == 1: # handling only one-var inequations
				var = list(vars1)[0]
				val = sympy.solve(expr1 - expr2,var)[0]
				if val == int(val): # if value to be removed is an integer
					if val in Vars[var]['possibles']:
						Vars[var]['possibles'].remove(val)

	def two_var_ineq(self):
		"""Make pairs of variables whose difference cannot be some specified difference
		Type I: where two vars are just themselves (i.e. coeff = 1 & const = 0) [a != b]
		Type II: where both vars have +ve coeffs and a specified difference [a - b != k]
		Type III: where var coeffs can be +1 or -1 [k1*a + k2*b + k0 != 0]"""
		pairs = itertools.combinations(self.cells,2)
		varList = list(Vars.keys())
		for pair in pairs:
			expr1, expr2 = pair[0].val, pair[1].val
			var1, var2 = expr1.free_symbols, expr2.free_symbols
			if len(var1) == len(var2) == 1: # if each expr contains just one var
				var1 = list(var1)[0]
				var2 = list(var2)[0]
				if var1 != var2: # if the two variable sets are not the same
					i1 = varList.index(var1) # index of var1 in main varList
					i2 = varList.index(var2)
					k1 = int(expr1.coeff(var1)) # coeff of var1 (+1/0/-1)
					k2 = int(expr2.coeff(var2))
					k0 = int((expr1-expr2).evalf(subs={var1:0, var2:0})) # constant
			# Type I treatment: a != b
					# if k1 == k2 == 0:
					# 	pairIneqs.append((varList.index(var1),varList.index(var2)))
			# Type II treatment: a - b != k
					# if k1 == k2 == 1: # if each var's coeff. is 1 (e.g.'a-4', but not '5-a')
					# 		pairIneqs.append((i1,i2,-1*k0))
			# Type III treatment: k1*a + k2*b + k0 != 0
					pairIneqs.append({'x1':i1, 'x2':i2, 'k1':k1, 'k2':k2, 'k0':k0})

# ********************** MAIN ***********************
def main(fileName, digits, print_stats=False):
	minDigit = digits[0] # least attainable value of cell (= 1)
	maxDigit = digits[-1] # "" highest "" (=9)
	digRange = maxDigit - minDigit + 1

	t0 = time.clock()
	t1 = time.clock()
	print('Puzzle',fileName.capitalize().split('.')[0])
	# Creating the main data & cell object arrays
	data = importData(fileName)
	N = len(data)
	A = [] 	# 2D list of cell objects
	for i in range(N):
		row = []
		for j in range(N):
			row.append(Cell(data,i,j))
		A.append(row)

	# Creating the rows & cols of editables
	Rows = [] # list of row (vector) objects
	Cols = [] # "" col ""
	for i in range(N):
		for j in range(N):
			R = A[i][j].make_row(data)
			if R and R not in Rows: 
				Rows.append(R)
			C = A[i][j].make_col(data)
			if C and C not in Cols: 
				Cols.append(C)

	# Working with symbolic expressions
	nx = 0 # no. of generated & assigned sym vars
	Vars = {} # main sym var container
	flatList = [j for i in A for j in i] # 1D (flattened) list of A
	for cell in flatList:
		nx = cell.populate(nx) # filling cells with sym vars
	for cell in flatList:
		cell.get_coeffs() # getting coeff. of vars in cell value expression

	# Getting bounds of variables from each cell
	for cell in flatList:
		if not cell.blank:
			var, bounds = cell.get_bounds()
			if var:
				Vars[var]['bounds'].append(bounds)

	# Getting possible values of each variable overall
	for key,val in Vars.items():
		lowLim = max([x[0] for x in val['bounds']])
		upLim = min([x[1] for x in val['bounds']])
		size = upLim - lowLim + 1
		# print(key,':',[lowLim,upLim],size)
		Vars[key]['possibles'] = list(range(lowLim,upLim+1))

	# Removing one-variable inequations within rows & cols and
	# reducing two-variable interdependencies
	pairIneqs = []
	for vector in Rows + Cols:
		vector.one_var_ineq() # modifies ['possibles'] in Vars{}
		vector.two_var_ineq()

	# Compare effectiveness of this algo
	ranges = [v['possibles'] for v in Vars.values()]
	rangeSizes = [len(item) for item in ranges]
	nComb = np.prod(rangeSizes, dtype=np.int64)
	# totalComb = digRange**len(Vars.keys()) # total initial combs,i.e. 9^(#vars)
	print('nEdit = {}, nVars = {}, nComb = {:,}'.format(
		len([a for a in flatList if a.editable]),len(Vars),nComb))

	tPre = time.clock() -t1

	# ***************************************
	# Solve different combinations (hit & trial numerical computation)
	comb = (itertools.product(*ranges)) # generator for a particular combination
	tOp = 0 	# time of each 1% completion of operations
	tPerc = 0.  # time taken in calculating % completed
	t2var = 0. 	# time taken in removing 2-var expression pairs (e.g. a-4 != 3+b)
	tNext = 0.  # time taken in calculating next combination through generator
	tTry = 0. 	# time taken in calculating trial values from coeffs & combi/ns
	tCheck = 0. # time taken in checking each matrix combination
	n2var = 0 	# no. of combi/ns trashed due to simult. equality of 2-var expr/ns
	perc = 0 	# percentage operations completed
	flatList = [j for i in A for j in i]

	def check_mat():
		"""Checks if current trial value matrix is the solution"""
		vecCorrect = []
		for vec in Rows + Cols:
			trials = [cell.trialVal for cell in vec.cells]
			if any([val not in digits for val in trials]): # any row/col member not in [1,9]
				vecCorrect.append(False)
			elif len(trials) != len(set(trials)): # duplicate values in row/col
				vecCorrect.append(False)
			elif sum(trials) != vec.Sum: # sum of vector members not equal to vector sum
				vecCorrect.append(False)
			else:
				vecCorrect.append(True)
		if all(vecCorrect): # if all rows & cols satisfy the design constraints
			return True
		else:
			return False

	tNC = time.clock()
	for j in range(nComb): # try each combination
		# Jugaadu progress bar
		t7 = time.clock()
		if int(j/nComb*100) == perc:
			tOp = time.clock() - tOp
			print('{:3}%: j= {:}  t= {:.4}  nCompl= {:}'.format(perc,j,tOp,n2var))
			perc += 1
		tPerc += time.clock() - t7

		t8 = time.clock()
		C = np.array(next(comb))
		pairBreak = False # break current loop if this is true
		tNext += time.clock() - t8
		t2 = time.clock()
	# Type I: a != b
		# for tup in pairIneqs: # remove combinations that make for simult. equality
		# 	if C[tup[0]] == C[tup[1]]:
	# Type II: a - b != k
		# for tup in pairIneqs:
		# 	if C[tup[0]] - C[tup[1]] == C[tup[2]]:
	# Type III: k1*a + k2*b + k0 != 0
		for T in pairIneqs: # two member expressions cannot be equal
			if C[T['x1']]*T['k1'] - C[T['x2']]*T['k2'] + T['k0'] == 0:
				pairBreak = True
				break
		t2var += time.clock() - t2
		if pairBreak:
			n2var += 1
			continue
		
		t3 = time.clock()
		for cell in flatList:
			cell.try_sol(C)
		tTry += time.clock() - t3
		t4 = time.clock()
		isCorrect = check_mat() # is this the required solution
		tCheck += time.clock() - t4
		if isCorrect:
			print_mat(A)
			print('Stopped at j =',j)
			break
	tNC = time.clock() - tNC
	tt = time.clock() - t0
	if print_stats:
		print('No. of operations removed due to compl. =',n2var)
		print('Preprocessing time = {:.4}'.format(tPre))
		print('Time in calculating %. completed = {:.4}'.format(tPerc))
		print('Time taken in removing 2-var dependencies = {:.4}'.format(t2var))
		print('Time in calculating next combination through generator = {:.4}'.format(tNext))
		print('Time taken in calculating trial values = {:.4}'.format(tTry))
		print('Time taken in checking matrix combination = {:.4}'.format(tCheck))
		print('Total time in numerical computation = {:.4}'.format(tNC))
		print('Total time = {:.4}'.format(tt))

# Globals
fileName = r'Problems\easy2.csv'
digits = list(range(1,10))
# main(fileName, digits, print_stats=True)
