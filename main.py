# Name: Sean Alfred A. Lipardo, 2018-04040 (Error Corrector)
# Partner: Jonathan Brito, 2019-21257 (Decoder)
# Section: CoE 164 WWX

import math
from copy import deepcopy

alpha = [1]
for i in range(0, 929):
	a = alpha[-1]*3
	if a > 928:
		alpha.append(a%929)
	else:
		alpha.append(a)



def inverse(number):
	t = 0
	next_t = 1
	p = 929

	while number > 0:
		q = p//number
		r = p%number
		t = t-mod(q*next_t)

		temp = next_t
		next_t = t
		t = temp

		p = number
		number = r

	return t%929


def mod(number):
	if number > 928 or number < 0:
		return number%929
	else:
		return number


def syndrome(SCV, ecc_level):
	global alpha
	numofecc = 2**(ecc_level+1)
	temp = []
	returnSyndrome = []

	for i in range(1, numofecc+1):
		exp = len(SCV)-1
		for j in range(0, len(SCV)):
			elementVal = mod(SCV[j]*mod(alpha[exp]**i))
			temp.append(elementVal)
			exp = exp-1

		sumTemp = mod(sum(temp))
		returnSyndrome.append(sumTemp)
		temp = []

	# print("Syndrome:", returnSyndrome)
	return returnSyndrome


def errorLocatorPolynomial(syndrome):
	global aplha
	lx = [1]
	lp = [1]
	ne = 0
	dp = 1
	m = 1

	for i in range(0, len(syndrome)):
		# Calculate discrepancy
		for j in range(0, ne+1):
			if j == 0:
				d = syndrome[i]
				continue
	
			evaluate = mod(syndrome[i-j] * lx[j])
			d = mod(d + evaluate)

		if d == 0:
			m += 1
		else:
			oldlx = deepcopy(lx)
			templp = deepcopy(lp)
			templp[:0] = [0]*(m)

			inverseDp = inverse(dp)
			ddp = mod(-d*inverseDp)
			ddpXlp = [mod(ddp*k) for k in templp]

			if len(ddpXlp)>len(lx):
				lx.extend([0]*(len(ddpXlp)-len(lx)))
			else:
				ddpXlp.extend([0]*(len(lx)-len(ddpXlp)))

			templx = [mod(l+n) for l,n in zip(lx, ddpXlp)]
			lx = deepcopy(templx)

			if 2*ne <= i:
				lp = oldlx
				ne = i+1-ne
				dp = d
				m = 1
			else:
				m += 1

	# print("LX:", lx)
	# print("ne:", ne)
	return lx


def elpRoots(LX):
	global alpha
	ne = len(LX)-1
	T = [alpha[i] for i in range(0, ne+1)]
	elpRoots = []

	for i in range(0, 929):
		if i == 0:
			le = deepcopy(LX)
		else:
			temple = [mod(j*k) for j,k in zip(le, T)]
			le = deepcopy(temple)

		elp_val = mod(sum(le))
		if elp_val == 0:
			elpRoots.append(i)
			# print("zero at:", i)
		
	# print("ELP Roots:", elpRoots)
	return elpRoots


def errorPolynomial(syndrome, lx, rootsLx):
	# Derivative
	ne = len(lx)-1
	dlx = [0]*ne
	result = [0]*(len(syndrome)+ne)

	for i in range(1, len(lx)):
		tempEval = mod(i*lx[i])
		dlx[i-1] = tempEval
	# print("dlx:", dlx)


	# Error Evaluator
	for i in range(len(syndrome)):
		for j in range(len(lx)):
			result[i+j] = mod(result[i+j]+syndrome[i]*lx[j])
	omega = result[0:ne]
	# print("omega:", omega)


	# Coefficients of e(x)
	e_coeffs = []
	for root in rootsLx:
		derivatives = []
		omegas = []
		tempne = 0
		for i,j in zip(dlx, omega):
			x = mod(alpha[root]**tempne)
			derivatives.append(mod(i*x))
			omegas.append(mod(j*x))
			tempne += 1

		# print("Derivatives:", derivatives)
		# print("Omegas:", omegas)

		sumOfDerivative = sum(derivatives)
		sumOfOmegas = sum(omegas)

		inverseOfDerivative = inverse(sumOfDerivative)
		modOfOmegas = mod(-sumOfOmegas)

		# print(inverseOfDerivative, modOfOmegas)

		e_coeffs.append(mod(modOfOmegas*inverseOfDerivative))


	# print("Error Polynomial:", e_coeffs)
	return e_coeffs


def trueMessage(m, rootsLx, e_coeffs):
	global alpha
	padded = [0]*len(m)

	for k, l in enumerate(rootsLx):
		rootsLx[k] = inverse(alpha[l])

	for i,j in zip(rootsLx, e_coeffs):
		# print("found:", i, "at index:", alpha.index(i))
		
		padded[alpha.index(i)] = j
		
	padded = padded[::-1]
	returnM = [mod(i+mod(-j)) for i,j in zip(m, padded)]

	return returnM


def hlValues(SCV, N):
	numofecc = 2**(N+1)
	num_of_data = len(SCV)-(numofecc)
	returnhlSCV = []
	for i in range(1, num_of_data):
		# print("i:", i, "Size of Data:", SCV[0])
		if SCV[i]==900:
			continue
		H = math.floor(SCV[i]/30)
		L = SCV[i] - 30*H
		returnhlSCV.extend([H, L])

	# print(returnhlSCV)
	return returnhlSCV


def barcodeDecoder(SCV, N):
	hlSCV = hlValues(SCV, N)
	codewords = {0: ["A", "a", "0", ";"], 1: ["B", "b", "1", "<"], 2: ["C", "c", "2", ">"], 3: ["D", "d", "3", "@"], 4: ["E", "e", "4", "["], 5: ["F", "f", "5", '\\'], 6: ["G", "g", "6", "]"], 7: ["H", "h", "7", "_"], 8: ["I", "i", "8", "`"], 9: ["J", "j", "9", "~"], 10: ["K", "k", "&", "!"], 11: ["L", "l", "\r", "\r"], 12: ["M", "m", "\t", "\t"], 13: ["N", "n", ",", ","], 14: ["O", "o", ":", ":"], 15: ["P", "p", "#", "\n"], 16: ["Q", "q", "-", "-"], 17: ["R", "r", ".", "."], 18: ["S", "s", "$", "$"], 19: ["T", "t", "/", "/"], 20: ["U", "u", "+", '"'], 21: ["V", "v", "%", "|"], 22: ["W", "w", "*", "*"], 23: ["X", "x", "=", "("], 24: ["Y", "y", "^", ")"], 25: ["Z", "z", "Lp", "?"], 26: [" ", " ", " ", "{"], 27: ["LL", "Sa", "LL", "}"], 28: ["Lm", "Lm", "La", "'"], 29: ["Sp", "Sp", "Sp", "La"]}
	x = 0
	temp = 0
	flag = 0
	returnDecoded = ""

	for i in hlSCV:
		if codewords[i][x] == "Sp":
			temp = x
			x = 3
			flag = 1
			continue
		elif codewords[i][x] == "Sa":
			temp = x
			x = 0
			flag = 1
			continue
		elif codewords[i][x] == "La":
			x = 0
			continue
		elif codewords[i][x] == "LL":
			x = 1
			continue
		elif codewords[i][x] == "Lm":
			x = 2
			continue
		elif codewords[i][x] == "Lp":
			x = 3
			continue

		if flag:
			returnDecoded += codewords[i][x]
			x = temp
			flag = 0
		else:
			returnDecoded += codewords[i][x]

	return returnDecoded

def master(message, N, size):
	synd = syndrome(message, N)
	elp = errorLocatorPolynomial(synd)
	elp_roots = elpRoots(elp)
	error_poly = errorPolynomial(synd, elp, elp_roots)
	true_msg = trueMessage(message, elp_roots, error_poly)
	# print("message:", message)
	# print("true message:", true_msg)
	decoded = barcodeDecoder(true_msg, N)
	num_of_changes = sum(1 for i, j in zip(message, true_msg) if i != j)
	
	return  num_of_changes, true_msg, decoded




input_T = int(input())

for i in range(input_T):
	x = input()
	ec_level, size_of_data = x.split(" ")
	y = input()
	input_message = y.split(" ")
	int_message = [int(k) for k in input_message]

	output = master(int_message, int(ec_level), int(size_of_data))
	num_of_changes = str(output[0])
	decoded_msg = output[2]
	true_message = ""
	for j in output[1]:
		true_message = true_message + str(j) + " "

	print("Case #" + str(i+1) + ":")
	print(num_of_changes, true_message)
	print(decoded_msg)

