import math
import sys
from itertools import izip
from math import fabs
import itertools


def normalize(poly):
    while poly and poly[-1] == 0:
        poly.pop()
    if poly == []:
        poly.append(0)


def poly_divmod(num, den):
    num = num[:]
    normalize(num)
    den = den[:]
    normalize(den)
    if len(num) >= len(den):
        shiftlen = len(num) - len(den)
        den = [0] * shiftlen + den
    else:
        return [0], num
    quot = []
    divisor = float(den[-1])
    for i in xrange(shiftlen + 1):
        mult = num[-1] / divisor
        quot = [mult] + quot
        if mult != 0:
            d = [mult * u for u in den]
            num = [u - v for u, v in zip(num, d)]
        num.pop()
        den.pop(0)
    normalize(num)
    return quot, num

 
def derivative(poly):
	dx = [poly[i] * i for i in range(1, len(poly))]
	return dx


def sturm_chain(poly):
	chain = []
	chain.append(poly)
	dx_poly =derivative(poly) 
	chain.append(dx_poly)
	remainder = []
	i = 0
	while (len(remainder) != 1):
		remainder = poly_divmod(chain[i],chain[i+1])[1]
		remainder = [ -x for x in remainder]
		chain.append(remainder)
		i+=1
	return chain


def evaluate_poly(poly, x): 
	n, tmp = 0, 0
	for a in poly:
		tmp = tmp + (a * (x**n))
		n += 1
	return tmp


def countSignChanges(seq):
	count = 0
	if (len(seq) > 1):
		for i in range(1,len(seq)):
			change = (seq[i] >= 0) ^ (seq[i-1] >= 0)   # xor para verificar se os sinais são iguais
			if (change == True):
				count += 1
	return count


def verify_root(seq):
	acc = []
	for i in range(1,len(seq)):
		change = seq[i-1][1] - seq[i][1] # xor para verificar se os sinais são iguais
		if (change > 0):
			pair = [seq[i-1][0] , seq[i][0]]
			acc.append(pair)
	return acc



def isolate_all_roots(poly, min, max):
	h = 1
	chain = sturm_chain(poly)
	point = min
	deltas = []
	while (point <= max):
		tmp = [evaluate_poly(chain_element,point) for chain_element in chain]
		pair = [point, countSignChanges(tmp)]
		deltas.append(pair)
		point += h
	intervals = verify_root(deltas)
	return intervals
