import math
import sys
from itertools import izip
from math import fabs
import itertools

# Polinomios sao representados como uma lista de coeficientes, exemplo: 2*x^3 + 3*x^2 + 1 := [1,0,3,2]

# Calcula divisao de dois polinomios num/den retorna o quociente e o resto em uma lista

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


# Testar se o reultado está correto: ([-4.0, 1.0], [-28.0, 7.0])
poly_divmod([-24,6,-4,1],[-1,0,1])



# Calcula a derivada formal de um polinomio
def derivative(poly):
	dx = [poly[i] * i for i in range(1, len(poly))]
	return dx

# Calcula a cadeia de sturm de um polinomio: cadeia composta dos restos da divisão do polinomio com a derivada dele http://www2.washjeff.edu/users/mwoltermann/Dorrie/24.pdf
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

# Verificar se o resultado está correto: [[-1, -3, 0, 0, 0, 1], [-3, 0, 0, 0, 5], [1.0, 2.4], [2.849295910493827]]
sturm_chain([-1,-3,0,0,0,1])


# Calcula polinomio em um ponto
def evaluate_poly(poly, x): 
	n, tmp = 0, 0
	for a in poly:
		tmp = tmp + (a * (x**n))
		n += 1
	return tmp

# Conta a quantidade de mudancas de sinais dos numeros em uma lista
def countSignChanges(seq):
	count = 0
	if (len(seq) > 1):
		for i in range(1,len(seq)):
			change = (seq[i] >= 0) ^ (seq[i-1] >= 0)   # xor para verificar se os sinais são iguais
			if (change == True):
				count += 1
	return count

# Verificar se o resultado está correto: no primeiro é 1 e no segundo é 2 
countSignChanges([0,-1])
countSignChanges([-1,0,-1])


# Recebe uma lista de tuplas, exemplo: [[-10, 1], [-9, 1], [-8, 0], [-7, 1], [-6, 1], [-5, 1], [-4, 1]] e verifica se teve mudança na segunda componente de cada tupla, por exemplo de [-9, 1] para [-8, 0] teve mudança de 1 para 0, verificar para todos os elementos e guardar as primeiras componentes, neste caso [-9,-8]
def verify_root(seq):
	acc = []
	for i in range(1,len(seq)):
		change = seq[i-1][1] - seq[i][1] 
		if (change > 0):
			pair = [seq[i-1][0] , seq[i][0]]
			acc.append(pair)
	return acc


# Isola todas as raizes de um polinomio em intervalos, procura no intervalo da reta [min,max] andando em sub-intervalos de h = 1
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


# Verificar se está correto: [[-2, -1], [1, 2]]
isolate_all_roots([-4,0,1],-10,10)
