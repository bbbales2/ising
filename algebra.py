import sympy
import itertools
import collections
from sympy.matrices import Matrix, zeros
from sympy import simplify, expand, pprint, Poly

N1 = 3
N2 = 3

q_ = zeros(N1, N2)

for i in range(N1):
    for j in range(N2):
        q_[i, j] = sympy.symbols("q{0}{1}".format(i, j))

def q(i, j):
    #if (i, j) not in idxs:
    #    return 0
    
    #if i < N1 and i >= 0 and j < N2 and j >= 0:
    return q_[(i + N1) % N1, (j + N2) % N2]
    #else:
    #    return 0

idxs = list(itertools.product(range(N1), range(N2)))
#idxs = [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)]

squares = []
for i, j in idxs:
    squares.append((q(i, j) * q(i, j), 1))

scalars = []
for i, j in idxs:
    scalars.append(q(i, j))

print q_, idxs
pairs = []
for i, j in idxs:
    pairs.append((q(i, j) * q(i + 1, j) +
                  q(i, j) * q(i + 1, j + 1) +
                  q(i, j) * q(i, j + 1) +
                  q(i, j) * q(i - 1, j) +
                  q(i, j) * q(i - 1, j - 1) +
                  q(i, j) * q(i, j - 1)) / 2)

triplets = []
for i, j in idxs:
    triplets.append((q(i, j) * q(i - 1, j) * q(i - 1, j - 1) +
                     q(i, j) * q(i, j + 1) * q(i - 1, j) +
                     q(i, j) * q(i + 1, j + 1) * q(i, j + 1) +
                     q(i, j) * q(i + 1, j) * q(i + 1, j + 1) +
                     q(i, j) * q(i, j - 1) * q(i + 1, j) +
                     q(i, j) * q(i - 1, j - 1) * q(i, j - 1)) / 3)

scalars = sum(scalars)
pairs = sum(pairs)
triplets = sum(triplets)

prod = expand(scalars * pairs)
prod = simplify(prod.subs(squares))

def split(eq):
    terms = collections.defaultdict(list)
    for exp in prod.args:
        poly = Poly(exp)
        terms[len(poly.free_symbols)].append(exp)

    for order in terms:
        terms[order] = sum(terms[order])
        
    return terms

print prod
print triplets
#print split(prod)
out = sympy.gcd(prod, triplets)
print "Number of triplets in prod: ", out
quit()
q2, r2 = sympy.div(r1, scalars)
print "Number of scalars in prod:", q2, r2

quit()

terms = collections.defaultdict(list)
for exp in prod.args:
    poly = Poly(exp)
    terms[poly.LC()].append(exp)
for key in sorted(terms.keys()):
    print "Monomials with coefficient {0}:".format(key)
    for poly in terms[key]:
        is_triplet = False
        for trip in triplets.args:
            if len((poly / trip).free_symbols) == 0:
                is_triplet = True
        print("{0} {1}".format("t" if is_triplet else " ", poly))
    print ""
#pprint(prod.args(), order = 'grlex')
