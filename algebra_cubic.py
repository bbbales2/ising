#%%
import sympy
import itertools
import collections
from sympy.matrices import Matrix, zeros
from sympy import simplify, expand, pprint, Poly

N1 = 4
N2 = 4

q_ = zeros(N1, N2)

for i in range(N1):
    for j in range(N2):
        q_[i, j] = sympy.symbols("q{0}{1}".format(i, j))

def q(i, j):
    return q_[(i + N1) % N1, (j + N2) % N2]

idxs = list(itertools.product(range(N1), range(N2)))

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
                  q(i, j) * q(i, j + 1) +
                  q(i, j) * q(i, j - 1) +
                  q(i, j) * q(i - 1, j)) / 2)

triplets = []
for i, j in idxs:
    triplets.append((q(i, j) * q(i - 1, j) * q(i - 1, j + 1) +
                     q(i, j) * q(i, j + 1) * q(i - 1, j + 1) +
                     q(i, j) * q(i, j + 1) * q(i + 1, j + 1) +
                     q(i, j) * q(i + 1, j) * q(i + 1, j + 1) +
                     q(i, j) * q(i + 1, j) * q(i + 1, j - 1) +
                     q(i, j) * q(i, j - 1) * q(i + 1, j - 1) +
                     q(i, j) * q(i, j - 1) * q(i - 1, j - 1) +
                     q(i, j) * q(i - 1, j) * q(i - 1, j - 1) +
                     q(i, j) * q(i - 1, j) * q(i, j + 1) +
                     q(i, j) * q(i, j + 1) * q(i + 1, j) +
                     q(i, j) * q(i + 1, j) * q(i, j - 1) +
                     q(i, j) * q(i, j - 1) * q(i - 1, j)) / 3)

scalars = sum(scalars)
pairs = sum(pairs)
triplets = sum(triplets)

prod = expand(scalars * pairs)
prod = simplify(prod.subs(squares))
#prod2 = simplify(expand(scalars * scalars * scalars).subs(squares))

#%%
def cdiv(a, b):
    for i in range(20, 0, -1):
        if(len((a - i * b).args) < len(a.args)):
            return i
    return 0

def split(eq):
    terms = collections.defaultdict(list)
    for exp in prod.args:
        poly = Poly(exp)
        terms[len(poly.free_symbols)].append(exp)

    for order in terms:
        terms[order] = sum(terms[order])

    return terms

def red(prod):
    q1 = cdiv(prod, triplets)
    r1 = prod - q1 * triplets
    print "Number of triplets in prod: ", q1
    q2 = cdiv(r1, scalars)
    r2 = r1 - q2 * scalars
    print "Number of scalars in prod:", q2
    return r2

prod = red(prod)
#prod2 = red(prod2)

#print "Remainder: ", r2
#%%

terms = collections.defaultdict(list)
for exp in prod.args:
    poly = Poly(exp)
    terms[poly.LC()].append(exp)

for key in sorted(terms.keys())[:2]:
    print "Monomials with coefficient {0}:".format(key)
    for poly in terms[key]:
        is_triplet = False
        for trip in triplets.args:
            if len((poly / trip).free_symbols) == 0:
                is_triplet = True
        print("{0} {1}".format("t" if is_triplet else " ", poly))
    print ""
#pprint(prod.args(), order = 'grlex')
