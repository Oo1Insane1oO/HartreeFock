import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from math import factorial
import scipy.linalg
import sys

from hermite import hermite

np.set_printoptions(linewidth=1000000000000)
    
def readBasis(filename, cut, dim):
    """ read in premade text file with quantum numbers in give basis """
    tmp = []
    with open(filename, "r") as ofile:
        for line in ofile:
            """ append only quantum numbers and energy to list """
            elements = line.split()
            elements = elements[:dim] + [elements[dim+2]]
            tmp.append([int(e) for e in elements])
        # end forline
    # end open filename
    tmp = tmp[:len(tmp)/2]

    tmptmp = [[tmp[0]]]
    currE = tmptmp[0][0][-1]
    j = 0
    for i in range(1,len(tmp)):
        if (currE == tmp[i][-1]):
            tmptmp[j].append(tmp[i])
        else:
            tmptmp.append([tmp[i]])
            currE = tmp[i][-1]
            j += 1
        # end ifelse
    # end forenis

    tmptmptmp = []
    i = 0
    while (i < cut+1):
        tmptmptmp += tmptmp[i]
        i += 1
    # end while

    return len(tmptmptmp), tmptmptmp
# end function readBasis

def sSquaredNorm(r):
    return sum([i**2 for i in r])

# end function sSquaredNorm

def gaussian(r, n, w):
    prod = 1
    for j,i in enumerate(r):
        prod *= (sp.sqrt(w)*i)**n[j]
    # end forji

    return prod * sp.exp(-sp.Rational(w,2)*sSquaredNorm(r))

# end function gaussian

def potential(r, w):
    return sp.Rational(w**2,2) * sSquaredNorm(r)
# end function potential

def laplacian(r, n, w):
    suml = 0
    sqrtw = sp.sqrt(w)
    for d in range(len(r)):
        prod = 1
        for dd in range(len(r)):
            if d != dd:
                prod *= (sqrtw*r[dd])**n[dd]
            else:
                prod *= w * (sqrtw*r[dd])**(n[dd]-2) * (n[dd]**2 - n[dd] -
                        2*n[dd]*w*r[dd]**2 + w**2*r[dd]**4
                                - w*r[dd]**2)
            # end ifelse
        # end fordd
        suml += prod
    # end ford

    return suml * sp.exp(-sp.Rational(w,2) * sSquaredNorm(r))
# end function laplacian

def rotateCols(C, colIdx, n=1):
    """ rotate columns in colIdx in matrix C by n """
    localC = np.copy(C)
    permIdx = colIdx[n:] + colIdx[:n]
    j = 0
    for i in range(len(C)):
        if j >= len(colIdx):
            break
        # end if

        if i == colIdx[j]:
            C[:,i] = localC[:,permIdx[j]]
            j += 1
        # end if
    # end fori

    return C
# end function rotateCols
    
def makeDiffs(psi, w, N, dim, states):
    """ make list of indices for which contracted is off """
    diffIdx = []
    r, h = makeHermites(w, N, dim, states)
    for i in range(len(h)):
        if np.linalg.norm(psi[i] - h[i]) > 1e-13: 
            diffIdx.append(i)
        # end if
    # end fori
    
    return diffIdx
# end function makeDiffs

def sortCoefficients(C, w, N, dim, states, norms):
    """ sort coefficients for degenerate eigenenergies """
    localC = np.copy(C)

    oldLength = 0
    rotateCount = 0
    while True:
        psi = buildContracted(w, localC, N, dim, states, norms)
        diffIdx = makeDiffs(psi, w, N, dim, states)

        if (len(diffIdx)==0):
            """ break if sorted """
            break
        # end if

        # reset rotation count if correct vector is found
        if (len(diffIdx) != oldLength):
            rotateCount = 0
        # end if

        # shift sign of eigenvector if rotated once 
        if rotateCount == len(diffIdx):
            for i in range(len(diffIdx)):
                localC[:,i] *= -1
            # end fori
            rotateCount = 0
        # end rotates
        
        # permute list once
        localC = rotateCols(localC, diffIdx)
        rotateCount += 1
        oldLength = len(diffIdx)
    # end while True

    return localC
# end function sortCoefficients

def buildContracted(w, C, N, dim, states, norms):
    r = np.array([np.linspace(-2.5,2.5,N) for i in range(dim)])
    g = lambda x,n: x**n * np.exp(-0.5*x**2)
    psi = np.zeros((len(C),N))
    for i in range(len(C)):
        for j in range(len(C)):
            gj = np.ones(N)
            for d in range(dim):
                gj *= g(np.sqrt(w)*r[d], states[j][d])
            # end ford
            psi[i] += C[j,i] * gj * norms[j]
        # end forj
    # end fori

    return psi
# end function buildContracted

def makeHermites(w, N, dim, states):
    r = np.array([np.linspace(-2.5,2.5,N) for i in range(dim)])
    hnorm = lambda n: (2**n*factorial(n)*(np.pi/w)**0.5)**0.5
    hermites = np.zeros((len(states),N))
    for i in range(len(states)):
        h = np.ones(N)
        for d in range(dim):
            h *= hermite(r[d]*np.sqrt(w), states[i][d]) *\
                    np.exp(-0.5*w*r[d]**2) / hnorm(states[i][d])
        # end ford
        hermites[i] = h
    # end fori

    return r, hermites
# end function makeHermites

def matZero(N,M):
    return sp.Matrix([[0 for j in range(N)] for i in range(M)])
# end function sZero

def vecZero(N):
    return sp.Matrix([0 for i in range(N)])
# end function vecZero

def matFill(N,M,val):
    return sp.Matrix([[val for j in range(N)] for i in range(M)])
# end function matFill

def vecFill(N,val):
    return sp.Matrix([val for i in range(N)])
# end function vecFill

def matDummy(N,M, D="D", r=False):
    A = []
    symbols = []
    for i in range(M):
        tmp = []
        for j in range(N):
            s = sp.symbols('%s_%i%i' % (D,i,j), real=r)
            tmp.append(s)
            symbols.append(s)
        # end forj
        A.append(tmp)
    # end fori
            
    return sp.Matrix(A), symbols
# end function matFill

def vecDummy(N, D="D", r=False):
    A = []
    symbols = []
    for i in range(M):
        s = sp.symbols('%s_%i' % (D,i), real=r)
        A.append(s)
        symbols.append(s)
    # end fori
            
    return sp.Matrix(A), symbols
# end function vecFill

def findCoefficients(w, N, dim, states):
    r = [sp.symbols('x%i' % i, real=True) for i in range(1,dim+1)]

#     H = np.zeros((numBasis,numBasis), dtype=np.longdouble)
#     G = np.zeros((numBasis,numBasis), dtype=np.longdouble)
# 
#     norms = np.zeros(numBasis, dtype=np.longdouble)

    numBasis = len(states)
    H = matZero(numBasis, numBasis)
    G = matZero(numBasis, numBasis)

    norms = vecFill(numBasis, 1)

    limits = tuple([(i, -sp.oo, sp.oo) for i in r])

    basisFuncs = [gaussian(r,states[i],w) for i in range(numBasis)]
    lapFuncs = [laplacian(r,states[i],w) for i in range(numBasis)]

    for i in range(numBasis):
        gii = sp.S(sp.integrate(basisFuncs[i]*basisFuncs[i], *limits))
        sqrtgii = sp.sqrt(gii)
        norms[i] /= sqrtgii
        basisFuncs[i] /= sqrtgii
        lapFuncs[i] /= sqrtgii
        
        # integrate over all dimensions
        lapii = -sp.Rational(1,2)*sp.S(sp.integrate(basisFuncs[i] *
            lapFuncs[i], *limits)) 
        potii = sp.S(sp.integrate(basisFuncs[i] * potential(r,w) *
            basisFuncs[i], *limits))

        H[i,i] = lapii + potii
        G[i,i] = sp.S(1.0)
    # end fori

    for i in range(numBasis):
        gi = basisFuncs[i]
        for j in range(i+1, numBasis):
            gj = basisFuncs[j]

            # integrate over all dimensions
            overlapij = sp.S(sp.integrate(gi * gj, *limits))
            lapij = -sp.Rational(1,2)*sp.S(sp.integrate(gi * lapFuncs[j],
                *limits))
            potij = sp.S(sp.integrate(gi * potential(r,w) * gj, *limits))

            Hij = lapij + potij
            H[i,j] = Hij
            G[i,j] = overlapij

            H[j,i] = Hij
            G[j,i] = overlapij
        # end forj
    # end fori

    print "H"
    sp.pprint(H)
    print "G"
    sp.pprint(G)

    Ginv = G**(-1)
    print "Ginv"
    sp.pprint(G**(-1))
    A = Ginv*H
    print "A"
    sp.pprint(A)
    sp.pprint(A.eigenvects())
# 
#     E, C = scipy.linalg.eigh(H, G)
#     print "Energy: ", E
#     print "C\n", C
    sys.exit(1)

    return E, C, norms
# end function findCoefficients

def swapCol(a, i, j):
    tmp = np.copy(a[:,i])
    a[:,i] = a[:,j]
    a[:,j] = tmp
    return a
# end function swap

filename = sys.argv[1]
dim = int(sys.argv[2])
w = float(sys.argv[3])
level = int(sys.argv[4])

numBasis, states = readBasis(filename, level, dim)

N = 1000

E, coeffs, norms = findCoefficients(w, N, dim, states);
psi = buildContracted(w, coeffs, N, dim, states, norms)
r, hermites = makeHermites(w, N, dim, states)

for i in range(len(coeffs)):
    plt.plot(np.linalg.norm(r, axis=0), psi[i], label="psi%i off:%f" % (i,
        np.linalg.norm(psi[i]-hermites[i])))
    plt.plot(np.linalg.norm(r, axis=0), hermites[i], label="H%i" % i,
            alpha=0.5)
    
plt.xlabel('$r$')
plt.ylabel('$\\psi(r)$')
# plt.legend(loc='upper right', bbox_to_anchor=(1.15,1.1))
plt.legend(loc='best')
plt.show()
