import sys
import numpy as np
import scipy.linalg
from math import factorial
from scipy.special import gamma
import itertools

from hermite import hermite

class generalizedFit:
    def __init__(self, filename, dim, cut):
        """ read in basis and set states """
        self.dim = dim
        self.states = []
        self.cut = cut
        self.size = 0
        self.readBasis(filename)
    # end __init__

    def readBasis(self, filename):
        """ read in premade text file with quantum numbers in give basis """
        tmp = []
        with open(filename, "r") as ofile:
            for line in ofile:
                """ append only quantum numbers and energy to list """
                elements = line.split()
                elements = elements[:self.dim] + [elements[self.dim+2]]
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
        while (i < self.cut+1):
            tmptmptmp += tmptmp[i]
            i += 1
        self.states = np.array(tmptmptmp, dtype=np.int64)
        self.size = len(self.states)
    # end function readBasis

    def overlap(self, w, n, m):
        """ return <g_n|g_m> """
        prod = 1
        for d in range(self.dim):
            prod = np.multiply(prod, self.overlapd(w, self.states[n,d],
                self.states[m,d]), dtype=np.longdouble)
        # end fori
        return prod
    # end function overlap

    def overlapd(self, w, n, m):
        """ return <g_n|g_m> in 1d"""
        s = np.int64(n+m)
        if s%2 or s<=-1:
            return np.longdouble(0.0)
        # end if

        return np.multiply(np.sqrt((np.pi/w), dtype=np.longdouble),
                np.multiply(np.divide(factorial(s),factorial(s/2),
                    dtype=np.longdouble), np.power(0.5,s, dtype=np.longdouble),
                    dtype=np.longdouble), dtype=np.longdouble)
#         return gamma((s+1)/2., dtype=float128)/np.sqrt(w, dtype=np.longdouble)
    # end function overlapSolution

    def laplacianOverlap(self, w, n, m):
        """ return <g_n|nabla|g_m> """
        sumsd = np.zeros(3, dtype=np.longdouble)
        for d in range(self.dim):
            tmpProdsdd = np.ones(3, dtype=np.longdouble)
            for dd in range(self.dim):
                ndd = self.states[n,dd]
                mdd = self.states[m,dd]
                if dd != d:
                    tmpProdsdd = np.multiply(tmpProdsdd, np.multiply(w,
                        self.overlapd(w,ndd,mdd), dtype=np.longdouble),
                        dtype=np.longdouble)
                else:
                    tmpProdsdd[0] = np.multiply(tmpProdsdd[0],
                            w*mdd*(mdd-1)*self.overlapd(w,ndd,mdd-2),
                            dtype=np.longdouble)
                    tmpProdsdd[1] = np.multiply(tmpProdsdd[1],
                            w*(2*mdd+1)*self.overlapd(w,ndd,mdd),
                            dtype=np.longdouble)
                    tmpProdsdd[2] = np.multiply(tmpProdsdd[2],
                            w*self.overlapd(w,ndd,mdd+2), dtype=np.longdouble)
                # end ifelse
            # end fordd
            sumsd += tmpProdsdd
        # end ford

        return np.sum(np.multiply(sumsd, np.array([1,-1,1],
            dtype=np.longdouble), dtype=np.longdouble), dtype=np.longdouble)
    # end function laplacianOverlap

    def potentialOverlap(self, w, n, m):
        """ calculate and return <g_n|V(r)=1/2w^2r^2|g_m> """
        sum1 = np.longdouble(0.0)
        for d in range(self.dim):
            tmpProd1 = 1
            for dd in range(self.dim):
                ndd = self.states[n,dd]
                mdd = self.states[m,dd]
                if dd != d:
                    tmpProd1 = np.multiply(tmpProd1, self.overlapd(w,ndd,mdd),
                            dtype=np.longdouble)
                else:
                    tmpProd1 = np.multiply(tmpProd1,
                            self.overlapd(w,ndd,mdd+2), dtype=np.longdouble)
                # end ifelse
            # end fordd
            sum1 += tmpProd1
        # end ford

        return np.multiply(np.longdouble(0.5), np.multiply(np.power(w,2,
            dtype=np.longdouble), sum1, dtype=np.longdouble),
            dtype=np.longdouble)
    # end function potentialOverlap

    def makeHermites(self, N):
        """ make hermite functions for each level """
        r = np.array([np.linspace(-2.5,2.5,N, dtype=np.longdouble) for i in
            range(dim)])
        hnorm = lambda n: (2**n*factorial(n)*(np.pi/w)**0.5)**0.5
        hermites = np.zeros((len(self.states),N), dtype=np.longdouble)
        for i in range(len(self.states)):
            h = np.ones(N, dtype=np.longdouble)
            for d in range(dim):
                h *= hermite(r[d]*np.sqrt(w, dtype=np.longdouble),
                        self.states[i,d]) * np.exp(-0.5*w*r[d]**2,
                                dtype=np.longdouble)/hnorm(self.states[i,d])
            # end ford
            hermites[i] = h
        # end forip
        
        return hermites
    # end function makeHermites

    def makeDiffs(self, psi, N):
        """ make list of indices for which contracted is off """
        diffIdx = []
        h = self.makeHermites(N)
        for i in range(len(h)):
            if np.linalg.norm(psi[i] - h[i], dtype=np.longdouble) > 5:
                diffIdx.append(i)
            # end if
        # end fori
        
        return diffIdx
    # end function makeDiffs

    def rotateCols(self, C, colIdx, n=1):
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

    def sortCoefficients(self, C, w, N):
        """ sort coefficients for degenerate eigenenergies """
        localC = np.copy(C)

        while True:
            psi = self.buildContracted(w, localC, N)
            diffIdx = self.makeDiffs(psi, N)

            if (len(diffIdx)==0):
                """ break if sorted """
                break
            # end if

            # permute list once
            localC = self.rotateCols(localC, diffIdx)
        # end while True

        return localC
    # end function sortCoefficients

    def buildContracted(self, w, C, N):
        """ build contracted function """
        r = np.array([np.linspace(-2.5,2.5,N, dtype=np.longdouble) for i in
            range(self.dim)])
        g = lambda x,n: np.multiply(np.power(x, n, dtype=np.longdouble),
                np.exp(-0.5*x**2, dtype=np.longdouble), dtype=np.longdouble)
        psi = np.zeros((len(C),N), dtype=np.longdouble)
        for i in range(len(C)):
            for j in range(len(C)):
                gj = np.ones(N, dtype=np.longdouble)
                for d in range(self.dim):
                    gj = np.multiply(gj, g(np.multiply(r[d], np.sqrt(w,
                        dtype=np.longdouble), dtype=np.longdouble),
                        self.states[j,d]), dtype=np.longdouble)
                # end ford
                psi[i] += np.divide(np.multiply(C[j,i], gj,
                    dtype=np.longdouble), gF.overlap(w,j,j)**0.5,
                    dtype=np.longdouble)
            # end forj
        # end fori

        return psi
    # end function buildContracted

    def contractedFunction(self, w, N):
        coeffs, states, epsilon = self.findCoefficients(w)
        coeffs = self.sortCoefficients(coeffs, w, N)
        return self.buildContracted(w, coeffs, N), coeffs, states, epsilon
    # end function contractedFunction

    def findCoefficients(self, w):
        """ find coefficients with equation HC=GCE, where C are the
        coefficients and the matrix elements are H_nm=<g_n|nabla+pot|g_m>,
        G_nm=<g_n|g_m>, E=diag({energy}) """

        H = np.zeros((self.size, self.size), dtype=np.longdouble)
        G = np.zeros((self.size, self.size), dtype=np.longdouble)
       
        norms = np.zeros(self.size, dtype=np.longdouble)
        for i in range(self.size):
            overlapii = self.overlap(w,i,i)
            norms[i] = np.divide(np.longdouble(1.0), np.sqrt((overlapii),
                dtype=np.longdouble), dtype=np.longdouble)
            H[i,i] = np.multiply(np.divide(np.multiply(-0.5,
                self.laplacianOverlap(w,i,i), dtype=np.longdouble) +
                self.potentialOverlap(w,i,i), w, dtype=np.longdouble),
                np.square(norms[i], dtype=np.longdouble), dtype=np.longdouble)
            G[i,i] = np.longdouble(1.0)
        # end fori

        for i in range(self.size):
            for j in range(i+1,self.size):
                normij = np.multiply(norms[i], norms[j], dtype=np.longdouble)
                H[i,j] = np.multiply(np.divide(np.multiply(-0.5,
                    self.laplacianOverlap(w,i,j), dtype=np.longdouble) +
                    self.potentialOverlap(w,i,j), w, dtype=np.longdouble),
                    normij, dtype=np.longdouble)
                H[j,i] = H[i,j]
                G[i,j] = np.multiply(self.overlap(w,i,j), normij,
                        dtype=np.longdouble)
                G[j,i] = G[i,j]
            # end forj
        # end fori

#         print "H:\n", H, "\n"
#         print "G:\n", G, "\n"

        # solve eigenvalue problem with scipy
        E, C = scipy.linalg.eigh(H, G)
        C = C.astype(np.longdouble)
        E = E.astype(np.longdouble)

        return C, self.states, E
    # end function findCoefficients
# end class generalizedFit

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt

    np.set_printoptions(linewidth=1000000000000, precision=19)

    try:
        filename = sys.argv[1]
        dim = np.int64(sys.argv[2])
        w = np.longdouble(sys.argv[3])
        cut = np.int64(sys.argv[4])
    except IndexError:
        print "USAGE: python generateHermiteGaussFit.py <output filename> <num"
        " dimensions> <w>"
        sys.exit(0)
    # end try-except

    def swapCol(a, i, j):
        tmp = np.copy(a[:,i])
        a[:,i] = a[:,j]
        a[:,j] = tmp
        return a
    # end function swap

    # read in basis
    gF = generalizedFit(filename, dim, cut)
    N = 10000
    coeffs, states, epsilon = gF.findCoefficients(w)
#     coeffs = swapCol(coeffs, 4, 5)
#     coeffs = swapCol(coeffs, 3, 4)
    coeffs = swapCol(coeffs, 6, 9)
    coeffs = swapCol(coeffs, 7, 8)
    coeffs = swapCol(coeffs, 7, 6)
    psi = gF.buildContracted(w, coeffs, N)
#     psi, coeffs, states, epsilon = gF.contractedFunction(w, N)
    print "Epsilon: ", epsilon, " Exact: ", states[:,-1]*w
    print "C:\n", coeffs, "\n"

    # plot contracted function and hermite function
    r = np.array([np.linspace(-2.5,2.5,N, dtype=np.longdouble) for i in
        range(dim)], dtype=np.longdouble)
    hermites = gF.makeHermites(N)

#     for i in range(len(coeffs)):
    for i in [6,7,8,9]:
        plt.plot(np.linalg.norm(r, axis=0), psi[i], label="psi%i off:%f" % (i,
            np.linalg.norm(psi[i]-hermites[i])))
        plt.plot(np.linalg.norm(r, axis=0), hermites[i], label="H%i" % i,
                alpha=0.5)
    # end fori
    plt.xlabel('$r$')
    plt.ylabel('$\\psi(r)$')
#     plt.legend(loc='upper right', bbox_to_anchor=(1.15,1.1))
    plt.legend(loc='best')
    plt.show()
# end ifmain
