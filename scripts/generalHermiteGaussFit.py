import sys
import numpy as np
import scipy.linalg
from scipy.misc import factorial, factorial2
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
                elements = elements[:dim] + [elements[self.dim+2]]
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
        self.states = np.array(tmptmptmp)
        self.size = len(self.states)
    # end function readBasis

    def overlap(self, w, n, m):
        """ return <g_n|g_m> """
        prod = np.longdouble(1.0)
        for d in range(self.dim):
            prod *= self.overlapd(w, self.states[n,d], self.states[m,d])
        # end fori
        return prod
    # end function overlap

    def overlapd(self, w, n, m):
        """ return <g_n|g_m> in 1d"""
        s = int(n+m)
        if s%2 or s<=-1:
            return 0.0
        # end if

        return np.sqrt((np.pi/w), dtype=np.longdouble) *\
                factorial(s)/factorial(s/2) * 0.5**s
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
                    tmpProdsdd *= w*self.overlapd(w,ndd,mdd)
                else:
                    tmpProdsdd[0] *= w*mdd*(mdd-1)*self.overlapd(w,ndd,mdd-2)
                    tmpProdsdd[1] *= w*(2*mdd+1)*self.overlapd(w,ndd,mdd)
                    tmpProdsdd[2] *= w*self.overlapd(w,ndd,mdd+2)
                # end ifelse
            # end fordd
            sumsd += tmpProdsdd
        # end ford

        return sumsd[0] - sumsd[1] + sumsd[2]
    # end function laplacianOverlap

    def potentialOverlap(self, w, n, m):
        """ calculate and return <g_n|V(r)=1/2w^2r^2|g_m> """
        sum1 = 0.0
        for d in range(self.dim):
            tmpProd1 = 1.0
            for dd in range(self.dim):
                ndd = self.states[n,dd]
                mdd = self.states[m,dd]
                if dd != d:
                    tmpProd1 *= self.overlapd(w,ndd,mdd)
                else:
                    tmpProd1 *= self.overlapd(w,ndd,mdd+2)
                # end ifelse
            # end fordd
            sum1 += tmpProd1
        # end ford

        return 0.5*w**2*sum1
    # end function potentialOverlap

    def makeHermites(self, N):
        """ make hermite functions for each level """
        r = np.array([np.linspace(-2.5,2.5,N) for i in range(dim)])
        hnorm = lambda n: (2**n*factorial(n)*(np.pi/w)**0.5)**0.5
        hermites = np.zeros((len(self.states),N))
        for i in range(len(self.states)):
            h = np.ones(N)
            for d in range(dim):
                h *= hermite(r[d]*np.sqrt(w), self.states[i,d]) *\
                        np.exp(-0.5*w*r[d]**2)/hnorm(self.states[i,d])
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
            if np.linalg.norm(psi[i] - h[i]) > 1e-13:
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

        oldLength = 0
        rotateCount = 0
        while True:
            psi = self.buildContracted(w, localC, N)
            diffIdx = self.makeDiffs(psi, N)

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
            localC = self.rotateCols(localC, diffIdx)
            rotateCount += 1
            oldLength = len(diffIdx)
        # end while True

        return localC
    # end function sortCoefficients

    def buildContracted(self, w, C, N):
        """ build contracted function """
        r = np.array([np.linspace(-2.5,2.5,N) for i in range(self.dim)])
        g = lambda x,n: x**n*np.exp(-0.5*x**2)
        psi = np.zeros((len(C),N))
        for i in range(len(C)):
            for j in range(len(C)):
                gj = np.ones(N)
                for d in range(self.dim):
                    gj *= g(r[d]*np.sqrt(w), self.states[j,d])
                # end ford
                psi[i] += C[j,i] * gj / gF.overlap(w,j,j)**0.5
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

        H = np.zeros((self.size, self.size))
        G = np.zeros((self.size, self.size))
       
        norms = np.zeros(self.size, dtype=np.longdouble)
        for i in range(self.size):
            overlapii = self.overlap(w,i,i)
            norms[i] = 1./np.sqrt((overlapii), dtype=np.longdouble)
            H[i,i] = (-0.5*self.laplacianOverlap(w,i,i) +
                    self.potentialOverlap(w,i,i))/w * np.square(norms[i],
                            dtype=np.longdouble)
            G[i,i] = 1.0
        # end fori

        for i in range(self.size):
            for j in range(i+1,self.size):
                normij = np.multiply(norms[i], norms[j], dtype=np.longdouble)
                H[i,j] = (-0.5*self.laplacianOverlap(w,i,j) +
                        self.potentialOverlap(w,i,j))/w * normij
                H[j,i] = H[i,j]
                G[i,j] = self.overlap(w,i,j) * normij 
                G[j,i] = G[i,j]
            # end forj
        # end fori

        with open("H.txt", "w") as oFile:
            for i in range(len(H)):
                for j in range(len(H)):
                    oFile.write(str(H[i,j]) + " ")
                oFile.write("\n")
        with open("G.txt", "w") as oFile:
            for i in range(len(H)):
                for j in range(len(H)):
                    oFile.write(str(G[i,j]) + " ")
                oFile.write("\n")

        print "H:\n", H, "\n"
        print "G:\n", G, "\n"

        # solve eigenvalue problem with scipy
        E, C = scipy.linalg.eigh(H, G)
#         E = scipy.linalg.eigvalsh(H, G)
#         C = scipy.linalg.solve(H, G*np.diag(E), sym_pos=True)

        return C, self.states, E
    # end function findCoefficients
# end class generalizedFit

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt

    np.set_printoptions(linewidth=1000000000000)

    try:
        filename = sys.argv[1]
        dim = int(sys.argv[2])
        w = float(sys.argv[3])
        cut = int(sys.argv[4])
    except IndexError:
        print "USAGE: python generateHermiteGaussFit.py <output filename> <num"
        " dimensions> <w>"
        sys.exit(0)
    # end try-except

    def swapCol(a, i, j):
        tmp = np.copy(a[:,i])
        a[:,i] = np.copy(a[:,j])
        a[:,j] = tmp
        return a
    # end function swap

    # read in basis
    gF = generalizedFit(filename, dim, cut)
    N = 1000
    coeffs, states, epsilon = gF.findCoefficients(w)
    coeffs = np.zeros((21,21))
    with open("C.txt", "r") as coeffFile:
        i = 0
        for line in coeffFile:
            row = line.split();
            for j,rw in enumerate(row):
                coeffs[i,j] = np.longdouble(rw);
            i += 1
    print coeffs
    coeffs = gF.sortCoefficients(coeffs, w, N)
    psi = gF.buildContracted(w, coeffs, N)
#     psi, coeffs, states, epsilon = gF.contractedFunction(w, N)
    print "Epsilon: ", epsilon, " Exact: ", states[:,-1]*w
    print "C:\n", coeffs, "\n"

    # plot contracted function and hermite function
    r = np.array([np.linspace(-2.5,2.5,N) for i in range(dim)])
    hermites = gF.makeHermites(N)

    for i in range(len(coeffs)):
        plt.plot(np.linalg.norm(r, axis=0), psi[i], label="psi%i off:%f" % (i,
            np.linalg.norm(psi[i]-hermites[i])))
#         plt.plot(np.linalg.norm(r, axis=0), hermites[i], label="H%i" % i,
#                 alpha=0.5)
    # end fori
    plt.xlabel('$r$')
    plt.ylabel('$\\psi(r)$')
    plt.legend(loc='upper right', bbox_to_anchor=(1.15,1.1))
#     plt.legend(loc='best')
    plt.show()
# end ifmain
