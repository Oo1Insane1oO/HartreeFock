import numpy as np
from scipy.linalg import eigh
from scipy.misc import factorial, factorial2

# import matplotlib.pyplot as plt

class generalizedFit:
    def __init__(self, filename, dim):
        """ read in basis and set states """
        self.dim = dim
        self.states = []
        self.size = 0
        self.readBasis(filename)
    # end __init__

    def readBasis(self, filename):
        """ read in premade text file with quantum numbers in give basis """
        with open(filename, "r") as ofile:
            for line in ofile:
                """ append only quantum numbers and energy to list """
                elements = line.split()
                elements = elements[:dim] + [elements[self.dim+2]]
                self.states.append([int(e) for e in elements])
            # end forline
        # end open filename
#         self.states = np.array(self.states[:len(self.states)/2])
        self.states = np.array(self.states[:len(self.states)/2])
        self.size = len(self.states)
    # end function readBasis

    def primitiveNorm(self, a, m):
        """ normalization factor for primitives """
        return np.sqrt(np.sqrt(2*a/np.pi) * (4**a)**m / factorial2(2*m-1))
    # end function norm

    def overlap(self, w, n, m):
        """ return <g_n|g_m> """
        prod = 1
        for d in range(self.dim):
            prod *= self.overlapd(w, self.states[n,d], self.states[m,d])
        # end fori
        return prod
    # end function overlap

    def overlapd(self, w, n, m):
        """ return <g_n|g_m> in 1d"""
        if ((n+m)%2):
            return 0
        s = n+m+1
#         return self.primitiveNorm(w/2., n) * self.primitiveNorm(w/2., m) * \
#                 1./np.sqrt(w) * 2. / s * factorial(s/2.)
        return 1./np.sqrt(w) * 2. / s * factorial(s/2.)
    # end function overlapSolution

    def laplacianOverlap(self, w, n, m):
        """ return <g_n|nabla|g_m> """
        sum1 = w*np.sum(self.states[m]*(self.states[m]+1)) * \
                self.overlap(w,n,m) 
        sum2 = 0
        sum3 = 0
        for d in range(self.dim):
            tmpsum2 = 1
            tmpsum3 = 1
            for dd in range(self.dim):
                if (dd != d):
                    md = self.states[m,d]
                    tmpsum2 *= md*(md-1) * self.overlapd(w, self.states[n,d],
                            md-2)
                    tmpsum3 *= self.overlapd(w, self.states[n,d], 1)
                else:
                    tmpoverlapd = self.overlapd(w, self.states[n,dd],
                            self.states[m,dd])
                    tmpsum2 *= tmpoverlapd
                    tmpsum3 *= tmpoverlapd
                # end ifelse
            sum2 += tmpsum2
            sum3 += tmpsum3
        # end ford

        return w * (sum1 - sum2 - sum3)
    # end function laplacianOverlap

    def findCoefficients(self, w):
        """ find coefficients with equation HC=GCE, where C are the
        coefficients and the matrix elements are H_nm=<g_n|nabla+pot|g_m>,
        G_nm=<g_n|g_m>, E=diag({energy}) """

        H = np.zeros((self.size, self.size))
        G = np.zeros((self.size, self.size))
#         exa = np.diagflat(self.states[:,-1])

        # calculate matrix elements
        for i in range(self.size):
            for j in range(self.size):
                H[i,j] = self.overlap(w, i, j)
                G[i,j] = self.laplacianOverlap(w, i, j)
            # end forj
        # end fori

        # solve eigenvalue problem with scipy
        E, C = eigh(H, G)
        print C
    # end function findCoefficients
# end class generalizedFit

if __name__ == "__main__":
    import sys

    try:
        filename = sys.argv[1]
        dim = int(sys.argv[2])
        w = float(sys.argv[3])
    except IndexError:
        print "USAGE: python generateHermiteGaussFit.py <output filename> <num"
        "dimensions> <w>"
        sys.exit(0)
    # end try-except

    # read in basis
    gF = generalizedFit(filename, dim)
    gF.findCoefficients(w)
# end ifmain
