import numpy as np
from scipy.linalg import eigh
from scipy.misc import factorial, factorial2
from scipy.special import gamma

from hermite import hermite

class generalizedFit:
    def __init__(self, filename, dim, cut):
        """ read in basis and set states """
        self.dim = dim
        self.states = []
        self.size = 0
        self.cut = cut
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
        self.states = np.array(self.states[:self.cut])
#         self.states = np.array(self.states[:len(self.states)/2])
        self.size = len(self.states)
    # end function readBasis

    def primitiveNorm(self, w, m):
        """ normalization factor for primitives """
        factor = np.sqrt(w/np.pi) * (2*w)**m / factorial2(2*m-1)
#         print factor, m
#         return (w/np.pi)**(0.5)
        return factor
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
        if (n+m <= -1) or ((n+m)%2):
            return 0.0
        s = n+m+1
#         return self.primitiveNorm(w, n) * self.primitiveNorm(w, m) * \
#                 (int((-1)**(float(n+m))) + 1) * factorial(s/2.) / (np.sqrt(w)*s)
        return self.primitiveNorm(w, n+m) * gamma(s/2.) / np.sqrt(w)
    # end function overlapSolution

    def laplacianOverlap(self, w, n, m):
        """ return <g_n|nabla|g_m> """
        sums = np.zeros(4)
        for d in range(self.dim):
            tmpProdsdd = np.ones(4)
            for dd in range(self.dim):
                tmpProdsddd = np.ones(4)
                for ddd in range(self.dim):
                    nddd = self.states[n,ddd]
                    mddd = self.states[m,ddd]
                    if ddd != d:
                        tmpProdsddd *= self.overlapd(w,nddd,mddd)
                    else:
                        tmpProdsddd[0] *= mddd*(mddd-1) * \
                                self.overlapd(w,nddd,mddd-2)
                        tmpProdsddd[1] *= self.overlapd(w,nddd,0)
                        tmpProdsddd[2] *= mddd*self.overlapd(w,nddd,mddd)
                        tmpProdsddd[3] *= self.overlapd(w,nddd,2)
                    # end ifelse
                # end forddd
                tmpProdsdd *= tmpProdsddd
            # end fordd
            sums += tmpProdsdd
        # end ford

        return -0.5*w**self.dim*np.sum(sums*np.array([1,-1,-1,-1]))
    # end function laplacianOverlap

    def potentialOverlap(self, w, n, m):
        """ calculate and return <g_n|V(r)|g_m> """
        sum1 = 0
        for d in range(self.dim):
            tmpProd1 = 1
            for dd in range(self.dim):
                tmpProd2 = 1
                for ddd in range(self.dim):
                    nddd = self.states[n,ddd]
                    mddd = self.states[m,ddd]
                    if ddd != d:
                        tmpProd2 *= self.overlapd(w,nddd,mddd)
                    else:
                        tmpProd2 *= self.overlapd(w,nddd,mddd+2)
                    # end ifelse
                # end forddd
                tmpProd1 *= tmpProd2
            # end fordd
            sum1 += tmpProd1
        # end ford

        return 0.5*w**2*sum1
    # end function potentialOverlap

    def findCoefficients(self, w):
        """ find coefficients with equation HC=GCE, where C are the
        coefficients and the matrix elements are H_nm=<g_n|nabla+pot|g_m>,
        G_nm=<g_n|g_m>, E=diag({energy}) """

        H = np.zeros((self.size, self.size))
        G = np.zeros((self.size, self.size))

        # calculate matrix elements
        for i in range(self.size):
            for j in range(self.size):
                H[i,j] = self.laplacianOverlap(w,i,j) + \
                        self.potentialOverlap(w,i,j)
                G[i,j] = self.overlap(w,i,j)
            # end forj
        # end fori

        # solve eigenvalue problem with scipy
        E, C = eigh(H, G)
        print E, self.states[:,-1]*w
        print
        print C

        return C
    # end function findCoefficients
# end class generalizedFit

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt

    try:
        filename = sys.argv[1]
        dim = int(sys.argv[2])
        w = float(sys.argv[3])
        cut = int(sys.argv[4])
    except IndexError:
        print "USAGE: python generateHermiteGaussFit.py <output filename> <num"
        "dimensions> <w>"
        sys.exit(0)
    # end try-except

    # read in basis
    gF = generalizedFit(filename, dim, cut)
    coeffs = gF.findCoefficients(w)

    # plot contracted function and hermite function
    x = np.linspace(-10, 10, 1000)
    y = np.linspace(-10,10,1000)
    e = np.exp(-w/2*(x**2+y**2))
    g = coeffs[0][0]*e
    h = hermite(np.sqrt(w)*x , 0)*hermite(np.sqrt(w)*y, 0)*e
    print np.linalg.norm(g-h)
#     plt.title("Error: %f" % (np.linalg.norm(g-h)))
#     plt.plot(x, g, label="Contracted")
#     plt.plot(x, h, label="Hermite")
#             
#     plt.legend(loc="best")
#     plt.show()
# end ifmain
