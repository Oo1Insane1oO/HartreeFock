import sys
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

    def primitiveNorm(self, w, m):
        """ normalization factor for primitives """
        if (m <= - 0.5):
            return 0.0
        return np.sqrt(np.sqrt(w) / gamma(m+0.5))
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
        s = n+m
        if (s <= -1) or s%2:
            return 0.0
        # end if

        return self.primitiveNorm(w,n)*self.primitiveNorm(w,m) *\
                gamma((s+1)/2.)/np.sqrt(w)
    # end function overlapSolution

    def laplacianOverlap(self, w, n, m):
        """ return <g_n|nabla|g_m> """
        sums = np.zeros(3)
        for d in range(self.dim):
            tmpProdsdd = np.ones(3)
            for dd in range(self.dim):
                ndd = self.states[n,dd]
                mdd = self.states[m,dd]
                if dd != d:
                    tmpProdsdd *= w*self.overlapd(w,ndd,mdd)

                else:
                    tmpProdsdd[0] *= w*mdd*(mdd-1)*self.overlapd(w,ndd,mdd-2)
                    tmpProdsdd[1] *= w*(2*mdd+1)*self.overlapd(w,ndd,mdd)
                    tmpProdsdd[2] *= w*self.overlapd(w,ndd,mdd+2)
                    if (n==5 and m==5):
                        print w*self.overlapd(w,ndd,mdd), w*self.overlapd(w,ndd,mdd+2)
                # end ifelse
            # end fordd
            sums += tmpProdsdd
        # end ford

        return np.sum(sums*np.array([1,-1,1]))
    # end function laplacianOverlap

    def potentialOverlap(self, w, n, m):
        """ calculate and return <g_n|V(r)=1/2w^2r^2|g_m> """
        sum1 = 0
        for d in range(self.dim):
            tmpProd1 = 1
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

    def findCoefficients(self, w):
        """ find coefficients with equation HC=GCE, where C are the
        coefficients and the matrix elements are H_nm=<g_n|nabla+pot|g_m>,
        G_nm=<g_n|g_m>, E=diag({energy}) """

        H = np.zeros((self.size, self.size))
        G = np.zeros((self.size, self.size))

        # calculate matrix elements
        for i in range(self.size):
            for j in range(self.size):
                H[i,j] = -0.5*self.laplacianOverlap(w,i,j) +\
                        self.potentialOverlap(w,i,j)
                G[i,j] = self.overlap(w,i,j)
            # end forj
        # end fori

        # solve eigenvalue problem with scipy
        E, C = eigh(H, G)
        print "Epsilon: ", E, " Exact: ", self.states[:,-1]*w
        print self.states
        print "H:\n", H, "\n"
        print "G:\n", G, "\n"
        print "C:\n", C, "\n"

        return C 
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

    # read in basis
    gF = generalizedFit(filename, dim, cut)
    coeffs = gF.findCoefficients(w)

    # plot contracted function and hermite function
    x = np.linspace(-10,10,10000)*np.sqrt(w)
    y = np.linspace(-10,10,10000)*np.sqrt(w)
    e = np.exp(-1/2*(x**2+y**2))
#     g = coeffs[0][0]*e
#     h = hermite(np.sqrt(w)*x , 0)*hermite(np.sqrt(w)*y, 0)*e
#     g = (coeffs[0][0] + coeffs[1][1]*y)*e

#     hnorm = lambda n: np.sqrt(2**n*factorial(n)*np.sqrt(np.pi/w))
#     h = hermite(np.sqrt(w)*x , 0)*hermite(np.sqrt(w)*y, 1)*e / (hnorm(0)*hnorm(1))
#     print np.linalg.norm(g-h)
#     plt.title("Error: %f" % (np.linalg.norm(g-h)))
#     plt.plot(x, g, label="Contracted")
#     plt.plot(x, h, label="Hermite")
#             
#     plt.legend(loc="best")
#     plt.show()
# end ifmain
