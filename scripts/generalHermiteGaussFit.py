import sys
import numpy as np
import scipy.linalg
from scipy.misc import factorial, factorial2
from scipy.special import gamma
import pickle

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
        self.states = np.array(tmptmptmp, dtype=np.int64)
        self.size = len(self.states)
    # end function readBasis

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

        return np.sum(np.multiply(sumsd, np.array([1,-1,1], dtype=np.longdouble),
            dtype=np.longdouble), dtype=np.longdouble)
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

        H = np.zeros((self.size, self.size), dtype=np.longdouble)
        G = np.zeros((self.size, self.size), dtype=np.longdouble)
       
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

        print "H:\n", H, "\n"
        print "G:\n", G, "\n"

        # solve eigenvalue problem with scipy
        E, C = scipy.linalg.eigh(H, G)

        return C, self.states, E
    # end function findCoefficients
# end class generalizedFit

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt

    np.set_printoptions(linewidth=1000000000000, precision=16)

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
    coeffs, states, epsilon = gF.findCoefficients(w)
    print "Epsilon: ", epsilon, " Exact: ", states[:,-1]*w

    tmp = np.copy(coeffs[:,3])
    coeffs[:,3] = coeffs[:,4]
    coeffs[:,4] = tmp
   
    if len(coeffs) > 5:
        tmp = np.copy(coeffs[:,6])
        coeffs[:,6] = coeffs[:,8]
        coeffs[:,8] = tmp

    print "C:\n", coeffs, "\n"

    # plot contracted function and hermite function
    N = 10000
    r = np.array([np.linspace(-2.5,2.5,N) for i in range(dim)])
    g = lambda x,n: x**n*np.exp(-0.5*x**2)
    psi = np.zeros((len(coeffs),N))
    for i in range(len(coeffs)):
        for j in range(len(coeffs)):
            gj = np.ones(N)
            for d in range(dim):
                gj *= g(r[d]*np.sqrt(w), states[j,d])
            # end ford
            psi[i] += coeffs[j,i] * gj / gF.overlap(w,j,j)**0.5
        # end forj
    # end fori

    hnorm = lambda n: (2**n*factorial(n)*(np.pi/w)**0.5)**0.5
    hermites = np.zeros((len(coeffs),N))
    for i,p in enumerate(psi):
        h = np.ones(N)
        for d in range(dim):
            h *= hermite(r[d]*np.sqrt(w), states[i,d]) * np.exp(-0.5*w*r[d]**2) /\
                    hnorm(states[i,d])
        # end ford
        hermites[i] = h
    # end forip

    for i in range(len(coeffs)):
#     for i in [1]:
        plt.plot(np.linalg.norm(r, axis=0), psi[i]**2, label="psi%i off:%f" %
                (i, np.linalg.norm(psi[i]**2-hermites[i]**2)))
        plt.plot(np.linalg.norm(r, axis=0), hermites[i]**2, 'o', markersize=0.15, label="H%i"
                % i, alpha=0.5)
    # end fori
    plt.xlabel('$r$')
    plt.ylabel('$|\\psi(r)|**2$')
    plt.legend(loc='upper right', bbox_to_anchor=(1.15,1.1))
    plt.show()
# end ifmain
