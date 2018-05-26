import numpy as np
import sympy as sp
from scipy.special import factorial, gamma
from scipy.linalg import eigh

# SPECIALIZED TO OMEGA=1.0

class Fit:
    def __init__(self, R, dim, basisFile, orderDW=1, perturbAxis=0):
        self.R = float(R)
        self.dim = int(dim)

        self.orderDW = int(orderDW)

        self.perturbAxis = int(perturbAxis)

        self.reinitialize(dim, basisFile)
    # end __init__

    def readBasis(self, filename):
        basisList = []
        with open(filename, "r") as basisFile:
            for line in basisFile:
                basisList.append(map(int, line.split()))
            # end for line
        # end with open

        return np.array(basisList[:len(basisList)/2], dtype=np.int64)
    # end functions readBasis
    
    def setNormalizationFactors(self):
        """ calculate normalizations for each basis function """
        self.normalizationFactors = np.ones(len(self.basisList))
        for i in range(len(self.basisList)):
            for d in range(self.dim):
                nd = self.basisList[i,d]
                self.normalizationFactors[i] *= 1.0 / np.sqrt(np.sqrt(np.pi) *\
                        2**nd * factorial(nd))
            # end ford
        # end fori
    # end function setNormalizationFactors

    def reinitialize(self, dim, basisFile):
        """ initialize everything """
        self.dim = dim
        self.basisList = self.readBasis(basisFile)
        self.coeffs = self.hermiteCoefficients(np.max(self.basisList[:,0]))
        self.setNormalizationFactors()
    # end function reinitialize 
    
    def hermiteCoefficients(self, n, physisists=True):
        """ use recursive relation for coefficients and return a table with
        coefficients each from 0 to n """
        table = [[1], [0,1]]
        while (n>=len(table)):
            i = len(table)
            preva = table[i-1]
            ppreva = table[i-2]
            tmp = [0 for i in range(i+1)]
            m = 1 - i
            tmp[0] =  m*ppreva[0]
            tmp[i-1] = preva[i-2]
            tmp[i] = preva[i-1]
            for k in range(1,i-1):
                tmp[k] = preva[k-1] + m*ppreva[k]
            # end fork
            table.append(tmp)
        # end while

        if physisists:
            """ convert to physicists' """
            for i,t in enumerate(table):
                tmp = sp.sympify(2**(sp.Rational(i,2)))
                for k,tt in enumerate(t):
                    table[i][k] = int(tmp*sp.sympify(2**(sp.Rational(k,2)))*tt)
                # end forktt
            # end forit
        # end if
        return table
    # end function hermiteCoefficients

    def overlapd(self, n, m):
        s = n+m
        if ((s<=-1) or (s%2==1)):
            """ not defined for negative values """
            return 0.0
        # end if

        return gamma((s+1)/2.)
    # end function overlapd

    def overlapdDW(self, n, m):
        s = n+m
        if ((s<=-1) or (s%2==1)):
            """ not defined for negative values """
            return 0.0
        # end if

        return gamma((s+self.orderDW+1.)/2.)
    # end function overlapdDW

    def ddexpr(self, nd, md, f):
        sums = 0.0
        for p in range(nd+1):
            for q in range(md+1):
                sums += self.coeffs[nd][p] * self.coeffs[md][q] * f(p,q)
            # end forq
        # end forp

        return sums
    # end function ddexpr

    def ddexprOverlap(self, p, q):
        return self.overlapd(p,q)
    # end function ddexprOverlap

    def ddexprPotDW(self, p, q):
        return self.overlapdDW(p,q)
    # end function ddexprOverlap
    
    def potDWElement(self, i, j):
        tmpProdsd = 1.0
        for d in range(self.dim):
            nd = self.basisList[i,d]
            md = self.basisList[j,d]
            if d != self.perturbAxis:
                tmpProdsd *= self.ddexpr(nd, md, self.ddexprOverlap)
            else:
                tmpProdsd *= self.ddexpr(nd, md, self.ddexprPotDW)
            # end ifelse
        # end ford

        res = - 0.5*self.R * tmpProdsd * self.normalizationFactors[i] *\
                self.normalizationFactors[j]
        if i==j:
            res += 1/8. * R**2
        # end if

        return res
    # end function potDWElement

    def HDWij(self, i, j):
        """ return <psi_i|h_DW|psi_j>, with psi_i and psi_j being basis
        function i and j respectively in basisList """
        res = self.potDWElement(i,j)
        if i==j:
            """ kinetic + HO potential part """
            if (self.dim == 2):
                res += self.basisList[i,-2]
            else:
                res += 3/2. + (self.basisList[i,-2] - 1)
        # end if

        return res
    # end function HDWij

    def findCoefficients(self):
        """ build H_ij and solve HC=EC using scipy eigh """
        n = len(self.basisList)
        H = np.zeros((n,n))
        for i in range(n):
            for j in range(i,n):
                H[i,j] = self.HDWij(i,j)
                H[j,i] = H[i,j]
            # end forj
        # end fori

        return eigh(H)
    # end function findCoefficients

    def writeToFile(self, C2D, C3D, E2D, E3D, trunc2D, trunc3D):
        n2D = len(C2D)
        n3D = len(C3D)
        with open("dwc.h", "a") as writeFile:
            writeFile.write("#ifndef DWC_H\n#define DWC_H\n")
            writeFile.write("\n#include <array>\n")
            writeFile.write("#include <Eigen/Sparse>\n")
            writeFile.write(
                    "\nclass DWC {\n"\
                    "   private:\n")
        # end with open

        self.writePrivate(C2D, E2D, "2D", trunc2D)
        self.writePrivate(C3D, E3D, "3D", trunc3D)

        with open("dwc.h", "a") as writeFile:
            writeFile.write("\n\n"
                    "   public:\n"
                    "       DWC(const unsigned int dim) {\n"
                    "           if (dim==2) {\n"
                    "               C = Eigen::Map<const Eigen::MatrixXd>(m_C2D.data(),%i,%i).sparseView(1,1e-8);\n"
                    "               epsDW = Eigen::Map<const Eigen::ArrayXd>(m_epsDW2D.data(),%i);\n"
                    "           } else if (dim==3) {\n"
                    "               C = Eigen::Map<const Eigen::MatrixXd>(m_C3D.data(),%i,%i).sparseView(1,1e-8);\n"
                    "               epsDW = Eigen::Map<const Eigen::ArrayXd>(m_epsDW3D.data(),%i);\n"
                    "           } // end ifeif\n\n"
                    "       } // end constructor\n" % (n2D, trunc2D, trunc2D, n3D, trunc3D, trunc3D))
        # end with open
            
        with open("dwc.h", "a") as writeFile:
            writeFile.write("       virtual ~DWC() {}\n\n"
                            "       Eigen::SparseMatrix<double> C;\n"
                            "       Eigen::ArrayXd epsDW;\n")
            writeFile.write("};\n\n#endif\n")
    # end function writeDefine

    def writePrivate(self, C, E, ext, trunc):
        """ write coefficients matrix C to C++ header DWC """
        with open("dwc.h", "a") as writeFile:
            n = len(C)
            writeFile.write("       static constexpr std::array<double,%i> m_epsDW%s = {\n"
                % (n,ext))
            for i in range(trunc):
                if i==n-1:
                    writeFile.write("                   " + str(E[i]))
                else:
                    writeFile.write("                   " + str(E[i]))
                    writeFile.write(",\n")
                # end ifelse
            # end fori
            writeFile.write("\n             };\n")
            writeFile.write("       static constexpr std::array<double,%i> m_C%s = {\n" % (n*trunc,ext))
            writeFile.write("               ")
            for i in range(n):
                for j in range(trunc):
                    if ((i==n-1) and (j==trunc-1)):
                        writeFile.write(str(C[i,j]))
                    else:
                        writeFile.write(str(C[i,j]))
                        writeFile.write(", ")
                    # end ifelse
                # end forj
                if i != n-1:
                    writeFile.write("\n               ")
                else:
                    writeFile.write("\n             };\n")
            # end fori
        # end with open
    # end function writePrivate
# end class Fit

if __name__ == "__main__":
    """ find coefficients """
    import sys
    import subprocess as subpr
    import os.path

    try:
        R = float(sys.argv[1])
        numBasis2D = int(sys.argv[2])
        numBasis3D = int(sys.argv[3])
        basisfname2D = sys.argv[4]
        basisfname3D = sys.argv[5]
    except IndexError:
        print "USAGE: python fit.py 'R' 'cutoff 2D' 'cutoff 3D' 'filename 2D' 'filename 3D'"
        sys.exit()
    # end try-except

    def compileBasis(n,d,basisfname):
        if not os.path.isfile(basisfname):
            """ check if basis file already exists """
            if not os.path.isfile("main"):
                """ compile first """
                print "Compiling..."
                args = ["g++", "-std=c++17", "-O3", "cartesian.cpp" , "main.cpp" ,
                        "-o", "main"]
                print "".join(map(lambda v: v+" ", args))
                subpr.call(args)
            # end if
            print "\nCreating basis file with %s..." % basisfname
            args = ["./main", str(n), str(d), basisfname]
            print "".join(map(lambda v: v+" ", args))
            subpr.call(args)
        # end if
    # end function compileBasis

    def checkFullShell(n,fit):
        if n not in fit.basisList[:,-1]:
            """ make sure shell is filled """
            print "Shell not filled!\nTry one of these(increase beyond last " \
                    "option for more alternatives): " + "".join(map(lambda v:
                        str(v)+" ", sorted(set(fit.basisList[:,-1]))))
            sys.exit()
        # end if
    # end function checkFullShell

    def writeEigenValuesToFile(E2D, E3D, fname2D, fname3D):
        """ write energies to text file """
        with open(fname2D, "w") as file2d:
            for e2 in E2D:
                file2d.write(str(e2) + " ")
            # end fore2
        # end withopen
        with open(fname3D, "w") as file3d:
            for e3 in E3D:
                file3d.write(str(e3) + " ")
            # end fore2
        # end withopen
    # end function writeEigenValuesToFile
   
    # write both 2D and 3D
    # write basis file first
    compileBasis(numBasis2D, 2, basisfname2D)
    fit = Fit(R, 2, basisfname2D)
    checkFullShell(numBasis2D, fit)

    print "\nFinding coefficients 2D..."
    e2D, C2D = fit.findCoefficients()

    print "\nEigenenergies 2D: ", e2D
    print "Coefficients matrix 2D:"
    print C2D

    compileBasis(numBasis3D, 3, basisfname3D)
    fit.reinitialize(3, basisfname3D)
    checkFullShell(numBasis3D, fit)

    print "\nFinding coefficients 3D..."
    e3D, C3D = fit.findCoefficients()

    print "\nEigenenergies 3D: ", e3D
    print "Coefficients matrix 3D:"
    print C3D

    fit.writeToFile(C2D, C3D, e2D, e3D, int(sys.argv[8]), int(sys.argv[9]))

    if sys.argv[6] and sys.argv[7]:
        writeEigenValuesToFile(e2D, e3D, sys.argv[6], sys.argv[7])
    # end if
# end ifmain
