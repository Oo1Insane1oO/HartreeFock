import sys

def degeneracy(i, dim):
    if (dim == 2):
        return 2*(i+1)
    else:
        return (i+1)*(i+2)
    # end ifelse
# end function degeneracy

def createMagic(numParticlesMax, dim):
    """ return list of magic number below incluse numParticlesMax """
    magic = [2]
    n = magic[0]
    while (n < numParticlesMax):
        magic.append(degeneracy(len(magic), dim))
        n += magic[-1]
    # end while
    for i in range(1,len(magic)):
        magic[i] += magic[i-1]
    # end fori
    return magic
# end function createMagic

def makeFiles(omegaList, numParticlesMax, numBasisMax, dim, dName):
    particles = createMagic(numParticlesMax, dim)
    basisSizes = createMagic(numBasisMax, dim)
    for w in omegaList:
        for i in particles:
            for l in basisSizes:
                fname = "w%.2f_N%i_L%i" % (w, i, l)
                with open(dName + "/" + fname + "_hfrun.yaml", "w") as openFile:
                    openFile.write(("omega: %.2f\n"
                                   "numParticles: %i\n"
                                   "dim: %i\n"
                                   "numBasis: %i\n"
                                   "maxIter: 1000\n"
                                   "progress: true\n"
                                   "filename: "+'"'+fname+".yaml"+'"\n') % (w,
                                       i, dim, l))
                # end with open openFile
            # end forl
        # end fori
    # end forw
# end function makeFiles

if __name__ == "__main__":
    try:
        numParticlesMax = int(sys.argv[1])
        numBasisMax = int(sys.argv[2])
        dim = int(sys.argv[3])
        dirName = sys.argv[4]
        omegaList = map(float, sys.argv[5:])
    except IndexError:
        print("USAGE: 'Nmax' 'Lmax' 'dim' 'dir' 'list omega'")
        sys.exit(1)
    # end try-except

    makeFiles(omegaList, numParticlesMax, numBasisMax, dim, dirName)
# end ifmain
