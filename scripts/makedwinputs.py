import sys

def makeFiles(omegaList, numParticlesMax, numBasisMax, dim, dName, dataDirName):
    particles = (range(numParticlesMax+1)[::2])[1:]
    basisSizes = (range(numBasisMax+1)[::2])[1:]
    for w in omegaList:
        for p,i in enumerate(particles):
            for j,l in enumerate(basisSizes):
                if l >= i:
                    fname = "w%.2f_D%i_N%i_L%i" % (w, dim, i, l)
                    with open(dName + "/" + fname + "_DWhfrun.yaml", "w") as openFile:
                        openFile.write(("R: %.2f\n"
                                       "numParticles: %i\n"
                                       "dim: %i\n"
                                       "numBasis: %i\n"
                                       "maxIter: 1000\n"
                                       "progress: true\n"
                                       "filename: "+'"' + dataDirName + '/' +
                                       fname + ".yaml" + '"\n') % (w,i,dim,l))
                    # end with open openFile
                # end if
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
        dataDirName = sys.argv[5]
        omegaList = map(float, sys.argv[6:])
    except IndexError:
        print("USAGE: 'Nmax' 'Lmax' 'dim' 'dir' 'datadir' 'list omega'")
        sys.exit(1)
    # end try-except

    makeFiles(omegaList, numParticlesMax, numBasisMax, dim, dirName,
            dataDirName)
# end ifmain
