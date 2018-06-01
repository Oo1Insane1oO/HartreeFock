import sys

def makeFiles(omegaList, particlesList, basisList, dim, dName, dataDirName):
    for w in omegaList:
        for p,i in enumerate(particlesList):
            for j,l in enumerate(basisList):
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
        particlesList = map(int, sys.argv[1].strip('[]').split(','))
        basisList = map(int, sys.argv[2].strip('[]').split(','))
        dirName = sys.argv[3]
        dataDirName = sys.argv[4]
        omegaList = map(float, sys.argv[5:])
    except IndexError:
        print("USAGE: 'Nmax' 'Lmax' 'dim' 'dir' 'datadir' 'list omega'")
        sys.exit(1)
    # end try-except

    makeFiles(omegaList, particlesList, basisList, 3, dirName, dataDirName)
# end ifmain
