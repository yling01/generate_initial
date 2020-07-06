import optparse
from chimScriptMaker import *
from check import *
import os
import shutil
if __name__ == "__main__":
    print("\n!!!NOTE: Test Version, Use With Causion!!!\n")
    parser = optparse.OptionParser()
    parser.add_option("--seq", dest = "seq",
                      default = "NO_INPUT",
                      help = "The sequence to be generated")
    parser.add_option("--cutoffRMSD", dest = "cutoffRMSD",
                      default = "1.0",
                      help = "The RMSD cutoff")
    parser.add_option("--cutoffOmega", dest = "cutoffOmega",
                      default = "150")
    parser.add_option("--cutoffBond", dest = "cutoffBond",
                      default = "1.4")
    parser.add_option("--n", dest = "n",
                      default = "2",
                      help = "Number of structure to be generated")
    parser.add_option('--cyclic', dest = 'cyclic',
                      default = 'True')

    (options, args) = parser.parse_args()

    seq = options.seq
    n = int(options.n)
    cutoffOmega = float(options.cutoffOmega)
    cutoffRMSD = float(options.cutoffRMSD)
    cutoffBond = float(options.cutoffBond)

    if options.cyclic.upper()[0] == "T":
        cyclic = True
    else:
        cyclic = False

    if seq == "NO_INPUT":
        seq = input("\nNo sequence was declared, please enter your sequence:\n")

    assert n >= 1

    counter = attempt = 0

    if os.path.exists(seq):
        shutil.rmtree(seq)
    os.mkdir(seq)


    while counter < n:
        counter += 1
        attempt += 1
        print("\nAttempt %d...\n" % attempt)
        structure = generate_structure(counter, seq)
        while not check_structure(structure, seq, cutoffOmega, cutoffRMSD, cutoffBond, cyclic):
            attempt += 1
            print("\nAttempt %d...\n" % attempt)
            structure = generate_structure(counter, seq)
            if attempt == 50:
                shutil.rmtree(seq)
                print("\nFailed...\n!!!Consider Losen Constraint!!!\nExiting...\n")

                os.remove("chimScript.py")
                os.remove("gly.pdb")
                os.remove("default.profraw")
                shutil.rmtree("__pycache__")

                sys.exit()

        print("\ns%d is successfully generated!\n" % counter)

    print("\nA total of %s structures are generated and are placed in %s/.\n" % (counter, seq))
    print("\n!!!Check chirality, peptide bond, and head-to-tail distance!!!\n")

    os.remove("chimScript.py")
    os.remove("gly.pdb")
    os.remove("default.profraw")
    shutil.rmtree("__pycache__")
