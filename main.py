import optparse
from chimScriptMaker import *
from check import *
import os
import shutil
from pathlib import Path

def run(seq, cutoffRMSD, cutoffOmega, cutoffBond, n, cyclic, result_path):

    print("=" * 80)
    print("\nMaking structures for %s\n" % seq)
    os.mkdir("%s/%s" % (result_path, seq))
    counter = attempt = 0
    while counter < n:
        counter += 1
        attempt += 1
        print("\nAttempt %d...\n" % attempt)
        structure = generate_structure(counter, seq, result_path)
        while not check_structure(structure, seq, cutoffOmega, cutoffRMSD, cutoffBond, cyclic, result_path):
            attempt += 1
            print("\nAttempt %d...\n" % attempt)
            structure = generate_structure(counter, seq, result_path)
            if attempt == 50:
                shutil.rmtree(seq)
                print("\nFailed...\n!!!Consider Losen Constraint!!!\nExiting...\n")

                os.remove("chimScript.py")
                os.remove("gly.pdb")
                os.remove("default.profraw")
                shutil.rmtree("__pycache__")

                sys.exit()

        print("s%d is successfully generated!\n" % counter)

    print("\nA total of %s structures are generated and are placed in %s/%s/.\n" % (counter, result_path, seq))
    print("\n!!!Check chirality, peptide bond, and head-to-tail distance!!!\n")


if __name__ == "__main__":
    print("\n!!!NOTE: Test Version, Use With Causion!!!\n")
    print("\n\
*\t                                       /;    ;\\\n\
*\t                                   __  \\____//\n\
*\t                                  /{_\_/   `'\____\n\
*\t                                  \___   (o)  (o  }\n\
*\t       _____________________________/          :--'\n\
*\t   ,-,'`@@@@@@@@       @@@@@@         \_    `__\\\n\
*\t  ;:(  @@@@@@@@@        @@@             \___(o'o)\n\
*\t  :: )  @@@@          @@@@@@        ,'@@(  `===='\n\
*\t  :: : @@@@@:          @@@@         `@@@:\n\
*\t  :: \  @@@@@:       @@@@@@@)    (  '@@@'\n\
*\t  ;; /\      /`,    @@@@@@@@@\   :@@@@@)\n\
*\t  ::/  )    {_----------------:  :~`,~~;\n\
*\t ;;'`; :   )                  :  / `; ;\n\
*\t;;;; : :   ;                  :  ;  ; :\n\
*\t`'`' / :  :                   :  :  : :\n\
*\t    )_ \__;      \";\"          :_ ;  \_\       `,','\n\
*\t    :__\  \    * `,'*         \  \  :  \   *  8`;'*  *\n\
*\t        `^'     \ :/           `^'  `-^-'   \\v/ :  \/\n\
*\t          \n")
    parser = optparse.OptionParser()
    parser.add_option("--seq", dest = "seq",
                      default = "NO_INPUT",
                      help = "The sequence to be generated or the file that contains sequences")
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
        seq = input("\nNo sequence was declared, please enter your sequence or the file that contains sequences:\n")

    assert n >= 1

    result_path = "Result"

    if os.path.exists(result_path):
        shutil.rmtree(result_path)
    os.mkdir(result_path)

    my_file = Path(seq)
    if my_file.is_file():
        fo = open(seq)
        seq_list = map(str.strip, list(fo))
        fo.close()
    else:
        seq_list = [seq]

    for seq in seq_list:
        run(seq, cutoffRMSD, cutoffOmega, cutoffBond, n, cyclic, result_path)

    os.remove("chimScript.py")
    os.remove("gly.pdb")
    os.remove("default.profraw")
    shutil.rmtree("__pycache__")

    print("\n\
   ***\n\
  ** **\n\
 **   **\n\
 **   **         ****\n\
 **   **       **   ****\n\
 **  **       *   **   **\n\
  **  *      *  **  ***  **\n\
   **  *    *  **     **  *\n\
    ** **  ** **        **\n\
    **   **  **\n\
   *           *\n\
  *             *\n\
 *    0     0    *\n\
 *   /   @   \   *\n\
 *   \__/ \__/   *\n\
   *     W     *\n\
     **     **\n\
       *****\n")
    print("=" * 80)

