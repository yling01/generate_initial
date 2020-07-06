import random
import os

def write_gly_pdb():
    with open("gly.pdb", "w+") as fo:
        fo.write("ATOM      1  N   GLY A   1      -1.550  -1.067   0.000  1.00  0.00           N\n")
        fo.write("ATOM      2  CA  GLY A   1      -0.101  -1.067   0.000  1.00  0.00           C\n")
        fo.write("ATOM      3  C   GLY A   1       0.422   0.362  -0.000  1.00  0.00           C\n")
        fo.write("ATOM      4  O   GLY A   1       1.315   0.629   0.801  1.00  0.00           O\n")
        fo.write("ATOM      5  OXT GLY A   1      -0.087   1.143  -0.801  1.00  0.00           O\n")
        fo.write("CONECT    3    2    4    5\n")
        fo.write("CONECT    2    1    3\n")
        fo.write("CONECT    1    2\n")
        fo.write("CONECT    4    3\n")
        fo.write("CONECT    5    3\n")
        fo.write("END\n")

def createChimScript(seq, counter):
    length = len(seq)
    phi = generate_angles(length)
    psi = generate_angles(length)
    write_gly_pdb()
    chim = open("chimScript.py", "w+")
    chim.write("import chimera\n")
    chim.write("from chimera import runCommand\n")
    chim.write("runCommand('open gly.pdb')\n")
    chim.write("runCommand('rotation 1 :1@C :1@CA')\n")
    chim.write("runCommand('rotation 1 %s')\n" % phi[0])

    chim.write("runCommand('~select :1@C :1@CA')\n")
    for index in range(1, length):
        chim.write("runCommand('addaa gly,%d,%s,%s :%d.a')\n" % ((index + 1),
                                                              phi[index],
                                                              psi[index],
                                                              index))

    chim.write("runCommand('delete :%s@OXT')\n" % length)
    chim.write("runCommand('bond :%s@C :1@N')\n" % length)

    chim.write("runCommand('minimize nogui True, nsteps 1000')\n")

    for index, amino in enumerate(seq):
        toInvert = False
        isPro = False
        if amino.islower():
            toInvert = True
            amino = amino.upper()
        if amino == "P":
            amino = "A"
            isPro = True
        chim.write("runCommand('swapaa %s #:%s.a')\n" % (oneToThree(amino), index + 1))
        if toInvert:
            chim.write("runCommand('invert :%s@ca')\n" % (index + 1))
            if amino == "V":
                chim.write("runCommand('invert :%s@cb')\n" % (index + 1))
        if isPro:
            chim.write("runCommand('swapaa pro #:%s.a')\n" % (index + 1))

    chim.write("runCommand('minimize nogui True, nsteps 5000')\n")
    chim.write("runCommand('select element.H')\n")
    chim.write("runCommand('~select :1.A@H')\n")
    chim.write("runCommand('delete sel')\n")

    file_name = "%s/s%d.pdb" % (seq, counter)
    chim.write("runCommand('write #0 %s')\n" % file_name)

    chim.close()

    return file_name

def oneToThree(one):
    aa_dic = {"A": "ala",
              "C": "cys",
              "D": "asp",
              "E": "glu",
              "F": "phe",
              "G": "gly",
              "H": "his",
              "I": "ile",
              "K": "lys",
              "L": "leu",
              "M": "met",
              "N": "asn",
              "P": "pro",
              "Q": "gln",
              "R": "arg",
              "S": "ser",
              "T": "thr",
              "V": "val",
              "Y": "tyr"}

    return aa_dic[one]

def generate_angles(length):
    return [random.randint(-180, 180) for _ in range(length)]

def generate_structure(counter, seq):
    file_name = createChimScript(seq, counter)
    os.system("/Applications/Chimera.app/Contents/MacOS/chimera --script chimScript.py --nogui --silent &> /dev/null")
    return file_name
