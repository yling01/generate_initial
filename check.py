import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.dihedrals import Dihedral
import os
import sys
import shutil
def get_ref(seq, file_name, result_path):
    seq = "/".join([result_path, seq])
    return ["/".join([seq, i]) for i in os.listdir(seq) if i[-3:] == "pdb" and "/".join([seq, i]) != file_name]

def calculate_omega(u, start, stop, cyclic):
    u_protein = u.select_atoms('protein')
    num_res = len(u_protein.residues)
    ags_omega = [res.omega_selection() for res in u_protein.residues[:-1]]
    if cyclic:
        ags_omega += [u_protein.select_atoms("resid %d and name CA" % num_res, \
                                             "resid %d and name C" % num_res, \
                                             "resid 1 and name N", \
                                             "resid 1 and name CA")]
    R_omega = Dihedral(ags_omega).run(start=start, stop=stop)
    return R_omega.angles

def pick_out_trans(omega_angles, cutoff):
    res = np.ones(len(omega_angles), dtype=bool)
    for i, omega_angle in enumerate(omega_angles):
        if np.all(np.absolute(omega_angle) > cutoff):
            continue
        else:
            res[i] = False

    return res

def pick_out_chirality(chirality_set, seq):
    assert len(seq) == len(chirality_set[0])
    res = np.ones(len(chirality_set), dtype=bool)
    for j, frame in enumerate(chirality_set):
        for i, chirality in enumerate(frame):
            if chirality == "L" and seq[i].isupper():
                continue
            elif chirality == "D" and seq[i].islower():
                continue
            else:
                res[j] = False
    return res

def calculate_chirality(u, start, stop):
    u_protein = u.select_atoms('protein')
    num_frame = len(u.trajectory)
    sequence = [str(i)[9:12] for i in u_protein.residues]
    chirality = np.ones((num_frame, 1))
    for i, aa in enumerate(sequence, start=1):
        if aa == "GLY":
            chirality_temp = np.empty((num_frame, 1), dtype=str)
            chirality_temp.fill("L")
            chirality = np.hstack((chirality, chirality_temp))
        else:
            ags_improper = u_protein.select_atoms("resid %d and name CA" % i, \
                                                  "resid %d and name N" % i, \
                                                  "resid %d and name C" % i, \
                                                  "resid %d and name CB" % i)
            improper_angles = Dihedral([ags_improper]).run(start=start, stop=stop).angles
            chirality_temp = np.empty((num_frame, 1), dtype=str)
            for i, angle in enumerate(improper_angles):
                if angle == 0:
                    sys.exit("\nExiting...Ambigious Improper Dihedral Angle...\n")
                elif angle > 0:
                    chirality_temp[i] = "L"
                else:
                    chirality_temp[i] = "D"
            chirality = np.hstack((chirality, chirality_temp))
    return chirality[:, 1:]

def check_rmsd(u, ref_structures, cutoff_rmsd):
    if not ref_structures:
        return True
    atomGroup = "backbone and not name O"
    for structure in ref_structures:
        u_ref = mda.Universe(structure)
        rmsValue = rms.rmsd(u.select_atoms(atomGroup).positions, \
                            u_ref.select_atoms(atomGroup).positions, \
                            superposition=True)
        if rmsValue < cutoff_rmsd:
            return False
    return True

def check_cyclization(u, cutoff_bond, cyclic):
    if cyclic:
        return True
    nTerm = u.select_atoms("resid 1 and name N")
    cTerm = u.select_atoms("resid %d and name C" % len(u.residues))
    bond_distance = mda.analysis.distances.dist(nTerm, cTerm)[2][0]
    return bond_distance <= cutoff_bond

def recycle():
    trash_path = "Discarded"
    if os.path.exists(trash_path):
        trash_files = [int(i.split("_")[0]) for i in os.listdir(trash_path) if i[-3:] == "pdb"]
        trash_files.sort()
        return trash_files[-1] if trash_files else 0
    else:
        os.mkdir(trash_path)
        return 0

def check_structure(file_name, seq, cutoff_omega, cutoff_rmsd, cutoff_bond, cyclic, result_path):
    u = mda.Universe(file_name)
    omega_angles = calculate_omega(u, 0, 1, cyclic)
    chirality = calculate_chirality(u, 0, 1)
    omega_mask = pick_out_trans(omega_angles, cutoff_omega)
    chirality_mask = pick_out_chirality(chirality, seq)
    rmsd_mask = check_rmsd(u, get_ref(seq, file_name, result_path), cutoff_rmsd)
    cyclization_mask = check_cyclization(u, cutoff_bond, cyclic)
    if not (np.all(chirality_mask & omega_mask) and rmsd_mask and cyclization_mask):
        if not np.all(chirality_mask & omega_mask):
            print("chirality or trans/cis problem")
        if not rmsd_mask:
            print("rmsd problem")
        if not cyclization_mask:
            print("cyclization problem")
        file_dest = "Discarded/%d_%s.pdb" % (recycle() + 1, seq)
        print("\nFile moved to %s\n" % file_dest)

        shutil.move(file_name, file_dest)
        return False
    return True
