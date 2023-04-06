from pyteomics import mgf, mass
import numpy as np
import sys
import math
from numpy.linalg import norm
from tabulate import tabulate


precursor_mass_tol = 0.1
bin_size = 1.0

class Peptide:
    def __init__(self, seq, mass, tag):
        self.seq = seq
        self.mass = mass
        self.tag = tag

class MS:
    def __init__(self, pep_mass, charge, mass_arr, intensity_arr):
        self.pep_mass = pep_mass
        self.charge = charge
        self.mass_arr = mass_arr
        self.intensity_arr = np.sqrt(intensity_arr)
        self.acc_mass = (pep_mass - 1.0073) * charge

class Spec_pep_pair: 
    def __init__(self, spec, pep, score):
        self.spec = spec
        self.pep = pep
        self.score = score 
        
    
def fragments(peptide, types=('b', 'y'), maxcharge=1):
    for i in range(1, len(peptide)):
        for ion_type in types:
            for charge in range(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    yield mass.fast_mass(
                            peptide[:i], ion_type=ion_type, charge=charge)
                else:
                    yield mass.fast_mass(
                            peptide[i:], ion_type=ion_type, charge=charge)



def parse_pep(target_file, decoy_file):
    pep_lst = []
    f1 = open(target_file, 'r')
    target_lines = f1.readlines()
    for line in target_lines:
        pep, mass = line.strip().split()
        pep_lst.append(Peptide(pep, float(mass), "t"))

    f2 = open(decoy_file, 'r')
    decoy_lines = f2.readlines()
    for line in decoy_lines:
        pep, mass = line.strip().split()
        pep_lst.append(Peptide(pep, float(mass), "d"))
    return pep_lst


def parse_mgf(mgf_file):
    spec_lst = []
    raw_specs = mgf.read(mgf_file)
    for raw_spec in raw_specs:
        spec_lst.append(MS(float(raw_spec["params"]["pepmass"][0]), float(raw_spec["params"]["charge"][0]), raw_spec["m/z array"], raw_spec["intensity array"]))
    return spec_lst



def calc_similarity(spec, pep):
    # Discretize
    pep_frags = list(fragments(pep.seq))
    pep_frags.sort()

    spec_mass_arr = spec.mass_arr
    spec_intensity_arr = spec.intensity_arr

    min_mass = math.floor(min(min(pep_frags), min(spec_mass_arr)))
    max_mass = math.ceil(max(max(pep_frags), max(spec_mass_arr)))

    dis_mass = []
    dis_spec_intensity = []
    dis_pep_intensity = []

    cur_spec_mass_index = 0
    cur_pep_mass_index = 0

    cur_bin_val = min_mass + bin_size

    while cur_bin_val < (max_mass+bin_size):
        dis_mass.append(cur_bin_val)
        cur_bin_intensity = [0]
        if cur_spec_mass_index != len(spec_mass_arr):
            while spec_mass_arr[cur_spec_mass_index] < cur_bin_val:
                cur_bin_intensity.append(spec_intensity_arr[cur_spec_mass_index])
                cur_spec_mass_index =  cur_spec_mass_index + 1
                if cur_spec_mass_index == len(spec_mass_arr):
                    break
        dis_spec_intensity.append(max(cur_bin_intensity))
    
        if cur_pep_mass_index != len(pep_frags):
            if pep_frags[cur_pep_mass_index] < cur_bin_val:
                dis_pep_intensity.append(1.0)
                cur_pep_mass_index = cur_pep_mass_index + 1
            else:
                dis_pep_intensity.append(0.0)
        else:
            dis_pep_intensity.append(0.0)
        
        cur_bin_val = cur_bin_val + bin_size

    spec_intensity_norm = norm(dis_spec_intensity)

    normalized_spec_intensity_arr = dis_spec_intensity / spec_intensity_norm

    similarity_score = np.inner(normalized_spec_intensity_arr, dis_pep_intensity)
    return similarity_score

def calc_max_fdr_line(pair_lst):
    num_target = 0
    num_decoy = 0        

    n = 0

    for i in range(len(pair_lst)):
        if pair_lst[i].pep.tag == "t":
            num_target = num_target + 1
        else:
            num_decoy = num_decoy + 1 

        if num_decoy / num_target > 0.01:
            n = i + 1
            break

    return n

def main(): 
    spec_lst = parse_mgf(sys.argv[1])
    pep_lst = parse_pep(sys.argv[2], sys.argv[3])

    spec_pep_pair_lst = []

    for spec in spec_lst:
        max_pep_score = 0
        pep_index = -1
        for pep in pep_lst:
            if abs(spec.acc_mass - pep.mass) < precursor_mass_tol:
                cur_score = calc_similarity(spec, pep)
                if cur_score > max_pep_score: 
                    max_pep_score = cur_score
                    pep_index = pep_lst.index(pep)

        if pep_index != -1:
            pair = Spec_pep_pair(spec, pep_lst[pep_index], max_pep_score)
            spec_pep_pair_lst.append(pair)

    spec_pep_pair_lst.sort(key=lambda x:x.score, reverse=True)

    n_line = calc_max_fdr_line(spec_pep_pair_lst)
    print(n_line)

    output = []
    for pair in spec_pep_pair_lst:
        line = [spec_lst.index(pair.spec), pair.spec.pep_mass, pair.spec.charge, pair.score, pair.pep.seq]
        output.append(line)

    tabulated_output = tabulate(output, headers=["id", "m/z", "z", "score", "peptide"])

    output_f = open(sys.argv[4], 'w')
    output_f.write(tabulated_output)

            
if __name__ == "__main__":
    main()