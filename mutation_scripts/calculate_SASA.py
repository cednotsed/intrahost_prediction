import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.Data.PDBData import residue_sasa_scales

# Create three letter dictionary
low_dict = protein_letters_3to1
prot_dict = dict()
for key in low_dict.keys():
    prot_dict[key.upper()] = low_dict[key]

# Get maxASA
max_dict = residue_sasa_scales['Wilke']

# Read PDB file
p = PDBParser(QUIET=1)
struct = p.get_structure("sc2", "/mnt/c/git_repos/intrahost_prediction/data/external_datasets/6vxx.filt.pdb")

# Calculate SASA
sr = ShrakeRupley()
sr.compute(struct, level="R")

# Save SASA
my_list = []
for chain in struct[0]:
    for res in chain:
        res_name = res.get_resname()
        if res_name != "NAG":
            codon_number = res.get_id()[1]
            SASA = res.sasa
            max_ASA = max_dict[res_name]
            chain = res.get_full_id()[2]

            my_list.append((prot_dict[res_name], codon_number, SASA, max_ASA, chain))

df = pd.DataFrame(my_list, columns = ['ref_AA', 'codon_number', 'SASA', 'maxASA', 'chain'])
df.to_csv('/mnt/c/git_repos/intrahost_prediction/data/external_datasets/6vxx.filt.SASA.csv',
          sep = ',',
          index = False)
