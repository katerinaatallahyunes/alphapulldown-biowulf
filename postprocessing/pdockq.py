import argparse
import sys
import os
import numpy as np
import pandas as pd
from collections import defaultdict
import pdb
import glob

parser = argparse.ArgumentParser(description='Calculate a predicted DockQ score for a predicted structure.')
parser.add_argument('--pdbfile', nargs=1, type=str, default=sys.stdin, help='Path to pdbfile to be scored. Note that this file needs to contain at least three chains. The B-factor column is assumed to contain the plDDT score from AlphaFold.')

#####################FUNCTIONS#########################
def parse_atm_record(line):
    '''Get the atm record
    '''
    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['atm_alt'] = line[17]
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])

    return record

def read_pdb(pdbfile):
    '''Read a pdb file predicted with AF and rewritten to contain all chains
    '''

    chain_coords, chain_plddt = {}, {}
    with open(pdbfile, 'r') as file:
        for line in file:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            # Get CB - CA for GLY
            if record['atm_name'] == 'CB' or (record['atm_name'] == 'CA' and record['res_name'] == 'GLY'):
                if record['chain'] in [*chain_coords.keys()]:
                    chain_coords[record['chain']].append([record['x'], record['y'], record['z']])
                    chain_plddt[record['chain']].append(record['B'])
                else:
                    chain_coords[record['chain']] = [[record['x'], record['y'], record['z']]]
                    chain_plddt[record['chain']] = [record['B']]

    # Convert to arrays
    for chain in chain_coords:
        chain_coords[chain] = np.array(chain_coords[chain])
        chain_plddt[chain] = np.array(chain_plddt[chain])

    return chain_coords, chain_plddt

def calc_pdockq(chain_coords, chain_plddt, t):
    '''Calculate the pDockQ scores
    pdockQ = L / (1 + np.exp(-k*(x-x0)))+b
    L= 0.724 x0= 152.611 k= 0.052 and b= 0.018
    '''

    chains = [*chain_coords.keys()]
    num_chains = len(chains)
    if num_chains != 3:
        print('This script requires PDB files with exactly three chains.')
        return None, None

    coords = [chain_coords[ch] for ch in chains]
    plddts = [chain_plddt[ch] for ch in chains]

    # Calculate distances between all atoms
    a_min_b = coords[0][:, np.newaxis, :] - coords[1][np.newaxis, :, :]
    dists_ab = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
    a_min_c = coords[0][:, np.newaxis, :] - coords[2][np.newaxis, :, :]
    dists_ac = np.sqrt(np.sum(a_min_c.T ** 2, axis=0)).T
    b_min_c = coords[1][:, np.newaxis, :] - coords[2][np.newaxis, :, :]
    dists_bc = np.sqrt(np.sum(b_min_c.T ** 2, axis=0)).T

    contact_dists = [dists_ab, dists_ac, dists_bc]
    contacts = [np.argwhere(d <= t) for d in contact_dists]
    
    if all(len(contact) == 0 for contact in contacts):
        pdockq = 0
        ppv = 0
    else:
        all_if_plddt = []
        for chain_idx in range(len(plddts)):
            unique_contact_indices = np.unique(contacts[chain_idx][:, 0])
            if len(unique_contact_indices) > 0:
                # Ensure that indices are within bounds
                valid_indices = unique_contact_indices[unique_contact_indices < len(plddts[chain_idx])]
                if len(valid_indices) > 0:
                    all_if_plddt.append(plddts[chain_idx][valid_indices])
    
        if len(all_if_plddt) > 0:
            avg_if_plddt = np.average(np.concatenate(all_if_plddt))
    
            n_if_contacts = sum(len(contact) for contact in contacts)
            x = avg_if_plddt * np.log10(n_if_contacts)
            pdockq = 0.724 / (1 + np.exp(-0.052 * (x - 152.611))) + 0.018
        else:
            pdockq = 0
    
        PPV = np.array([0.98128027, 0.96322524, 0.95333044, 0.9400192,
                        0.93172991, 0.92420274, 0.91629946, 0.90952562, 0.90043139,
                        0.8919553, 0.88570037, 0.87822061, 0.87116417, 0.86040801,
                        0.85453785, 0.84294946, 0.83367787, 0.82238224, 0.81190228,
                        0.80223507, 0.78549007, 0.77766077, 0.75941223, 0.74006263,
                        0.73044282, 0.71391784, 0.70615739, 0.68635536, 0.66728511,
                        0.63555449, 0.55890174])
    
        pdockq_thresholds = np.array([0.67333079, 0.65666073, 0.63254566, 0.62604391,
                                    0.60150931, 0.58313803, 0.5647381, 0.54122438, 0.52314392,
                                    0.49659878, 0.4774676, 0.44661346, 0.42628389, 0.39990988,
                                    0.38479715, 0.3649393, 0.34526004, 0.3262589, 0.31475668,
                                    0.29750023, 0.26673725, 0.24561247, 0.21882689, 0.19651314,
                                    0.17606258, 0.15398168, 0.13927677, 0.12024131, 0.09996019,
                                    0.06968505, 0.02946438])
        
        inds = np.argwhere(pdockq_thresholds >= pdockq)
        if len(inds) > 0:
            ppv = PPV[inds[-1]][0]
        else:
            ppv = PPV[0]
    

    return pdockq, ppv


#################MAIN####################

# Get a list of pdb files using glob
pdb_files = glob.glob('*_model_*.pdb')

# Create a list to store the scores
scores = []

# Loop over the pdb files
for pdb_file in pdb_files:
    # Read chains
    chain_coords, chain_plddt = read_pdb(pdb_file)

    # Check chains
    if len(chain_coords.keys()) < 2:
        print('Only one chain in pdb file:', pdb_file)
        continue

    # Calculate pdockq
    t = 8  # Distance threshold, set to 8 Å
    pdockq, ppv = calc_pdockq(chain_coords, chain_plddt, t)
    
    # Store the scores in a dictionary
    score_dict = {
        'pdb_file': pdb_file,
        'pdockq': pdockq,
        'ppv': ppv
    }
    scores.append(score_dict)

# Create a DataFrame from the scores list
scores_df = pd.DataFrame(scores)

# Sort the scores based on the pdb_file name
scores_df = scores_df.sort_values(by='pdb_file')

# Save the scores to a CSV file
scores_df.to_csv('pdockq.csv', index=False)
