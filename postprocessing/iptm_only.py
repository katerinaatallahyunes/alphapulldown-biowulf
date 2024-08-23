#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Katerina Atallah-Yunes

This script was adapted from Chop Yan Lee's original code. The implemented changes including removing additional data output that was 
not needed for my analysis pipeline.

The original code source by Chop Yan Lee: https://github.com/KatjaLuckLab/AlphaFold_manuscript/blob/main/scripts/calculate_template_independent_metrics.py
"""

# This script contains generic functions that extract and manipulate metrics and information that can be obtained from predicted models without the need of a template.
# Author: Chop Yan Lee
# pDockQ code source: https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py
# iPAE code source: https://github.com/fteufel/alphafold-peptide-receptors/blob/main/qc_metrics.py

from pymol import cmd
import numpy as np
import pandas as pd
import json, os, pickle, argparse, sys
from collections import defaultdict

class Prediction_folder:
    """Class that stores prediction folder information"""
    def __init__(self,prediction_folder,num_model=5,project_name=None):
        """Initialize an instance of Prediction

        Args:
            prediction_folder (str): absolute path to the prediction folder
        """
        self.prediction_folder = prediction_folder
        self.num_model = num_model
        self.path_to_prediction_folder = os.path.split(self.prediction_folder)[0]
        self.prediction_name = os.path.split(self.prediction_folder)[1]
        self.rank_to_model = {}
        self.model_confidences = {}
        self.fasta_sequence_dict = {'A':'','B':''}
        # instantiate the amount of Predicted_model according to the number of models given as argument, otherwise 5
        self.model_instances = {}
        if project_name is not None:
            self.project_name = project_name
        # need an attribute to annotate if a prediction folder has been successfully predicted without internal AlphaFold error
        self.predicted = True

    def parse_ranking_debug_file(self):
        """Read the ranking_debug_file and save relevant information into attribute of self
        """
        if not os.path.exists(os.path.join(self.prediction_folder,'ranking_debug.json')):
            self.predicted = False
            return
        else:
            with open(os.path.join(self.prediction_folder,'ranking_debug.json'), 'r') as f:
                data = json.load(f)
            self.rank_to_model = {f'ranked_{i}':model for i, model in enumerate(data.get("order"))}
            sorted_model_confidence = sorted(data.get("iptm+ptm").values(),reverse=True)
            self.model_confidences = {f'ranked_{i}':float(confidence) for i, confidence in enumerate(sorted_model_confidence)}
        
    def parse_prediction_fasta_file(self):
        """Read the fasta file of the prediction to retrieve information on chain and sequence identity
        """
        fasta_path = f'{self.prediction_folder}.fasta'
        with open(fasta_path, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip() != '']
        chain_id = 0
        for line in lines:
            if line[0] == '>':
                chain = list(self.fasta_sequence_dict)[chain_id]
                chain_id += 1
                continue
            self.fasta_sequence_dict[chain] += line

    def instantiate_predicted_model(self):
        """Initialize the amount of Predicted_model instance according to the number of model specified and save it in the dict self.model_instances
        """
        self.model_instances = {f'ranked_{i}':Predicted_model(f'ranked_{i}') for i in range(self.num_model)}

    def assign_model_info(self):
        """Assign information stored in the prediction folder to their corresponding predicted model
        """
        for model_id, model_inst in self.model_instances.items():
            model_inst.model_confidence = self.model_confidences.get(model_id)
            model_inst.multimer_model = self.rank_to_model.get(model_id)
            model_inst.path_to_model = self.prediction_folder

    def process_all_models(self):
        """Use the instances of Predicted_model and run the wrapper function Predicted_model.get_model_independent_metrics function on themselves
        """
        self.parse_ranking_debug_file()
        self.parse_prediction_fasta_file()
        if self.predicted:
            self.instantiate_predicted_model()
            self.assign_model_info()
            for model_id, model_inst in self.model_instances.items():
                model_inst.get_model_independent_metrics()
    
    def write_out_calculated_metrics(self, project_name=None):
        """
        Write out the information that has been processed for every predicted model.
    
        Args:
            project_name (str): a project name given to model contacts dataframe as a key identifier
    
        Returns:
            template_indep_info.tsv: A tsv file with the calculated template independent metrics
        """
        metrics_out_path = os.path.join(self.path_to_prediction_folder, 'template_indep_info.tsv')
    
        # Define the columns and their data types for the DataFrame
        metrics_columns_dtype = {'project_name': str, 'prediction_name': str, 'chain_A_length': int, 'chain_B_length': int, 'model_id': str, 'model_confidence': float}
    
        # Check if template_indep_info.tsv already exists
        if os.path.exists(metrics_out_path):
            metrics_df = pd.read_csv(metrics_out_path, sep='\t', index_col=0)
            metrics_df.reset_index(drop=True, inplace=True)
        else:
            metrics_df = pd.DataFrame(columns=metrics_columns_dtype.keys())
            metrics_df = metrics_df.astype(dtype=metrics_columns_dtype)
    
        if self.project_name is not None:
            common_info = [self.project_name]
        else:
            common_info = []
        common_info += [self.prediction_name, len(self.fasta_sequence_dict.get('A')), len(self.fasta_sequence_dict.get('B'))]
    
        # check if the prediction folder has been predicted successfully without internal error from AlphaFold
        if not self.predicted:
            row = common_info + ['Prediction failed'] + [None] * 3 + [None] * (len(metrics_columns_dtype) - 6)  # Adjust the None values based on column count
            metrics_df.loc[len(metrics_df)] = row
        else:
            # insert metric info in a row-wise manner
            for model_id, model_inst in self.model_instances.items():
                row = common_info + [model_id, model_inst.model_confidence]  # Add other model-specific metrics as needed
                metrics_df.loc[len(metrics_df)] = row
    
        # Filter the DataFrame based on specific criteria (example: model_confidence >= 0.5)
        filtered_df = metrics_df[metrics_df['model_confidence'] >= 0.5]
    
        # Write out the filtered DataFrame to a new file
        filtered_out_path = os.path.join(self.path_to_prediction_folder, 'filtered_template_indep_info.tsv')
        filtered_df.to_csv(filtered_out_path, sep='\t', index=False)
    
        # Write out the original metrics DataFrame
        metrics_df.to_csv(metrics_out_path, sep='\t')
        print(f'Calculated metrics saved in {metrics_out_path}!')
        print(f'Filtered metrics saved in {filtered_out_path}!')
    



class Predicted_model:
    """Class that stores predicted model"""
    def __init__(self,predicted_model):
        """Initialize an instance of Predicted_model
        
        Args:
            predicted_model (str): name of the predicted model like ranked_0
        """
        self.predicted_model = predicted_model
        self.path_to_model = None
        self.multimer_model = None
        self.chain_coords = None
        self.chain_plddt = None
        self.pickle_data = None
        self.model_confidence = None

    def check_chain_id(self):
        """Some models have their chain ids start from B instead of A. As the code requires the chain ids to be consistent (start from chain A), this function checks and rename the chain ids if necessary
        """
        model_path = os.path.join(self.path_to_model,f'{self.predicted_model}.pdb')
        cmd.load(model_path)
        chains = cmd.get_chains(f'{self.predicted_model}')
        if 'C' in chains:
            # change the chain id into a temporary arbitrary name
            cmd.alter(f'{self.predicted_model} and chain B', 'chain="tempA"')
            cmd.sort()
            cmd.alter(f'{self.predicted_model} and chain C', 'chain="tempB"')
            cmd.sort()
            # change the chain id into A and B
            cmd.alter(f'{self.predicted_model} and chain tempA', 'chain="A"')
            cmd.sort()
            cmd.alter(f'{self.predicted_model} and chain tempB', 'chain="B"')
            cmd.sort()
            # save and overwrite the predicted model
            cmd.save(model_path,f'{self.predicted_model}')
        cmd.reinitialize()

    def read_pickle(self):
        """Read in the pickle file of multimer model

        Returns:
            self.pickle_data (dict): Pickle data of multimer model
        """
        multimer_model_pickle = os.path.join(self.path_to_model,f'result_{self.multimer_model}.pkl')
        with open(multimer_model_pickle, 'rb') as f:
            self.pickle_data = pickle.load(f)

    def parse_atm_record(self,line):
        """Get the atm record from pdb file

        Returns:
            record (dict): Dict of parsed pdb information from .pdb file
        """
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
    
    def read_pdb(self):
        """Read a pdb file predicted with AF and rewritten to conatin all chains

        Returns:
            self.chain_coords (dict): Dict of chain coordination (x,y,z)
            self.chain_plddt (dict): Dict of chain id as key and plddt array as value
        """

        chain_coords, chain_plddt = {}, {}
        model_path = os.path.join(self.path_to_model,f'{self.predicted_model}.pdb')

        with open(model_path, 'r') as file:
            for line in file:
                if not line.startswith('ATOM'):
                    continue
                record = self.parse_atm_record(line)
                #Get CB - CA for GLY
                if record['atm_name']=='CB' or (record['atm_name']=='CA' and record['res_name']=='GLY'):
                    if record['chain'] in [*chain_coords.keys()]:
                        chain_coords[record['chain']].append([record['x'],record['y'],record['z']])
                        chain_plddt[record['chain']].append(record['B'])
                    else:
                        chain_coords[record['chain']] = [[record['x'],record['y'],record['z']]]
                        chain_plddt[record['chain']] = [record['B']]

        #Convert to arrays
        for chain in chain_coords:
            chain_coords[chain] = np.array(chain_coords[chain])
            chain_plddt[chain] = np.array(chain_plddt[chain])

        self.chain_coords = chain_coords
        self.chain_plddt = chain_plddt

    def parse_ptm_iptm(self):
        """Parse the ptm and iptm of a predicted model by using the pickle file of the multimer model where the ptm and iptm can be found
        
        Returns:
            ptm (float): the parsed ptm, saved as attribute of self
            iptm (float): the parsed iptm, saved as attribute of self
        """
        self.ptm = float(self.pickle_data['ptm'])
        self.iptm = float(self.pickle_data['iptm'])

    def get_model_independent_metrics(self):
        """Wraps all the functions together to process a predicted model

        Returns:
            None
        """
        # self.parse_ptm_iptm() # skipped for now because the pickle file has JAX dependency and I am not sure what to do with it
        self.check_chain_id()
        if 'multimer_v2' in self.multimer_model:
            if os.path.exists(os.path.join(self.path_to_model,f'result_{self.multimer_model}.pkl')):
                self.read_pickle()
                self.calculate_iPAE()
        print(f'{os.path.join(self.path_to_model,self.predicted_model)} processed!')

def main():
    """Parse arguments and wraps all functions into main for executing the program in such a way that it can handle multiple run ids given to it
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-run_ids', type=str, help='Run IDs for metrics calculation', dest='run_ids')
    parser.add_argument('-path_to_run', type=str, help='Either provide a path to a folder where multiple runs of AlphaFold predictions are contained and specify the run_ids to be processed or use -path_to_prediction to specify a folder that you want to process, include "/" at the end', dest='path_to_run')
    parser.add_argument('-path_to_prediction', type=str, help='Path to the prediction folder "/" at the end', dest='path_to_prediction')
    parser.add_argument('-project_name', type=str, help='Optional name for the project', dest='project_name')
    parser.add_argument('-skip_write_out_contacts', action='store_true', help='Exclude writing out  found in predicted models', dest='skip_write_out_contacts')
    args = parser.parse_args()
    run_ids = vars(args)['run_ids']
    path_to_run = vars(args)['path_to_run']
    path_to_prediction = vars(args)['path_to_prediction']
    project_name = vars(args)['project_name']
    skip_contacts = vars(args)['skip_write_out_contacts']

    # a list to contains already processed files
    calculated_files = []

    # check which argument, -path_to_run or -path_to_prediction, is provided
    if (path_to_run is None) and (path_to_prediction is None):
        print('Please provide either -path_to_run or -path_to_prediction and try again!')
        sys.exit()
    elif path_to_prediction is not None:
        if os.path.exists(f'{path_to_prediction}template_indep_info.tsv'):
            temp = pd.read_csv(f'{path_to_prediction}template_indep_info.tsv',sep='\t',index_col=0)
            calculated_files = temp['prediction_name'].unique()
        for file in os.listdir(path_to_prediction):
            if file in calculated_files:
                # print(file)
                continue
            file_abs = os.path.join(path_to_prediction,file)
            if os.path.isdir(file_abs):
                folder = Prediction_folder(file_abs,num_model=5,project_name=project_name)
                folder.process_all_models()
                folder.write_out_calculated_metrics()
    else:
        for run_id in run_ids.split(','):
            if os.path.exists(f'{path_to_run}run{run_id}/template_indep_info.tsv'):
                temp = pd.read_csv(f'{path_to_run}run{run_id}/template_indep_info.tsv',sep='\t',index_col=0)
                calculated_files = temp['prediction_name'].unique()
            run_path = f'{path_to_run}run{run_id}'
            for file in os.listdir(run_path):
                if file in calculated_files:
                    # print(file)
                    continue
                file_abs = os.path.join(run_path,file)
                if os.path.isdir(file_abs):
                    if not os.path.exists(os.path.join(run_path,f"{file}.fasta")):
                        print(f"Skipping the folder named {file}")
                        continue
                    folder = Prediction_folder(file_abs,num_model=5,project_name=project_name)
                    folder.process_all_models()
                    folder.write_out_calculated_metrics()
                    if skip_contacts:
                        continue
                    folder.write_out_contacts()

if __name__ == '__main__':
    main()
