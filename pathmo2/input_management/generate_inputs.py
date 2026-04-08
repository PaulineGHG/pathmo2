import csv
import json
import os.path
from utils import mol_to_asp, rxn_to_asp

from rdkit.Chem import AllChem
from rdkit.Chem import Draw

METACYC_MAPPING_FILE = '../data/atom-mappings-smarts.json'
METACYC_REACTIONS_FILE = '../data/reactions.json'
METACYC_COMPOUNDS_FILE = '../data/compounds.json'


def generate_input(run_name, output_path, source, target, metacyc_ref_reactions):
    run_path = os.path.join(output_path, run_name)
    if not os.path.exists(run_path):
        os.mkdir(run_path)
        os.mkdir(os.path.join(run_path, 'Inputs'))
        os.mkdir(os.path.join(run_path, 'Outputs'))

    reactions_file = os.path.join(run_path, 'Inputs', 'Reactions_references.tsv')

    if metacyc_ref_reactions:

        with open('input.lp', 'w') as flp:
            rxn_data = {}
            with open(METACYC_REACTIONS_FILE, 'r') as f:
                r_f = json.load(f)
                for rxn in metacyc_ref_reactions:
                    rxn_data[rxn] = r_f[rxn]

            all_cpd = set()
            for rxn, data in rxn_data.items():
                for c in data['left']:
                    all_cpd.add(c)
                for c in data['right']:
                    all_cpd.add(c)

            cpd_data = {}
            with open(METACYC_COMPOUNDS_FILE, 'r') as f:
                c_f = json.load(f)
                for cpd in all_cpd:
                    cpd_data[cpd] = c_f[cpd]

            write_chemicals(source, target, cpd_data, flp)
            write_reactions(rxn_data, flp)

            mc_ref_mappings = {}
            with open(METACYC_MAPPING_FILE, 'r') as f:
                mc_map = json.load(f)
                for r in metacyc_ref_reactions:
                    mc_ref_mappings[r] = mc_map[r]

            flp.write(f'%*\nMAPPINGS\n{100 * "="}\n*%\n')
            for rxn_id, mapping in mc_ref_mappings.items():
                rxn = AllChem.ReactionFromSmarts(mapping)
                reactants = rxn.GetReactants()
                products = rxn.GetProducts()
                write_mappings(reactants, products, rxn_id, flp)

                draw_rxn(rxn, os.path.join(run_path, 'Outputs', rxn_id + '.png'))


def draw_rxn(rxn, output):
    d2d = Draw.MolDraw2DCairo(1600, 600)
    d2d.DrawReaction(rxn, highlightByReactant=False)
    png = d2d.GetDrawingText()
    open(output, 'wb+').write(png)


def write_chemicals(source, target, rxn_chemicals_lst, lp_f):
    lp_f.write(f'%*\nCHEMICALS\n{100 * "="}\n*%\n')
    for chem_id, smile in {**source, **target, **rxn_chemicals_lst}.items():
        lp_f.write(f'\n% {chem_id}\n')
        asp_atoms = mol_to_asp(mol_name=chem_id, mol_code=smile, encoding='SMILES')
        for atom in asp_atoms:
            lp_f.write(f'{atom}\n')

    lp_f.write(f'\n%*\nSOURCE - GOAL\n{100*"="}\n*%\n')
    source = list(source)[0]
    lp_f.write(f'\n% SOURCE\nsource("{source}").\n')
    lp_f.write(f'\n% GOAL\ngoal(pathway("{source}","{list(target)[0]}")).\n')


def write_reactions(reactions_data, lp_f):
    lp_f.write(f'\n%*\nREACTIONS\n{100 * "="}\n*%\n')
    for rxn, data in reactions_data.items():
        if data['direction'] == 'LEFT-TO-RIGHT':
            reactants = data['left']
            products = data['right']
            direction = 'uni'
        if data['direction'] == 'RIGHT-TO-LEFT':
            reactants = data['right']
            products = data['left']
            direction = 'uni'
        if data['direction'] == 'REVERSIBLE':
            reactants = data['left']
            products = data['right']
            direction = 'bi'
        lp_f.write(f'\n% {rxn}\n')
        rxn_to_asp(rxn, reactants, products, direction, lp_f)



def write_mappings(reactants, products, rxn, lp_f):
    lp_f.write(f'\n% {rxn}\n')
    for r in reactants:
        asso_if = dict()
        for atom in r.GetAtoms():
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            amap = atom.GetAtomMapNum()
            lp_f.write(f'atomMappingReactant("{rxn}","{symbol}",{amap}).\n')
            asso_if[idx] = amap
        for bond in r.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            bond_type = str(bond.GetBondType()).lower()
            lp_f.write(f'bondMappingReactant("{rxn}","{asso_if[a1]}","{asso_if[a2]}","{bond_type}").\n')

    for p in products:
        asso_if = dict()
        for atom in p.GetAtoms():
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            amap = atom.GetAtomMapNum()
            lp_f.write(f'atomMappingProduct("{rxn}","{symbol}",{amap}).\n')
            asso_if[idx] = amap
        for bond in p.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            bond_type = str(bond.GetBondType()).lower()
            lp_f.write(
                f'bondMappingProduct("{rxn}","{asso_if[a1]}","{asso_if[a2]}","{bond_type}").\n')


RUN_NAME = 'OXY'
RUN_PATH = '/home/phamongi/Documents/Dev/pathmodel/Files'
SOURCE = {'linoleate': 'CCCCC\C=C/C\C=C/CCCCCCCC([O-])=O'}
TARGET = {'12_hydroxy_13_glutation_OME': 'CCCCCC(SCC(NC(=O)CCC(N)C(=O)O)C(=O)NCC(=O)O)C(O)C/C=C/CCCCCCCC(=O)O'}
MC_REF_RXN = ['LEUKOTRIENE-C4-SYNTHASE-RXN', 'RXN-8495', 'LIPOXYGENASE-RXN']

generate_input(RUN_NAME, RUN_PATH, SOURCE, TARGET, MC_REF_RXN)