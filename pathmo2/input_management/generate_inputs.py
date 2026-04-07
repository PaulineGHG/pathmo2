import csv
import os.path

from rdkit.Chem import AllChem
from rdkit.Chem import Draw

METACYC_MAPPING_FILE = '../data/atom-mappings-smiles.dat'


def generate_input(run_name, output_path, source, target, metacyc_ref_reactions):
    run_path = os.path.join(output_path, run_name)
    if not os.path.exists(run_path):
        os.mkdir(run_path)
        os.mkdir(os.path.join(run_path, 'Inputs'))
        os.mkdir(os.path.join(run_path, 'Outputs'))

    reactions_file = os.path.join(run_path, 'Inputs', 'Reactions_references.tsv')
    create_reactions_file(reactions_file, metacyc_ref_reactions, run_path)


def create_reactions_file(reactions_file, metacyc_ref_reactions, run_path):
    if metacyc_ref_reactions:
        mc_ref_mappings = {}
        with open(METACYC_MAPPING_FILE, 'r') as f:
            mc_map = csv.reader(f, delimiter='\t')
            for l in mc_map:
                if l[0] in metacyc_ref_reactions:
                    mc_ref_mappings[l[0]] = l[1]

        for rxn_id, mapping in mc_ref_mappings.items():
            rxn = AllChem.ReactionFromSmarts(mapping)
            reactants = rxn.GetReactants()
            products = rxn.GetProducts()
            i = 0
            for r in reactants:
                i += 1
                r_id = f'M{i}'
                for atom in r.GetAtoms():
                    idx = atom.GetIdx()
                    symbol = atom.GetSymbol()
                    amap = atom.GetAtomMapNum()
                    print(r_id, idx, symbol, amap)
                for bond in r.GetBonds():
                    a1 = bond.GetBeginAtomIdx()
                    a2 = bond.GetEndAtomIdx()
                    bond_type = bond.GetBondType()
                    print(a1, a2, bond_type)
            print(products)

            draw_rxn(rxn, os.path.join(run_path, 'Outputs', rxn_id + '.png'))


def draw_rxn(rxn, output):
    d2d = Draw.MolDraw2DCairo(1600, 600)
    d2d.DrawReaction(rxn, highlightByReactant=False)
    png = d2d.GetDrawingText()
    open(output, 'wb+').write(png)



RUN_NAME = 'OXY'
RUN_PATH = '/home/phamongi/Documents/Dev/pathmodel/Files'
SOURCE = ('linoleate', 'CCCCC\C=C/C\C=C/CCCCCCCC([O-])=O')
TARGET = ('12_hydroxy_13_glutation_OME', 'CCCCCC(SCC(NC(=O)CCC(N)C(=O)O)C(=O)NCC(=O)O)C(O)C/C=C/CCCCCCCC(=O)O')
MC_REF_RXN = ['LEUKOTRIENE-C4-SYNTHASE-RXN', 'RXN-8495', 'LIPOXYGENASE-RXN']

generate_input(RUN_NAME, RUN_PATH, SOURCE, TARGET, MC_REF_RXN)