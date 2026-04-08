import csv
import json

MC_AMS_O = 'atom-mappings-smarts.json'
MC_R_O = 'reactions.json'
MC_C_O = 'compounds.json'


def generate_metacyc_data(atom_mappings_smiles_dat, reactions_dat, compounds_dat):
    # atom-mappings-smarts.json
    with open(atom_mappings_smiles_dat, 'r') as f, open(MC_AMS_O, 'w') as o:
        ams_d = {}
        ams_f = csv.reader(f, delimiter='\t')
        for l in ams_f:
            ams_d[l[0]] = l[1]
        json.dump(ams_d, o, indent=1)
    # reactions.json
    with open(reactions_dat, 'r') as f, open(MC_R_O, 'w') as o:
        r_d = {}
        for l in f:
            if l == '//\n':
                if r_id:
                    r_d[r_id] = {'left': left, 'right': right, 'direction': direction}
            if l.startswith('UNIQUE-ID'):
                r_id = l.strip().split(' - ')[1]
                left = []
                right = []
                direction = 'REVERSIBLE'
            if l.startswith('LEFT'):
                left.append(l.strip().split(' - ')[1])
            if l.startswith('RIGHT'):
                right.append(l.strip().split(' - ')[1])
            if l.startswith('REACTION-DIRECTION'):
                direction = l.strip().split(' - ')[1].replace('PHYSIOL-', '')
        json.dump(r_d, o, indent=1)
    # compounds.json
    with open(compounds_dat, 'r') as f, open(MC_C_O, 'w') as o:
        c_d = {}
        for l in f:
            if l.startswith('UNIQUE-ID'):
                c_id = l.strip().split(' - ')[1]
            if l.startswith('SMILES'):
                c_d[c_id] = l.strip().split(' - ')[1]
        json.dump(c_d, o, indent=1)


MC_AMSD = 'atom-mappings-smiles.dat'
MC_RD = 'reactions.dat'
MC_CD = 'compounds.dat'

generate_metacyc_data(MC_AMSD, MC_RD, MC_CD)