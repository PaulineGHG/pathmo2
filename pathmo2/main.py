import os
import csv
import clyngor
import rdkit.Chem

from pathmo2.input_management.parse_input import *


ROOT = os.path.dirname(__file__)


def generate_transformations(run_path):
    """
    Detect reaction sites by comparing molecules implied in a reaction.
    Return the result as a string.

    """
    print('~~~~~Creation of Reaction~~~~~')
    transformations_script = os.path.join(*[ROOT, 'asp', 'TransformationsExtraction.lp'])
    input_file = os.path.join(run_path, INPUT_DIR, LP_INPUT)

    reaction_solver = clyngor.solve([input_file, transformations_script], use_clingo_module=False)
    result_atoms = []
    transformations = {}
    reactant_transformation_centers = {}
    product_transformation_centers = {}
    for atom in next(reaction_solver.parse_args.int_not_parsed.sorted):
        atom_str = atom[0] + '(' + ', '.join(atom[1]) + ')'
        result_atoms.append(atom_str)
        # Extract transformations
        if atom[0].endswith('Difference'):
            if atom[1][0] not in transformations:
                transformations[atom[1][0]] = []
            transformations[atom[1][0]].append(atom)
        # Extract Reactant Transformation Centers
        if atom[0].startswith('reactantTransformationCenter'):
            if atom[1][0] not in reactant_transformation_centers:
                reactant_transformation_centers[atom[1][0]] = []
            reactant_transformation_centers[atom[1][0]].append(atom)
        # Extract Product Transformation Centers
        if atom[0].startswith('productTransformationCenter'):
            if atom[1][0] not in product_transformation_centers:
                product_transformation_centers[atom[1][0]] = []
            product_transformation_centers[atom[1][0]].append(atom)

    print(result_atoms)
    print(transformations)
    print(product_transformation_centers)
    print(reactant_transformation_centers)
    return result_atoms, transformations, product_transformation_centers, reactant_transformation_centers

    # output_folder = os.path.join(run_path, OUTPUT_DIR)
    # pathmodel_output_transformation_path = os.path.join(output_folder, 'pathmodel_data_transformations.tsv')
    # with open(pathmodel_output_transformation_path, 'w') as transformation_file:
    #     csvwriter = csv.writer(transformation_file, delimiter = '\t')
    #     csvwriter.writerow(['reaction_id', 'reactant_sbustructure', 'product_substructure'])
    #     for reaction in reactions:
    #         if reaction in transformation_reactants:
    #             reactant = transformation_reactants[reaction]
    #         else:
    #             reactant = []
    #         if reaction in transformation_products:
    #             product = transformation_products[reaction]
    #         else:
    #             product = []
    #         csvwriter.writerow([reaction, reactant, product])
    #
    #
    # reaction_result = '\n'.join([atom+'.' for atom in reaction_results])
    #
    # return reaction_result


def export_transformations_patterns():
    product_transformation_pattern = rdkit.Chem.Mol


# ==================================================================================================

RUN_PATH = '/home/phamongi/Documents/Dev/pathmodel/Files'
# RUN_PATH = 'C:\\Users\Octav\PycharmProjects\pathmodel\Files'
RUN_NAME = 'ToyExemple'

# generate_lp_input(os.path.join(RUN_PATH, RUN_NAME), True)
generate_transformations(os.path.join(RUN_PATH, RUN_NAME))


