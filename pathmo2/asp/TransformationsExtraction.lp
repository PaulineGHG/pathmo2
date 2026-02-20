%*
This script will create the reference transformation extracting a reference reaction center.
CHANGE ? -> To do this it compares links and atoms between two molecules that are implied in a reaction.
When it finds difference, the script will extract the links and atoms and write the result in a file.
*%

%*
Chemical variable:
valence(Atom type, Number of valence)
*%
valence(c, 4). valence(n, 3). valence(o, 2). valence(p, 5).

bond(Molecule, Atom2Position, Atom1Position, BondOrder) :- bond(Molecule, Atom1Position, Atom2Position, BondOrder).

%*
1. ETAPE PRELIMINAIRE : ALIGNEMENT DES GRAPHES DE MOLECULES (AAM)
========================================================================================================================
*%


%*
2. REACTANT / PRODUCT DIFFERENCES EXTRACTION
========================================================================================================================

The idea is to model a Condensed Graph of Reaction (CGR) extracting only the differences between the reactant and the
product as the reaction center :
What was added : atomAddition / bondAddition
What was removed : atomRemoval / bondRemoval
*%


    % 2.1 Compare atoms between molecules.
    % ------------------------------------

        % 2.1.1 Atoms added in the transformation = present in Product but not in Reactant.
atomDifference(Transformation, AtomElement, AtomPosition, addition) :- reaction(Transformation, Reactant, Product),
                                                                       not atom(Reactant, AtomElement, AtomPosition),
                                                                       atom(Product, AtomElement, AtomPosition).

        % 2.1.2 Atoms removed in the transformation = present in Reactant but not in Product.
atomDifference(Transformation, AtomElement, AtomPosition, removal) :- reaction(Transformation, Reactant, Product),
                                                                      not atom(Product, AtomElement, AtomPosition),
                                                                      atom(Reactant, AtomElement, AtomPosition).


    % 2.2 Compare bonds between molecules.
    % ------------------------------------

        % 2.2.1 Bonds added in the transformation = present in Product but not in Reactant.
bondDifference(Transformation, Atom1Position, Atom2Position, BondOrder, addition) :- reaction(Transformation, Reactant, Product),
                                                                                     bond(Product, Atom1Position, Atom2Position, BondOrder),
                                                                                     not bond(Reactant, Atom1Position, Atom2Position, BondOrder).

        % 2.2.2 Bonds removed in the transformation = present in Reactant but not in Product.
bondDifference(Transformation, Atom1Position, Atom2Position, BondOrder, removal) :- reaction(Transformation, Reactant, Product),
                                                                                    bond(Reactant, Atom1Position, Atom2Position, BondOrder),
                                                                                    not bond(Product, Atom1Position, Atom2Position, BondOrder).


%*
3. Identification du bloc lié à la transformation
========================================================================================================================
L'objectif est d'intentifier le bloc du CGR regroupant les atomes et liaisons directement liées au centre de réaction
*%

%   3.1 Identifier les atomes du reactant impliqués dans la transformation
%   ----------------------------------------------------------------------

    % 3.1.1 REACTANT

        % Lien modifié entre ATOME et Atom2 dans le reactant
reactantAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondDifference(Transformation, AtomPosition, Atom2Position, BondOrder, Modification),
                                                                      reaction(Transformation, Reactant, Product),
                                                                      atom(Reactant, AtomElement, AtomPosition).


    % 3.1.2 PRODUCT

        % Lien modifié entre ATOME et Atom2 dans le produit
productAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondDifference(Transformation, AtomPosition, Atom2Position, BondOrder, Modification),
                                                                     reaction(Transformation, Reactant, Product),
                                                                     atom(Product, AtomElement, AtomPosition).



%   3.2 Identifier les liens liés à ces atomes
%   ------------------------------------------

    % 3.2.1 REACTANT

        % Lien entre un atome identifié et un autre atome dans le réactant
reactantTransformationCenterBond(Transformation, Atom1Position, Atom2Position, BondOrder) :- reactantAtomTransformed(Transformation, AtomElement, Atom1Position),
                                                                                             bond(Reactant, Atom1Position, Atom2Position, BondOrder),
                                                                                             reaction(Transformation, Reactant, Product).
        % Symétrie
reactantTransformationCenterBond(Transformation, Atom2Position, Atom1Position, BondOrder) :- reactantTransformationCenterBond(Transformation, Atom1Position, Atom2Position, BondOrder).

    % 3.2.2 PRODUCT

        % Lien entre un atome identifié et un autre atome dans le produit
productTransformationCenterBond(Transformation, Atom1Position, Atom2Position, BondOrder) :- productAtomTransformed(Transformation, AtomElement, Atom1Position),
                                                                                            bond(Product, Atom1Position, Atom2Position, BondOrder),
                                                                                            reaction(Transformation, Reactant, Product).
        % Symétrie
productTransformationCenterBond(Transformation, Atom2Position, Atom1Position, BondOrder) :- productTransformationCenterBond(Transformation, Atom1Position, Atom2Position, BondOrder).



%   3.3 Identifier les atomes liés à ces liens
%   ------------------------------------------

    % 3.3.1 REACTANT

reactantTransformationCenterAtom(Transformation, AtomElement, AtomPosition) :- reactantTransformationCenterBond(Transformation, AtomPosition, Atom2Position, BondOrder),
                                                                               atom(Reactant, AtomElement, AtomPosition),
                                                                               reaction(Transformation, Reactant, Product).

    % 3.3.2 PRODUCT

productTransformationCenterAtom(Transformation, AtomElement, AtomPosition) :- productTransformationCenterBond(Transformation, AtomPosition, Atom2Position, BondOrder),
                                                                              atom(Product, AtomElement, AtomPosition),
                                                                              reaction(Transformation, Reactant, Product).



%#show atomAddition/3.
%#show atomRemoval/3.
%#show bondAddition/4.
%#show bondRemoval/4.
#show atomDifference/4.
#show bondDifference/5.
#show reactantAtomTransformed/3.
#show productAtomTransformed/3.
#show reactantTransformationCenterBond/4.
#show productTransformationCenterBond/4.
#show reactantTransformationCenterAtom/3.
#show productTransformationCenterAtom/3.