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
atomAddition(Transformation, AtomElement, AtomPosition) :- reaction(Transformation, Reactant, Product),
                                                           not atom(Reactant, AtomElement, AtomPosition),
                                                           atom(Product, AtomElement, AtomPosition).

        % 2.1.2 Atoms removed in the transformation = present in Reactant but not in Product.
atomRemoval(Transformation, AtomElement, AtomPosition) :- reaction(Transformation, Reactant, Product),
                                                          not atom(Product, AtomElement, AtomPosition),
                                                          atom(Reactant, AtomElement, AtomPosition).


    % 2.2 Compare bonds between molecules.
    % ------------------------------------

        % 2.2.1 Bonds added in the transformation = present in Product but not in Reactant.
bondAddition(Transformation, Atom1Position, Atom2Position, BondOrder) :- reaction(Transformation, Reactant, Product),
                                                                         bond(Product, Atom1Position, Atom2Position, BondOrder),
                                                                         not bond(Reactant, Atom1Position, Atom2Position, BondOrder).

        % 2.2.2 Bonds removed in the transformation = present in Reactant but not in Product.
bondRemoval(Transformation, Atom1Position, Atom2Position, BondOrder) :- reaction(Transformation, Reactant, Product),
                                                                        bond(Reactant, Atom1Position, Atom2Position, BondOrder),
                                                                        not bond(Product, Atom1Position, Atom2Position, BondOrder).


%*
3. Identification du bloc lié à la transformation
========================================================================================================================
L'objectif est d'intentifier le bloc du CGR regroupant les atomes et liaisons directement liées au centre de réaction
*%

%   3.1 Identifier les atomes du reactant impliqués dans la transformation
%   ----------------------------------------------------------------------

    % REACTANT

% Lien créé entre ATOME et Atom2
reactantAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondAddition(Transformation, AtomPosition, Atom2Position, BondOrder),
                                                                      reaction(Transformation, Reactant, Product),
                                                                      atom(Reactant, AtomElement, AtomPosition).
% Lien créé entre ATOME et Atom1
reactantAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondAddition(Transformation, Atom1Position, AtomPosition, BondOrder),
                                                                      reaction(Transformation, Reactant, Product),
                                                                      atom(Reactant, AtomElement, AtomPosition).
% Lien supprimé entre ATOME et Atom2
reactantAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondRemoval(Transformation, AtomPosition, Atom2Position, BondOrder),
                                                                      reaction(Transformation, Reactant, Product),
                                                                      atom(Reactant, AtomElement, AtomPosition).
% Lien supprimé entre ATOME et Atom2
reactantAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondRemoval(Transformation, Atom1Position, AtomPosition, BondOrder),
                                                                      reaction(Transformation, Reactant, Product),
                                                                      atom(Reactant, AtomElement, AtomPosition).

    % PRODUCT

% Lien créé entre ATOME et Atom2
productAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondAddition(Transformation, AtomPosition, Atom2Position, BondOrder),
                                                                     reaction(Transformation, Reactant, Product),
                                                                     atom(Product, AtomElement, AtomPosition).
% Lien créé entre ATOME et Atom1
productAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondAddition(Transformation, Atom1Position, AtomPosition, BondOrder),
                                                                     reaction(Transformation, Reactant, Product),
                                                                     atom(Product, AtomElement, AtomPosition).
% Lien supprimé entre ATOME et Atom2
productAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondRemoval(Transformation, AtomPosition, Atom2Position, BondOrder),
                                                                     reaction(Transformation, Reactant, Product),
                                                                     atom(Product, AtomElement, AtomPosition).
% Lien supprimé entre ATOME et Atom2
productAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondRemoval(Transformation, Atom1Position, AtomPosition, BondOrder),
                                                                     reaction(Transformation, Reactant, Product),
                                                                     atom(Product, AtomElement, AtomPosition).

%   3.2 Identifier les liens liés à ces atomes
%   ------------------------------------------

% Lien entre un atome identifié et un autre atome dans le réactant
reactantTransformationCenterBond(Transformation, Atom1Position, Atom2Position, BondOrder) :- reactantAtomTransformed(Transformation, AtomElement, Atom1Position),
                                                                                             bond(Reactant, Atom1Position, Atom2Position, BondOrder),
                                                                                             reaction(Transformation, Reactant, Product).

% Lien entre un autre atome et un atome identifié dans le réactant
reactantTransformationCenterBond(Transformation, Atom1Position, Atom2Position, BondOrder) :- reactantAtomTransformed(Transformation, AtomElement, Atom2Position),
                                                                                             bond(Reactant, Atom1Position, Atom2Position, BondOrder),
                                                                                             reaction(Transformation, Reactant, Product).

% Lien entre un atome identifié et un autre atome dans le produit
productTransformationCenterBond(Transformation, Atom1Position, Atom2Position, BondOrder) :- productAtomTransformed(Transformation, AtomElement, Atom1Position),
                                                                                            bond(Product, Atom1Position, Atom2Position, BondOrder),
                                                                                            reaction(Transformation, Reactant, Product).

% Lien entre un autre atome et un atome identifié dans le produit
productTransformationCenterBond(Transformation, Atom1Position, Atom2Position, BondOrder) :- productAtomTransformed(Transformation, AtomElement, Atom2Position),
                                                                                            bond(Product, Atom1Position, Atom2Position, BondOrder),
                                                                                            reaction(Transformation, Reactant, Product).


%   3.3 Identifier les atomes liés à ces liens
%   ------------------------------------------

reactantTransformationCenterAtom(Transformation, AtomElement, AtomPosition) :- reactantTransformationCenterBond(Transformation, AtomPosition, Atom2Position, BondOrder),
                                                                               atom(Reactant, AtomElement, AtomPosition),
                                                                               reaction(Transformation, Reactant, Product).

reactantTransformationCenterAtom(Transformation, AtomElement, AtomPosition) :- reactantTransformationCenterBond(Transformation, Atom2Position, AtomPosition, BondOrder),
                                                                               atom(Reactant, AtomElement, AtomPosition),
                                                                               reaction(Transformation, Reactant, Product).

productTransformationCenterAtom(Transformation, AtomElement, AtomPosition) :- productTransformationCenterBond(Transformation, AtomPosition, Atom2Position, BondOrder),
                                                                              atom(Product, AtomElement, AtomPosition),
                                                                              reaction(Transformation, Reactant, Product).

productTransformationCenterAtom(Transformation, AtomElement, AtomPosition) :- productTransformationCenterBond(Transformation, Atom2Position, AtomPosition, BondOrder),
                                                                              atom(Product, AtomElement, AtomPosition),
                                                                              reaction(Transformation, Reactant, Product).

#show atomAddition/3.
#show atomRemoval/3.
#show bondAddition/4.
#show bondRemoval/4.
#show reactantAtomTransformed/3.
#show productAtomTransformed/3.
#show reactantTransformationCenterBond/4.
#show productTransformationCenterBond/4.
#show reactantTransformationCenterAtom/3.
#show productTransformationCenterAtom/3.