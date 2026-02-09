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

        % 2.1.1 To extract atoms added in the transformation we look for atoms present in Reactant but not in Product.
atomAddition(Transformation, AtomElement, AtomPosition) :- reaction(Transformation, Reactant, Product),
                                                           atom(Reactant, AtomElement, AtomPosition),
                                                           not atom(Product, AtomElement, AtomPosition).

        % 2.1.2 Same as above but for atoms present in Product and not in Reactant.
atomRemoval(Transformation, AtomElement, AtomPosition) :- reaction(Transformation, Reactant, Product),
                                                          atom(Product, AtomElement, AtomPosition),
                                                          not atom(Reactant, AtomElement, AtomPosition).


    % 2.2 Compare bonds between molecules.
    % ------------------------------------

        % 2.2.1 To extract bonds added in the transformation we look for bonds present in Reactant but not in Product.
bondAddition(Transformation, BondOrder, Atom1Position, Atom2Position) :- reaction(Transformation, Reactant, Product),
                                                                         bond(Reactant, BondOrder, Atom1Position, Atom2Position),
                                                                         not bond(Product, BondOrder, Atom1Position, Atom2Position).

        % 2.2.2 Same as above but for bonds present in Product and not in Reactant.
bondRemoval(Transformation, BondOrder, Atom1Position, Atom2Position) :- reaction(Transformation, Reactant, Product),
                                                                        bond(Product, BondOrder, Atom1Position, Atom2Position),
                                                                        not bond(Reactant, BondOrder, Atom1Position, Atom2Position).


%*
3. Identification du bloc lié à la transformation
========================================================================================================================
L'objectif est d'intentifier le bloc du CGR regroupant les atomes et liaisons directement liées au centre de réaction
*%

%   3.1 Identifier les atomes du reactant impliqués dans la transformation
%   ----------------------------------------------------------------------

% Lien créé entre ATOME et Atom2
reactantAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondAddition(Transformation, BondOrder, AtomPosition, Atom2Position),
                                                                      reaction(Transformation, Reactant, Product),
                                                                      atom(Reactant, AtomElement, AtomPosition).
% Lien créé entre ATOME et Atom1
reactantAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondAddition(Transformation, BondOrder, Atom1Position, AtomPosition),
                                                                      reaction(Transformation, Reactant, Product),
                                                                      atom(Reactant, AtomElement, AtomPosition).
% Lien supprimé entre ATOME et Atom2
reactantAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondRemoval(Transformation, BondOrder, AtomPosition, Atom2Position),
                                                                      reaction(Transformation, Reactant, Product),
                                                                      atom(Reactant, AtomElement, AtomPosition).
% Lien supprimé entre ATOME et Atom2
reactantAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondRemoval(Transformation, BondOrder, Atom1Position, AtomPosition),
                                                                      reaction(Transformation, Reactant, Product),
                                                                      atom(Reactant, AtomElement, AtomPosition).

%   3.2 Identifier les liens liés à ces atomes
%   ------------------------------------------


%   3.3 Identifier les atomes liés à ces liens
%   ------------------------------------------
