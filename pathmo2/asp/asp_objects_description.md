# Molecules representation

*A molecule "M1" is defined by a set of atoms and bonds*

### atom /3

`atom(MoleculeName, AtomElement, AtomPosition)`
- MoleculeName : (str) Molecule name containing the atom
- AtomElement : (cat) Element of the atom (ex : c for carbon, o for oxygen, etc...)
- AtomPosition : (int) Position of the atom in the 2D graph structure representation of the molecule

ex : `atom("M1", c, 1)`


### bond /4

`bond(MoleculeName, Atom1Position, Atom2Position, BondOrder)`
- MoleculeName : (str) Molecule name containing the bond
- Atom1Position : (int) Position of the 1st atom of the bond in the 2D graph structure 
representation of the molecule
- Atom2Position : (int) Position of the 2nd atom of the bond in the 2D graph structure 
representation of the molecule
- BondOrder : (cat) Order of the bond between the 2 atoms (ex : simple, double)

ex : `bond("M1", 1, 2, simple)`


# Reactions representation

*A reaction "R1" is defined by its reactant and its product*

### reaction /3
`reaction(ReactionName, Reactant, Product)`
- ReactionName : (str) The name of the reaction
- Reactant : (str) The molecule name (`MoleculeName`) of the reactant
- Product : (str) The molecule name (`MoleculeName`) of the product

ex : `reaction("R1", "M1", "M2")`

# Transformation representation

*A transformation "T1" is defined by the atoms and bonds differences (addition or removal)
in the reaction from which it derives*

### atomDifference /4

`atomDifference(ReactionName, AtomElement, AtomPosition, DifferenceType)`
- ReactionName : (str) Name of the reaction the transformation derives of
- AtomElement : (cat) Element of the atom (ex: c for carbon) of the atom different between
the product and the reactant
- AtomPosition : (int) Position of the atom different between the product and the reactant
- DifferenceType : (cat) Type of the difference, addition if present in product and not 
in reactant, removal if present in the reactant but not in the product

ex : `atomDifference("R1", c, 14, addition)`

### bondDifference /5

`bondDifference(ReactionName, Atom1Position, Atom2Position, BondOrder, DifferenceType)`
- ReactionName : (str) Name of the reaction the transformation derives of
- Atom1Position : (int) Position of the 1st atom of the bond
- Atom2Position : (int) Position of the 2nd atom of the bond
- BondOrder : (cat) Order of the bond between the 2 atoms 
- DifferenceType : (cat) Type of the difference, addition if present in product and not 
in reactant, removal if present in the reactant but not in the product

ex : `bondDifference("R1", 13, 14, , double, removal)`

# Reaction Center representation

*A Reaction Center is defined by a set of atoms and bonds directly linked with atoms implied in
the transformation from which it derives. We separate reaction center from the reactant and from
the product.*

## Reactant Reaction Center

### reactantAtomTransformed /3
*As intermediate, we extract atoms from the reactant implied in atoms or bonds modified in the 
transformation.*

`reactantAtomTransformed(ReactionName, AtomElement, AtomPosition)`
- ReactionName : (str) Name of the reaction the transformation derives of
- AtomElement : (cat) Element of the atom of the reactant implied
- AtomPosition : (str) Position of the atom of the reactant implied

### reactantTransformationCenterBond /4
*Extraction of all bonds implicating atoms transformed in the reactant (reactantAtomTransformed)*

`reactantTransformationCenterBond(ReactionName, Atom1Position, Atom2Position, BondOrder)`
- ReactionName : (str) Name of the reaction the transformation derives of
- Atom1Position : (int) Position of the 1st atom of the bond
- Atom2Position : (int) Position of the 2nd atom of the bond
- BondOrder : (cat) Order of the bond between the 2 atoms 

### reactantTransformationCenterAtom /3
*Extraction of all atoms implied in transformation center bonds of the reactant 
(reactantTransformationCenterBond)*

`reactantTransformationCenterAtom(ReactionName, AtomElement, AtomPosition)`
- ReactionName : (str) Name of the reaction the transformation derives of
- AtomElement : (cat) Element of the atom of the reactant implied
- AtomPosition : (str) Position of the atom of the reactant implied

## Product Reaction Center

### productAtomTransformed /3
*Same as reactantAtomTransformed for the product*

`productAtomTransformed(ReactionName, AtomElement, AtomPosition)`

### productTransformationCenterBond /4
*Same as reactantTransformationCenterBond for the product*

`productTransformationCenterBond(ReactionName, Atom1Position, Atom2Position, BondOrder)`

### productTransformationCenterAtom /3
*Same as reactantTransformationCenterAtom for the product*

`productTransformationCenterAtom(ReactionName, AtomElement, AtomPosition)`