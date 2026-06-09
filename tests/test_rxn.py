import pytest
from cvtypecorr.data.rxn import FormationReaction, reference_species


@pytest.mark.parametrize("name, formula, delta_g, el_by_priority", [
    ("C2H3OH", [("C", 2), ("O", 1), ("H", 5)], -10.0, None),
    ("C2H3OH", [("C", 2), ("O", 1), ("H", 5)], -10.0, ["C", "O", "H"]),
])
def test_formation_reaction(name, formula, delta_g, el_by_priority):
    test = FormationReaction(name, formula, delta_g, el_by_priority)
    assert test.name == name
    assert test.delta_g == delta_g
    assert test.el_by_priority == (el_by_priority if el_by_priority is not None else ["O", "C", "H"])
    # Check that the reactants and products have the correct formula
    reactant_formula = {}
    for specie_name, coeff in test.reactants:
        specie = reference_species[specie_name]
        for el, el_coeff in specie.formula:
            reactant_formula[el] = reactant_formula.get(el, 0) + el_coeff * coeff
    for el, el_coeff in formula:
        assert reactant_formula.get(el, 0) == el_coeff, f"Element {el} has incorrect coefficient in reactants (expected {el_coeff}, got {reactant_formula.get(el, 0)})"
    product_formula = {}
    for specie_name, coeff in test.products:
        specie = reference_species[specie_name]
        for el, el_coeff in specie.formula:
            product_formula[el] = product_formula.get(el, 0) + el_coeff * coeff
    for el, el_coeff in formula:
        assert product_formula.get(el, 0) == el_coeff, f"Element {el} has incorrect coefficient in products (expected {el_coeff}, got {product_formula.get(el, 0)})"
    el_by_priority_missing_one = test.el_by_priority[:-1]
    with pytest.raises(AssertionError):
        FormationReaction(name, formula, delta_g, el_by_priority_missing_one)