from pymatgen.core.structure import Element


class Reaction:
    def __init__(self, name, reactants, products):
        self.name = name
        self.reactants = reactants
        self.products = products

def get_specie_charge(specie_name):
    if "-" in specie_name:
        name, charge = specie_name.split("-")
        return -int(charge)
    elif "n" in specie_name:
        name, charge = specie_name.split("n")
        return -int(charge)
    elif "p" in specie_name:
        name, charge = specie_name.split("p")
        return int(charge)
    elif "+" in specie_name:
        name, charge = specie_name.split("+")
        return int(charge)
    else:
        return 0

def get_charge_balance(reactants, products):
    reactant_charge = sum(get_specie_charge(specie_name) * coeff for specie_name, coeff in reactants)
    product_charge = sum(get_specie_charge(specie_name) * coeff for specie_name, coeff in products)
    return reactant_charge - product_charge

def get_num_electrons(reactants, products):
    charge_balance = get_charge_balance(reactants, products)
    return -charge_balance

class RedoxReaction(Reaction):
    def __init__(self, name, reactants, products, equilibrium_potential):
        super().__init__(name, reactants, products)
        self.electrons = get_num_electrons(reactants, products)
        self.equilibrium_potential = equilibrium_potential

class ThermoReaction(Reaction):
    def __init__(self, name, reactants, products, delta_g):
        super().__init__(name, reactants, products)
        self.delta_g = delta_g

class ReferenceSpecie:
    def __init__(self, rep_el, name, formula: list[tuple[str, int | float]]):
        assert rep_el in [p[0] for p in formula], "Representative element must be in the formula"
        self.rep_el = rep_el
        self.name = name
        self.formula = formula
        self.depends_on = {p[0]: p[1] for p in formula if p[0] != rep_el}

reference_species = {
    "O": ReferenceSpecie("O", "H2O", [("H", 2), ("O", 1)]),
    "H": ReferenceSpecie("H", "H2", [("H", 2)]),
    "N": ReferenceSpecie("N", "NH3", [("N", 1), ("H", 3)]),
    "C": ReferenceSpecie("C", "CH4", [("C", 1), ("H", 4)]),
}

def get_cur_reactants_formula_dict(reactants: list[tuple[str, int | float]]):
    formula = {}
    for specie_name, coeff in reactants:
        ref_specie = [p for p in reference_species.values() if p.name == specie_name]
        if len(ref_specie) == 0:
            raise ValueError(f"Specie {specie_name} not found in reference species")
        ref_specie = ref_specie[0]
        for el, el_coeff in ref_specie.formula:
            formula[el] = formula.get(el, 0) + el_coeff * coeff
    return formula

def formula_dict_is_formula(formula_dict: dict[str, int | float], formula: list[tuple[str, int | float]]):
    for el, coeff in formula:
        if formula_dict.get(el, 0) != coeff:
            return False
    return True


def get_all_els(formula: list[tuple[str, int | float]]):
    els = set(p[0] for p in formula)
    all_els = set([el for el in els])
    for el in els:
        assert el in reference_species, f"Element {el} does not have a reference specie defined"
        depends_on = reference_species[el].depends_on
        all_els |= set(list(depends_on.keys()))
    for el in all_els - els:
        assert el in reference_species, f"Reference-dependent Eeement {el} does not have a reference specie defined"
    return all_els

class FormationReaction(ThermoReaction):
    def __init__(self, name, formula, delta_g, el_by_priority: list | None = None):
        all_els = get_all_els(formula)
        if el_by_priority is None:
            els = list(all_els)
            Zs = [Element(el).Z for el in els]
            el_by_priority = [x for _, x in sorted(zip(Zs, els), reverse=True)]
        else:
            assert set(el_by_priority) == all_els, f"Provided el_by_priority does not match the elements in the formula (missing {all_els - set(el_by_priority)}, extra {set(el_by_priority) - all_els})"
        self.el_by_priority = el_by_priority
        products = [(name, 1)]
        reactants = []
        cur_formula_dict = get_cur_reactants_formula_dict(reactants)
        finished = formula_dict_is_formula(cur_formula_dict, formula)
        while not finished:
            for el in el_by_priority:
                assert el in reference_species, f"Element {el} does not have a reference specie defined"
                target_coeff = next(p[1] for p in formula if p[0] == el)
                cur_formula_dict = get_cur_reactants_formula_dict(reactants)
                delta_coeff = target_coeff - cur_formula_dict.get(el, 0)
                ref_specie = reference_species[el]
                ref_formula = ref_specie.formula
                ref_coeff = delta_coeff / next(p[1] for p in ref_formula if p[0] == el)
                reactants.append((ref_specie.name, ref_coeff))
                finished = formula_dict_is_formula(get_cur_reactants_formula_dict(reactants), formula)
        super().__init__(name, reactants, products, delta_g)


    


test = FormationReaction("C2H3OH", [("C", 2), ("O", 1), ("H", 5)], -10.0)
print("name", test.name)
print("reactants", test.reactants)
print("products", test.products)
print(test.delta_g)
print(test.el_by_priority)

# test = FormationReaction("C2H3OH", [("C", 2), ("O", 1), ("H", 5)], -10.0, el_by_priority=["O"])

# redox_reactions = {
#     # "no3h_to_no_redox": RedoxReaction(
#     #     "no3h_to_no_redox",
#     #     ["NO3H"],
#     #     ["NO"],
#     # )
# }

# thermo_reactions = {

# }