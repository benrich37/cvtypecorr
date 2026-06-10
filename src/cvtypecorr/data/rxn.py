from pymatgen.core.structure import Element, Composition

v0_SHE_CANDLE = 4.66
def resolve_specie_set(specie_set: list[str | tuple[str, int]]):
    resolved_specie_set = []
    for specie in specie_set:
        if isinstance(specie, str):
            resolved_specie_set.append((specie, 1))
        elif isinstance(specie, list) or isinstance(specie, tuple):
            if len(specie) == 1:
                resolved_specie_set.append((specie[0], 1))
            else:
                resolved_specie_set.append(tuple(specie))
    return resolved_specie_set

class Reaction:
    def __init__(self, name, reactants, products, delta_E):
        self.name = name
        reactants = resolve_specie_set(reactants)
        products = resolve_specie_set(products)
        charge_balance = get_charge_balance(reactants, products)
        if charge_balance != 0:
            raise ValueError(f"Reaction is not charge balanced. Charge balance: {charge_balance}")
        atom_balances = get_atom_balances(reactants, products)
        if not all(balance == 0 for balance in atom_balances.values()):
            imbalanced_els = {el: balance for el, balance in atom_balances.items() if balance != 0}
            raise ValueError(f"Reaction is not atom balanced. Imbalanced elements: {imbalanced_els}")
        self.reactants = reactants
        self.products = products
        self.delta_E = delta_E

def get_name_and_charge(specie_name):
    if "-" in specie_name:
        name, charge = specie_name.split("-")
        charge = -1 if not len(charge) else -int(charge)
    elif "+" in specie_name:
        name, charge = specie_name.split("+")
        charge = 1 if not len(charge) else int(charge)
    else:
        name = specie_name
        charge = 0
    return name, charge

def get_specie_charge(specie_name):
    return get_name_and_charge(specie_name)[1]

def get_specie_name(specie_name):
    return get_name_and_charge(specie_name)[0]

def get_charge_balance(reactants, products):
    reactant_charge = sum(get_specie_charge(specie_name) * coeff for specie_name, coeff in reactants)
    product_charge = sum(get_specie_charge(specie_name) * coeff for specie_name, coeff in products)
    return reactant_charge - product_charge

def remove_atomless_species(resolved_specie_set):
    exclude_species = ["e-"]
    cleaned_specie_set = []
    for specie_name, coeff in resolved_specie_set:
        if specie_name not in exclude_species:
            cleaned_specie_set.append((specie_name, coeff))
    return cleaned_specie_set

def get_atom_balances(reactants, products):
    creactants = remove_atomless_species(reactants)
    cproducts = remove_atomless_species(products)
    reactant_compositions = [Composition(get_specie_name(specie_name)).as_dict() for specie_name, coeff in creactants]
    product_compositions = [Composition(get_specie_name(specie_name)).as_dict() for specie_name, coeff in cproducts]
    # all_els = set(list(rc.keys()) for rc in reactant_compositions) | set(list(pc.keys()) for pc in product_compositions)
    all_els = set()
    for compset in [reactant_compositions, product_compositions]:
        for rc in compset:
            all_els |= set(list(rc.keys()))
    # atom_balances = {el: sum(rc.as_dict().get(el, 0) for rc in reactant_compositions) - sum(pc.as_dict().get(el, 0) for pc in product_compositions) for el in all_els}
    atom_balances = {el: 0 for el in all_els}
    for i, pdt_comp_dict in enumerate(product_compositions):
        for el, coeff in pdt_comp_dict.items():
            atom_balances[el] += coeff * cproducts[i][1]
    for i, rct_comp_dict in enumerate(reactant_compositions):
        for el, coeff in rct_comp_dict.items():
            atom_balances[el] -= coeff * creactants[i][1]
    return atom_balances

def get_num_electrons(reactants, products):
    n_elec_reactants = sum(get_specie_charge(specie_name) * coeff for specie_name, coeff in reactants if specie_name == "e-")
    n_elec_products = sum(get_specie_charge(specie_name) * coeff for specie_name, coeff in products if specie_name == "e-")
    return n_elec_reactants - n_elec_products

class RedoxReaction(Reaction):
    def __init__(self, name, reactants, products, equilibrium_potential, v0: float | None = None):
        if v0 is None:
            v0 = v0_SHE_CANDLE
        self._v0 = v0
        super().__init__(name, reactants, products, 0.0)
        self.electrons = get_num_electrons(self.reactants, self.products)
        self.equilibrium_potential = equilibrium_potential
        self.delta_E = eq_pot_to_electronless_reduction_energy(self.equilibrium_potential, self.electrons, self.v0)

    # @property
    # def delta_E(self):
    #     return eq_pot_to_electronless_reduction_energy(self.equilibrium_potential, self.electrons, self.v0)

    @property
    def v0(self):
        return self._v0

    @v0.setter
    def v0(self, new_v0):
        self._v0 = new_v0
        self.delta_E = eq_pot_to_electronless_reduction_energy(self.equilibrium_potential, self.electrons, self._v0)

class ThermoReaction(Reaction):
    def __init__(self, name, reactants, products, delta_g):
        super().__init__(name, reactants, products, delta_g)
        self.delta_g = delta_g

class ReferenceSpecie:
    def __init__(self, rep_el, name, formula: list[tuple[str, int | float]] | None = None):
        if formula is None:
            formula_dict = Composition(name).as_dict()
            formula = [(el, coeff) for el, coeff in formula_dict.items()]
        assert rep_el in [p[0] for p in formula], "Representative element must be in the formula"
        self.rep_el = rep_el
        self.name = name
        self.formula = formula
        self.depends_on = {p[0]: p[1] for p in formula if p[0] != rep_el}

reference_species = {
    "O": ReferenceSpecie("O", "H2O"),
    "H": ReferenceSpecie("H", "H2"),
    "N": ReferenceSpecie("N", "NH3"),
    "C": ReferenceSpecie("C", "CH4"),
}

def get_cur_reactants_formula_dict(reactants: list[tuple[str, int | float]]):
    # reactants must be a list of species in reference_species
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


def mu_to_voltage(mu, v0=v0_SHE_CANDLE):
    voltage = (- mu) - v0
    return voltage

def voltage_to_mu(voltage, v0=v0_SHE_CANDLE):
    mu = - (voltage + v0)
    return mu

def eq_pot_to_electronless_reduction_energy(eq_pot, n_electrons, v0=v0_SHE_CANDLE):
    mu = voltage_to_mu(eq_pot, v0)
    electronless_reduction_energy = mu*n_electrons
    return electronless_reduction_energy

def electronless_reduction_energy_to_eq_pot(electronless_reduction_energy, n_electrons, v0=v0_SHE_CANDLE):
    mu = electronless_reduction_energy/n_electrons
    eq_pot = mu_to_voltage(mu, v0)
    return eq_pot