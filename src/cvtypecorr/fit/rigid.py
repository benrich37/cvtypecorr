from cvtypecorr.data.rxn import Reaction, RedoxReaction, ThermoReaction, v0_SHE_CANDLE, eq_pot_to_electronless_reduction_energy, electronless_reduction_energy_to_eq_pot

# v0_SHE_CANDLE = 4.66

class Collection:
    def __init__(self, corrections: dict[str, float] | None = None, v0: float | None = None):
        if corrections is None:
            corrections = {}
        if v0 is None:
            v0 = v0_SHE_CANDLE
        self.corrections = corrections
        self.v0 = v0

    def apply_correction(self, rxn: Reaction, dft_nrg: float):
        e_err = dft_nrg - rxn.delta_E
        missing_reactants = [specie for specie in rxn.reactants if not specie[0] in self.corrections]
        missing_products = [specie for specie in rxn.products if not specie[0] in self.corrections]
        n_missing = len(missing_products) + len(missing_reactants)
        if n_missing > 1:
            raise ValueError(f"Rigid correction can only apply fit correction at a time (currently missing {missing_reactants} and {missing_products})")
        missing_specie = (missing_products + missing_reactants)[0][0]
        missing_coef = (missing_products + missing_reactants)[0][1]
        coef = -1 if len(missing_products) else 1
        self.corrections[missing_specie] = coef*(1/missing_coef)*e_err

    def apply_redox_correction(self, rxn: RedoxReaction, dft_eq_pot: float):
        n_elec = rxn.electrons
        target_electronless_reduction_energy = eq_pot_to_electronless_reduction_energy(rxn.equilibrium_potential, n_elec, self.v0)
        rxn.delta_E = target_electronless_reduction_energy
        dft_electronless_reduction_energy = eq_pot_to_electronless_reduction_energy(dft_eq_pot, n_elec, self.v0)
        self.corrections["e-"] = 0.0
        self.apply_correction(rxn, dft_electronless_reduction_energy)