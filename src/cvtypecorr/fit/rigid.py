from cvtypecorr.data.rxn import (
    Reaction, RedoxReaction, ThermoReaction, 
    v0_SHE_CANDLE, eq_pot_to_electronless_reduction_energy, 
    electronless_reduction_energy_to_eq_pot, 
    AcidDeprotReaction,
    pka_to_delta_g,
)

# v0_SHE_CANDLE = 4.66

class Collection:
    def __init__(self, corrections: dict[str, float] | None = None, v0: float | None = None):
        if corrections is None:
            corrections = {}
        if v0 is None:
            v0 = v0_SHE_CANDLE
        self.corrections = corrections
        self.v0 = v0

    def apply_correction(self, rxn: Reaction, dft_nrg: float, overwrite: str | None = None):
        e_err = dft_nrg - rxn.delta_E
        if not overwrite is None:
            self.corrections.pop(overwrite, None)
        missing_reactants = [specie for specie in rxn.reactants if not specie[0] in self.corrections]
        missing_products = [specie for specie in rxn.products if not specie[0] in self.corrections]
        n_missing = len(missing_products) + len(missing_reactants)
        if n_missing > 1:
            raise ValueError(f"Rigid correction can only apply fit correction at a time (currently missing {missing_reactants} and {missing_products})")
        missing_specie = (missing_products + missing_reactants)[0][0]
        missing_coef = (missing_products + missing_reactants)[0][1]
        coef = -1 if len(missing_products) else 1
        self.corrections[missing_specie] = coef*(1/missing_coef)*e_err

    def apply_redox_correction(self, rxn: RedoxReaction, dft_eq_pot: float, overwrite: str | None = None):
        n_elec = rxn.electrons
        rxn.v0 = self.v0
        dft_electronless_reduction_energy = eq_pot_to_electronless_reduction_energy(dft_eq_pot, n_elec, rxn.v0)
        self.corrections["e-"] = 0.0
        self.apply_correction(rxn, dft_electronless_reduction_energy, overwrite=overwrite)

    def apply_acid_correction(self, rxn: AcidDeprotReaction, dft_pka: float, overwrite: str | None = None):
        dft_delta_g = pka_to_delta_g(dft_pka, n_protons=1, T=rxn.T)
        self.apply_correction(rxn, dft_delta_g, overwrite=overwrite)
        # n_elec = rxn.electrons
        # rxn.v0 = self.v0
        # dft_electronless_reduction_energy = eq_pot_to_electronless_reduction_energy(dft_eq_pot, n_elec, rxn.v0)
        # self.corrections["e-"] = 0.0
        # self.apply_correction(rxn, dft_electronless_reduction_energy)