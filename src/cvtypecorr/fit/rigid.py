from cvtypecorr.data.rxn import Reaction

class Collection:
    def __init__(self, corrections: dict[str, float] | None = None):
        if corrections is None:
            corrections = {}
        self.corrections = products

    def apply_correction(self, rxn: Reaction, dft_nrg: float):
        e_err = dft_nrg - rxn.delta_E
        missing_reactants = [specie[0] for specie in rxn.reactants if not specie in self.corrections]
        missing_products = [specie[0] for specie in rxn.products if not specie in self.corrections]
        n_missing = len(missing_products) + len(missing_reactants)
        if n_missing > 1:
            raise ValueError("Rigid correction can only apply fit correction at a time")
        missing_specie = (missing_products + missing_reactants)[0][0]
        missing_coef = (missing_products + missing_reactants)[0][1]
        coef = 1 if len(missing_products) else -1
        self.corrections[missing_specie] = coef*(1/missing_coef)*e_err


    