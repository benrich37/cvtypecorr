import pytest
from cvtypecorr.data.rxn import FormationReaction, reference_species, Reaction, RedoxReaction, AcidDeprotReaction, ThermoReaction
from cvtypecorr.data.pka import data_dict as pka_ref_dict, pka_to_delta_g, k_ev
from cvtypecorr.fit.rigid import Collection
import numpy as np


def test_collection_reg():
    rxn1 = Reaction("test1", ["H3O+", "NO3-"], ["H2O", "NO3H"], -1)
    test_collection = Collection(
        corrections={"H3O+": 0, "H2O": 0, "NO3H": 0}
    )
    assert all([specie in test_collection.corrections for specie in ["H3O+", "H2O", "NO3H"]])
    test_collection.apply_correction(rxn1, -2)
    assert test_collection.corrections["NO3-"] == -1

def test_collection_redox():
    nelec = 2
    delta_V = -1
    V_exp = -1
    rxn1 = RedoxReaction("test1", ["CO2", ("e-", nelec), "H+"], ["HCOO-"], V_exp)
    grounded_species = ["CO2", "H+"]
    test_collection = Collection(
        corrections={sp: 0 for sp in grounded_species}
    )
    assert all([specie in test_collection.corrections for specie in grounded_species])
    test_collection.apply_redox_correction(rxn1, V_exp + delta_V)
    assert test_collection.corrections["HCOO-"] == -nelec*delta_V

def test_collection_acid():
    acid = "NO3H"
    base = "NO3-"
    pka = pka_ref_dict[acid]
    delta_pka = -1
    rxn1 = AcidDeprotReaction("test1", acid, base, pka)
    grounded_species = [acid, "H3O+", "H2O"]
    test_collection = Collection(
        corrections={sp: 0 for sp in grounded_species}
    )
    assert all([specie in test_collection.corrections for specie in grounded_species])
    test_collection.apply_acid_correction(rxn1, pka + delta_pka)
    expected_correction = pka_to_delta_g(delta_pka, n_protons=1, T=rxn1.T)
    assert test_collection.corrections[base] == -expected_correction

# test_collection_reg()
# test_collection_redox()