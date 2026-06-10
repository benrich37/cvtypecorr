import pytest
from cvtypecorr.data.rxn import FormationReaction, reference_species, Reaction
from cvtypecorr.fit.rigid import Collection


def test_collection():
    rxn1 = Reaction("test1", ["H3O+", "NO3-"], ["H2O", "NO3H"], -1)
    test_collection = Collection(
        {"H3O+": 0, "H2O": 0, "NO3H": 0}
    )
    test_collection.apply_correction(rxn1, -2)
    assert test_collection.corrections["NO3-"] == -1

test_collection()