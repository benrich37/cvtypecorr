# Step free energies on FeN4-i (eV), pH 6.3

Written by `surface_steps.py`. dG = dG_raw + correction. Configurations and calc dirs: `FeN4-i_steps.csv`.

| step | equation | U | dG_raw | correction | dG | correction terms |
|---|---|---|---:|---:|---:|---|
| * -> NO3 | * + NO3- -> *NO3 | +0.00 V_SHE | +0.193 | +0.000 | +0.193 | -c(NO3-) +c(NO3-) |
| * -> NO3 | * + NO3- -> *NO3 | -1.00 V_SHE | +0.393 | +0.000 | +0.393 | -c(NO3-) +c(NO3-) |
| * -> O | * + NO3- -> *O + NO2- | +0.00 V_SHE | +0.054 | -0.254 | -0.200 | -c(NO3-) +c(O) +c(NO2-) |
| * -> O | * + NO3- -> *O + NO2- | -1.00 V_SHE | -0.451 | -0.254 | -0.705 | -c(NO3-) +c(O) +c(NO2-) |
| * -> NO2 | * + NO2- -> *NO2 | +0.00 V_SHE | -0.536 | +0.000 | -0.536 | -c(NO2-) +c(NO2-) |
| * -> NO2 | * + NO2- -> *NO2 | -1.00 V_SHE | -0.446 | +0.000 | -0.446 | -c(NO2-) +c(NO2-) |
| * -> H | * + H3O+ -> *H + H2O | +0.00 V_SHE | +0.876 | -0.291 | +0.585 | -c(H3O+) +c(H) +c(H2O) |
| * -> H | * + H3O+ -> *H + H2O | -1.00 V_SHE | -0.265 | -0.291 | -0.557 | -c(H3O+) +c(H) +c(H2O) |
| NO3 -> NO3H | *NO3 + H3O+ -> *NO3H + H2O | +0.00 V_SHE | +0.597 | -0.074 | +0.523 | -c(NO3-) -c(H3O+) +c(NO3H) +c(H2O) |
| NO3 -> NO3H | *NO3 + H3O+ -> *NO3H + H2O | -1.00 V_SHE | +0.418 | -0.074 | +0.344 | -c(NO3-) -c(H3O+) +c(NO3H) +c(H2O) |
| NO3 -> O | *NO3 -> *O + NO2- | +0.00 V_SHE | -0.139 | -0.254 | -0.393 | -c(NO3-) +c(O) +c(NO2-) |
| NO3 -> O | *NO3 -> *O + NO2- | -1.00 V_SHE | -0.844 | -0.254 | -1.098 | -c(NO3-) +c(O) +c(NO2-) |
| NO3H -> NO2 | *NO3H + H3O+ -> *NO2 + 2 H2O | +0.00 V_SHE | -1.428 | -0.762 | -2.190 | -c(NO3H) -c(H3O+) +c(NO2-) +2c(H2O) |
| NO3H -> NO2 | *NO3H + H3O+ -> *NO2 + 2 H2O | -1.00 V_SHE | -3.359 | -0.762 | -4.122 | -c(NO3H) -c(H3O+) +c(NO2-) +2c(H2O) |
| O -> OH | *O + H3O+ -> *OH + H2O | +0.00 V_SHE | -0.240 | -0.819 | -1.059 | -c(O) -c(H3O+) +c(OH-) +c(H2O) |
| O -> OH | *O + H3O+ -> *OH + H2O | -1.00 V_SHE | -1.055 | -0.819 | -1.874 | -c(O) -c(H3O+) +c(OH-) +c(H2O) |
| OH -> * | *OH + H3O+ -> * + 2 H2O | +0.00 V_SHE | +0.085 | +0.236 | +0.321 | -c(OH-) -c(H3O+) +2c(H2O) |
| OH -> * | *OH + H3O+ -> * + 2 H2O | -1.00 V_SHE | -0.596 | +0.236 | -0.360 | -c(OH-) -c(H3O+) +2c(H2O) |
| H -> * | *H + H3O+ -> * + H2O + H2 | +0.00 V_SHE | +0.457 | -0.291 | +0.165 | -c(H) -c(H3O+) +c(H2O) +c(H2) |
| H -> * | *H + H3O+ -> * + H2O + H2 | -1.00 V_SHE | -0.402 | -0.291 | -0.694 | -c(H) -c(H3O+) +c(H2O) +c(H2) |
