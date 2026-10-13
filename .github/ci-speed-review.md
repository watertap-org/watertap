# CI test runtime review

The saved timing reports are single runs, so these numbers are leads rather than speedup evidence. The table covers every file whose setup, call, and teardown time reached five seconds in one of the five proposed PR environments. The environment shown is the slowest saved PR cell for that file. Timings exclude collection, process startup, and CI setup.

| Test file | Slowest PR cell | Phase seconds (setup / call) | Review |
| --- | --- | ---: | --- |
| `watertap/unit_models/tests/test_membrane_distillation_1D.py` | linux-py3.10 | 16.5 / 445.9 | Shared harness rebuilds an unused model in two checks; harness change applies. |
| `watertap/flowsheets/dye_desalination/tests/test_dye_sweep.py` | linux-py3.10 | 0.0 / 217.5 | Each parameter selects a distinct sweep; preserve all 19 cases. |
| `watertap/unit_models/tests/test_electrodialysis_bipolar_1D.py` | linux-py3.10 | 4.4 / 210.5 | Distinct configurations and six pressure-drop methods; class fixtures already share setup. |
| `watertap/unit_models/tests/test_osmotically_assisted_reverse_osmosis_1D.py` | linux-py3.10 | 2.7 / 99.0 | Shared harness rebuilds an unused model in two checks; harness change applies. |
| `watertap/flowsheets/full_water_resource_recovery_facility/test/test_BSM2_P_extension.py` | linux-py3.10 | 63.4 / 25.3 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/flowsheets/MD/tests/test_MD.py` | linux-py3.10 | 0.4 / 85.7 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/flowsheets/full_water_resource_recovery_facility/test/test_full_WRRF_with_ASM1_ADM1.py` | linux-py3.10 | 30.5 / 51.4 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/unit_models/tests/test_electrodialysis_1D.py` | linux-py3.10 | 2.7 / 68.8 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/unit_models/tests/test_reverse_osmosis_1D.py` | linux-py3.10 | 1.1 / 59.8 | Shared harness rebuilds an unused model in two checks; harness change applies. |
| `watertap/flowsheets/lsrro/tests/test_lssro_multi_sweep.py` | win64-py3.12 | 0.0 / 49.1 | Distinct end-to-end run(s); no reusable setup shown. |
| `watertap/unit_models/tests/test_gac.py` | linux-py3.10 | 5.9 / 42.0 | Shared harness rebuilds an unused model in two checks; harness change applies. |
| `watertap/unit_models/tests/test_nanofiltration_DSPMDE_0D.py` | linux-py3.10 | 1.1 / 45.9 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/flowsheets/dye_desalination/tests/test_dye_desalination.py` | linux-py3.10 | 2.3 / 41.2 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/flowsheets/flex_desal/tests/test_pricetaker_workflow.py` | linux-py3.10 | 18.7 / 23.0 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/flowsheets/lsrro/tests/test_lssro_paper_analysis.py` | win64-py3.12 | 0.0 / 40.6 | Each CSV row is a distinct baseline; collect all 15 on PR. |
| `watertap/unit_models/tests/test_membrane_distillation_0D.py` | linux-py3.10 | 2.0 / 38.0 | Shared harness rebuilds an unused model in two checks; harness change applies. |
| `watertap/flowsheets/electrodialysis/tests/test_electrodialysis_1stack.py` | linux-py3.10 | 1.0 / 33.7 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/flowsheets/mvc/tests/test_mvc_single_stage.py` | win64-py3.12 | 0.2 / 30.9 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/unit_models/tests/test_reverse_osmosis_0D.py` | linux-py3.10 | 0.5 / 30.4 | Shared harness rebuilds an unused model in two checks; harness change applies. |
| `watertap/flowsheets/oaro/tests/test_oaro_multi.py` | linux-py3.10 | 0.6 / 28.5 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/flowsheets/oaro/tests/test_oaro.py` | linux-py3.10 | 0.3 / 27.7 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/flowsheets/electrodialysis/tests/test_electrodialysis_1stack_conc_recirc.py` | linux-py3.10 | 0.7 / 27.1 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/flowsheets/gac/tests/test_gac.py` | linux-py3.10 | 0.1 / 26.2 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/flowsheets/RO_with_energy_recovery/tests/test_monte_carlo_sampling_RO_ERD.py` | linux-py3.10 | 0.0 / 15.5 | Distinct end-to-end run(s); no reusable setup shown. |
| `watertap/flowsheets/RO_with_energy_recovery/tests/test_RO_with_energy_recovery_simulation.py` | linux-py3.10 | 0.3 / 13.3 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/unit_models/tests/test_thickener_unit.py` | win64-py3.12 | 0.2 / 13.2 | Cost is in case calls; no safe reusable setup shown. |
| `watertap/unit_models/tests/test_dewatering_unit.py` | linux-py3.10 | 0.2 / 12.8 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/flowsheets/electroNP/tests/test_BSM2_electroNP_no_bioP.py` | linux-py3.10 | 12.1 / 0.6 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/unit_models/tests/test_electrodialysis_0D.py` | linux-py3.10 | 1.1 / 10.3 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/flowsheets/activated_sludge/tests/test_asm2_noPHA_flowsheet.py` | linux-py3.10 | 11.1 / 0.3 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/flowsheets/activated_sludge/tests/test_asm2d_flowsheet.py` | linux-py3.10 | 9.4 / 0.4 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/unit_models/tests/test_osmotically_assisted_reverse_osmosis_0D.py` | linux-py3.10 | 0.6 / 7.8 | Shared harness rebuilds an unused model in two checks; harness change applies. |
| `watertap/flowsheets/activated_sludge/tests/test_modified_asm2d_flowsheet.py` | linux-py3.10 | 7.0 / 1.3 | Class fixture already shares model setup; retain distinct initialization/solve checks. |
| `watertap/property_models/tests/test_multicomp_aq_sol_prop_pack.py` | win64-py3.12 | 0.3 / 7.6 | Cost is in case calls; no safe reusable setup shown. |
| `watertap/flowsheets/seawater_RO_desalination/tests/test_seawater_RO_desalination.py` | linux-py3.10 | 0.0 / 7.3 | Distinct end-to-end run(s); no reusable setup shown. |
| `watertap/flowsheets/nf_dspmde/tests/test_nf_with_bypass.py` | linux-py3.10 | 0.0 / 6.3 | Distinct end-to-end run(s); no reusable setup shown. |
| `watertap/unit_models/tests/test_ion_exchange_0D.py` | win64-py3.12 | 0.4 / 5.7 | Shared harness rebuilds an unused model in two checks; harness change applies. |
| `watertap/flowsheets/nf_dspmde/tests/test_nf.py` | linux-py3.10 | 0.0 / 6.0 | Distinct end-to-end run(s); no reusable setup shown. |
| `watertap/unit_models/tests/test_crystallizer.py` | linux-py3.10 | 0.3 / 4.8 | Shared harness rebuilds an unused model in two checks; harness change applies. |

The common `UnitTestHarness` change removes two unused model builds per harness class. All other cases, assertions, solver checks, tolerances, and failure paths remain on PRs. Files without a measured, repeatable improvement remain unchanged.

## Local before and after check

On the same macOS `watertap-dev` environment, the full membrane-distillation 1D test file passed all 96 tests before and after the shared harness edit. The JUnit manifests have identical test IDs and outcomes (96 passed). Pytest elapsed time went from 136.064 to 123.359 seconds; process wall time went from 144.82 to 126.26 seconds. The two affected checks went from 5.944 to 0.001 seconds for conservation and 17.852 to 11.365 seconds for unit solutions. These are one paired run, so they establish an observed saving and preserve outcomes but do not establish its run-to-run distribution. Initialization remained the largest cost, and its checks were retained.

The GAC and reverse-osmosis 1D files also passed together (99 tests) after the change. They exercise configured tolerances, condition-number checks, and scaling exceptions stored by the harness fixture.

The saved eight-cell union contains 3,562 test IDs. Before the matrix change, the five proposed PR cells cover 3,557 IDs; the five missing IDs all belong to LSRRO paper analysis. With `WATERTAP_LSRRO_ALL_CASES=1`, local collection finds all 15 paper-analysis cases, including those five. The saved CSV reports were produced on the prior commit and cannot establish new CI outcomes or the 30-minute job target.

A same-commit, same-environment coverage comparison using the crystallizer test file ran the existing `--cov` command and the new `--cov=watertap` command. Both passed the same 18 tests and recorded the identical 2,556 executed WaterTAP lines in 78 files. The scoped command also lists unexecuted WaterTAP modules in its denominator, so its whole-package percentage differs for this partial suite. One run per command does not establish a coverage-job speedup.
