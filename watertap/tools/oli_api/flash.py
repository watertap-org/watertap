#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

###############################################################################
#
# OLI Systems, Inc. Copyright © 2022, all rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.
#
# 3. Neither the name of OLI Systems, Inc. nor the names of any contributors to
# the software made available herein may be used to endorse or promote products derived
# from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
# SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
# OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the
# features, functionality or performance of the source code ("Enhancements") to anyone; however,
# if you choose to make your Enhancements available either publicly, or directly to OLI Systems, Inc.,
# without imposing a separate written license agreement for such Enhancements, then you hereby grant
# the following license: a non-exclusive, royalty-free perpetual license to install, use, modify, prepare
# derivative works, incorporate into other computer software, distribute, and sublicense such enhancements
# or derivative works thereof, in binary and source code form.
###############################################################################
__author__ = "Oluwamayowa Amusat, Alexander Dudchenko, Paul Vecchiarelli, Adam Atia"


import logging

import json
from pathlib import Path

from copy import deepcopy
from itertools import product

from watertap.tools.oli_api.util.watertap_to_oli_helper_functions import (
    get_oli_name,
    get_charge,
    get_charge_group,
)
from watertap.tools.oli_api.util.fixed_keys_dict import (
    optional_properties,
    input_unit_set,
    output_unit_set,
)

from numpy import reshape, sqrt

_logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
    "OLIAPI - %(asctime)s - %(levelname)s - %(message)s", "%H:%M:%S"
)
handler.setFormatter(formatter)
_logger.addHandler(handler)
_logger.setLevel(logging.DEBUG)


class Flash:
    """
    A class to execute OLI Cloud flash calculations.

    :param optional_properties: dictionary for optional properties to attach to OLI calls, all True by default
    :param input_unit_set: dictionary for conversions between OLI and Pyomo unit names
    :param output_unit_set: dictionary for preferred output units
    :param relative_inflows: bool switch for surveys - true to add specified value to initial value, false to replace initial value with specified value
    :param debug_level: string defining level of logging activity
    """

    def __init__(
        self,
        optional_properties=optional_properties,
        input_unit_set=input_unit_set,
        output_unit_set=output_unit_set,
        relative_inflows=True,
        debug_level="INFO",
    ):
        self.optional_properties = optional_properties
        self.input_unit_set = input_unit_set
        self.output_unit_set = output_unit_set
        self.relative_inflows = relative_inflows

        if debug_level == "INFO":
            _logger.setLevel(logging.INFO)
        else:
            _logger.setLevel(logging.DEBUG)

    def configure_water_analysis(
        self,
        inflows=None,
        temperature=None,
        pressure=None,
        reconciliation=None,
        electroneutrality=None,
        makeup_ion=None,
        ph=None,
        acid_titrant=None,
        base_titrant=None,
        alkalinity=None,
        alkalinity_ph=None,
        alkalinity_titrant=None,
        tic=None,
        allow_solids=False,
        included_solids=None,
        excluded_solids=None,
        calc_alkalinity=False,
        use_scaling_rigorous=True,
        file_name=None,
    ):
        """
        Configure Water Analysis JSON input.

        :param inflows: dictionary of solutes
        :param temperature: float for temperature in Kelvins
        :param pressure: float for pressure in Pascals
        :param reconciliation: string for method of reconciliation: "EquilCalcOnly" (default), "ReconcilePh", "ReconcilePhAndAlkalinity", or "ReconcilePhAndAlkalinityAndTic"; "ReconcileCo2Gas" not supported currently.
        :param electroneutrality: string for method of electroneutrality calculation: "DominantIon", "ProrateCations", "ProrateAnions", "Prorate", "AutoNACL", or "MakeupIon" are supported
        :param makeup_ion: string for ion to use for electroneutrality balance, if "MakeupIon,
        :param ph: float for pH to reconcile solution to, required for pH based reconciliation
        :param acid_titrant: string for acidification titrant, used in pH based reconciliation
        :param base_titrant: string for basification titrant, used in pH based reconciliation
        :param alkalinity: float for alkalinity to reconcile solution to, required for Alk based reconciliation
        :param alkalinity_ph: float for alkalinity endpoint ph, used in Alk based reconciliation
        :param alkalinity_titrant: string for alkalinity titration species, used in Alk based reconciliation
        :param tic: float for total inorganic carbon concentration to reconcile solution to, required for TIC based reconcilation
        :param allow_solids: bool to enable solid phase formation
        :param included_solids: list of solids to include in analysis
        :param excluded_solids: list of solids to exclude from analysis
        :param calc_alkalinity: bool to calculate alkalinity of solution
        :param use_scaling_rigorous: bool to switch between Rigorous (default) and Estimated scaling computations
        :param file_name: string for file to write, if any
        :param mesh_grid: if True (default) the input array will be combined to generate combination of all possible samples
            if False, the direct values in survey_arrays will be used

        :return json_input: JSON for Water Analysis
        """

        _logger.info("Configuring Water Analysis JSON ...")
        input_list = []

        if not inflows:
            raise RuntimeError("Inflows must be defined for Water Analysis.")

        temp_input = {
            "group": "Properties",
            "name": "Temperature",
            "unit": self.input_unit_set["temperature"]["oli_unit"],
            "value": 273.15,
        }
        if temperature is not None:
            if float(temperature):
                temp_input.update({"value": float(temperature)})
            else:
                raise ValueError(f"Invalid temperature: {temperature}. Expected number")
        input_list.append(temp_input)

        pres_input = {
            "group": "Properties",
            "name": "Pressure",
            "unit": self.input_unit_set["pressure"]["oli_unit"],
            "value": 101325,
        }
        if pressure is not None:
            if float(pressure):
                pres_input.update({"value": float(pressure)})
            else:
                raise ValueError(f"Invalid pressure: {pressure}. Expected number")
        input_list.append(pres_input)

        reconciliation_options = [
            "EquilCalcOnly",
            "ReconcilePh",
            "ReconcilePhAndAlkalinity",
            "ReconcilePhAndAlkalinityAndTic",
        ]
        rec_input = {
            "group": "Calculation Options",
            "name": "CalcType",
            "value": "EquilCalcOnly",
        }
        if reconciliation is not None:
            if reconciliation in reconciliation_options:
                rec_input.update({"value": reconciliation})
            else:
                raise RuntimeError(
                    f"Invalid reconciliation option: {reconciliation}."
                    + f" Use one of {reconciliation_options}"
                )
        else:
            reconciliation = "EquilCalcOnly"

        input_list.append(rec_input)
        additional_req_input = []
        additional_req_args = []
        if "Ph" in reconciliation:
            additional_req_args.append([ph, acid_titrant, base_titrant])
            if not acid_titrant:
                acid_titrant = "HCl"
            if not base_titrant:
                base_titrant = "NaOH"
            additional_req_input.extend(
                [
                    {
                        "group": "Properties",
                        "name": "pH",
                        "value": ph,
                    },
                    {
                        "group": "Calculation Options",
                        "name": "PhAcidTitrant",
                        "value": get_oli_name(acid_titrant),
                    },
                    {
                        "group": "Calculation Options",
                        "name": "PhBaseTitrant",
                        "value": get_oli_name(base_titrant),
                    },
                ]
            )
        if "Alk" in reconciliation:
            additional_req_args.append([alkalinity, alkalinity_ph, alkalinity_titrant])
            if not alkalinity_titrant:
                alkalinity_titrant = "H2SO4"
            if not alkalinity_ph:
                alkalinity_ph = 4.5
                _logger.info("No alkalinity endpoint pH specified. Assuming 4.5.")
                additional_req_input.extend(
                    [
                        {
                            "group": "Properties",
                            "name": "Alkalinity",
                            "unit": self.input_unit_set["alkalinity"]["oli_unit"],
                            "value": alkalinity,
                        },
                        {
                            "group": "Properties",
                            "name": "AlkalinityTitrationEndPointpH",
                            "value": alkalinity_ph,
                        },
                        {
                            "group": "Calculation Options",
                            "name": "AlkalinityPhTitrant",
                            "value": alkalinity_titrant,
                        },
                    ]
                )
        if "Tic" in reconciliation:
            additional_req_args.append([tic])
            additional_req_input.append(
                {
                    "group": "Properties",
                    "name": "TIC",
                    "unit": self.input_unit_set["TIC"]["oli_unit"],
                    "value": tic,
                }
            )
        missing_keys = [arg for arg in additional_req_args if arg is None]
        if missing_keys:
            raise RuntimeError(f"Missing keys for {reconciliation}: {missing_keys}")
        input_list.extend(additional_req_input)

        electroneutrality_options = [
            "DominantIon",
            "ProrateCations",
            "ProrateAnions",
            "Prorate",
            "AutoNACL",
            "MakeupIon",
        ]
        elec_input = {
            "group": "Electroneutrality Options",
            "name": "ElectroNeutralityBalanceType",
            "value": "DominantIon",
        }

        if electroneutrality is not None:
            if electroneutrality in electroneutrality_options:
                elec_input.update({"value": electroneutrality})
            else:
                raise RuntimeError(
                    f"Invalid reconciliation option: {electroneutrality}."
                    + f" Use one of {electroneutrality_options}"
                )
        input_list.append(elec_input)
        if electroneutrality == "MakeupIon":
            if makeup_ion is not None:
                input_list.append(
                    {
                        "group": "Electroneutrality Options",
                        "name": "MakeupIonBaseTag",
                        "value": get_oli_name(makeup_ion),
                    }
                )

        input_list.extend(
            [
                {
                    "group": "Calculation Options",
                    "name": "AllowSolidsToForm",
                    "value": bool(allow_solids),
                },
                {
                    "group": "Calculation Options",
                    "name": "CalcAlkalnity",
                    "value": bool(calc_alkalinity),
                },
            ]
        )
        conc_unit = self.input_unit_set["inflows"]["oli_unit"]
        _logger.info(f"Using {conc_unit} for inflows input")
        for k, v in inflows.items():
            charge = get_charge(k)
            input_list.append(
                {
                    "group": get_charge_group(charge),
                    "name": get_oli_name(k),
                    "unit": conc_unit,
                    "value": v,
                    "charge": charge,
                }
            )

        json_input = self._add_to_json(
            "wateranalysis",
            input_list,
            included_solids,
            excluded_solids,
            use_scaling_rigorous,
            file_name,
        )
        return json_input

    def configure_flash_analysis(
        self,
        inflows=None,
        flash_method=None,
        temperature=None,
        pressure=None,
        calculated_variable=None,
        enthalpy=None,
        vapor_amount=None,
        vapor_fraction=None,
        volume=None,
        ph=None,
        acid_titrant=None,
        base_titrant=None,
        formed_solid=None,
        precipitant_inflow=None,
        included_solids=None,
        excluded_solids=None,
        contact_surface=None,
        flow_type=None,
        diameter=None,
        liq_velocity=None,
        gas_velocity=None,
        rot_velocity=None,
        shear_stress=None,
        roughness=None,
        nonaqueous_visc=None,
        water_cut_inversion=None,
        relative_visc_inversion=None,
        use_scaling_rigorous=True,
        file_name=None,
    ):
        """
        Configure Flash Analysis JSON input.

        :param inflows: dictionary of solutes, of the form {"unit": unit, "values": {solute_name: concentration_value}}; otherwise, {solute_name: concentration_value}
        :param flash_method: string for flash calculation name
        :param temperature: float for temperature in Kelvins
        :param pressure: float for pressure in Pascals
        :param calculated_variable: string for variable to calculate, such as temperature or pressure, used in 'bubblepoint', 'dewpoint', 'vapor-amount', 'vapor-fraction', and 'isochoric' flashes
        :param enthalpy: float for total enthalpy in Joules, used in 'isenthalpic' flash
        :param vapor_amount: float for vapor phase Moles, used in 'vapor-amount' flash
        :param vapor_fraction: float for vapor phase in Mole %, used in 'vapor-fraction' flash
        :param volume: float for total volume in Cubic Meters, used in 'isochoric' flash
        :param ph: float for target pH, used in 'setph' flash
        :param acid_titrant: string for acidification titrant, used in 'setph' flash
        :param base_titrant: string for basification titrant, used in 'setph' flash
        :param formed_solid: string for solid species to precipitate based on inflow sweep, used in 'precipitation-point'
        :param precipitant_inflow: string for inflow species to sweep, used in 'precipitation-point'
        :param included_solids: list of solids to include in analysis
        :param excluded_solids: list of solids to exclude from analysis
        :param contact_surface: string for contact surface metal name
        :param flow_type: string for flow configuration
        :param diameter: float for diameter of surface (i.e., pipe or rotor)
        :param liq_velocity: float for velocity of liquid flow
        :param gas_velocity: float for velocity of vapor flow, used in 'approximateMultiPhaseFlow'
        :param rot_velocity: float for rotational velocity
        :param shear_stress: float for defined shear stress, used in 'definedShearStress'
        :param roughness: float for pipe roughness, used in 'approximateMultiPhaseFlow'
        :param nonaqueous_visc: float for absolute viscosity of nonaqueous phase, used in 'approximateMultiPhaseFlow'
        :param water_cut_inversion: float for water cut at point of dispersion inversion, used in 'approximateMultiPhaseFlow'
        :param relative_visc_inversion: float for maximum relative viscosity of dispersion at inversion, used in 'approximateMultiPhaseFlow'
        :param use_scaling_rigorous: bool to switch between Rigorous (default) and Estimated scaling computations
        :param file_name: string for file to write, if any

        :return json_input: JSON for Water Analysis
        """

        _logger.info(f"Configuring {flash_method} Flash JSON ...")

        if flash_method not in [
            "isothermal",
            "isenthalpic",
            "bubblepoint",
            "dewpoint",
            "vapor-amount",
            "vapor-fraction",
            "isochoric",
            "setph",
            "precipitation-point",
            "corrosion-rates",
        ]:
            raise RuntimeError(
                f"Failed to configure Flash. Invalid method: {flash_method}"
            )

        if not inflows:
            raise RuntimeError("Inflows must be defined for Flash Analysis.")

        input_dict = {}
        temp_input = {
            "unit": self.input_unit_set["temperature"]["oli_unit"],
            "value": 273.15,
        }
        if temperature:
            if float(temperature):
                temp_input.update({"value": float(temperature)})
            else:
                raise ValueError(f"Invalid temperature: {temperature}. Expected number")
        input_dict["temperature"] = temp_input

        pres_input = {
            "unit": self.input_unit_set["pressure"]["oli_unit"],
            "value": 101325,
        }
        if pressure:
            if float(pressure):
                pres_input.update({"value": float(pressure)})
            else:
                raise ValueError(f"Invalid pressure: {pressure}. Expected number")
        input_dict["pressure"] = pres_input

        if flash_method in [
            "bubblepoint",
            "dewpoint",
            "vapor-amount",
            "vapor-fraction",
            "isochoric",
        ]:

            if calculated_variable is not None:
                if calculated_variable not in ["temperature", "pressure"]:
                    raise RuntimeError(
                        f"Invalid input for 'calculated_variable': {calculated_variable}; 'temperature' or 'pressure' supported."
                    )
                _logger.info(
                    f"{flash_method} will calculate {calculated_variable} as its variable"
                )
                input_dict["calculatedVariable"] = calculated_variable
            else:
                raise RuntimeError(
                    f"Missing argument for {flash_method}: 'calculated_variable'"
                )

        if flash_method == "isenthalpic":
            enth_input = {
                "unit": self.input_unit_set["enthalpy"]["oli_unit"],
                "value": None,
            }
            if float(enthalpy):
                enth_input.update({"value": float(enthalpy)})
            else:
                raise ValueError(f"Invalid enthalpy: {enthalpy}. Expected number")
            input_dict["enthalpy"] = enth_input

        if flash_method == "vapor-amount":
            vapor_amount_input = (
                {
                    "unit": input_unit_set["vaporAmountMoles"]["oli_unit"],
                    "value": None,
                },
            )
            if float(vapor_amount):
                vapor_amount_input.update({"value": float(vapor_amount)})
            else:
                raise ValueError(
                    f"Invalid vapor amount: {vapor_amount}. Expected number"
                )
            input_dict["vaporAmountMoles"] = vapor_amount_input

        if flash_method == "vapor-fraction":
            vapor_fraction_amount = (
                {
                    "unit": input_unit_set["vaporMolFrac"]["oli_unit"],
                    "value": None,
                },
            )
            if float(vapor_fraction):
                vapor_fraction_amount.update({"value": float(vapor_fraction)})
            else:
                raise ValueError(
                    f"Invalid vapor fraction: {vapor_fraction}. Expected number"
                )
            input_dict["vaporMolFrac"] = vapor_fraction_amount

        if flash_method == "isochoric":
            volume_input = {
                "unit": input_unit_set["totalVolume"]["oli_unit"],
                "value": None,
            }
            if float(volume):
                volume_input.update({"value": float(volume)})
            else:
                raise ValueError(f"Invalid volume: {volume}. Expected number")
            input_dict["totalVolume"] = volume_input

        if flash_method == "setph":
            ph_input = {
                "targetPH": {
                    "unit": "",
                    "value": None,
                },
            }
            if float(ph):
                ph_input["targetPH"].update({"value": float(ph)})
            else:
                raise ValueError(f"Invalid ph: {ph}. Expected number")
            input_dict["targetPH"] = ph_input
            if not acid_titrant:
                acid_titrant = "HCl"
            input_dict["pHAcidTitrant"] = get_oli_name(acid_titrant)
            if not base_titrant:
                base_titrant = "NaOH"
            input_dict["pHBaseTitrant"] = get_oli_name(base_titrant)

        if flash_method == "precipitation-point":
            missing_args = [
                arg for arg in [formed_solid, precipitant_inflow] if arg is None
            ]
            if missing_args:
                raise RuntimeError(
                    f"Missing argument(s) for {flash_method}: {missing_args}"
                )
            else:
                input_dict.update(
                    {
                        "solidToPrecipitate": formed_solid,
                        "inflowToAdjust": precipitant_inflow,
                    }
                )

        if ("unit" in inflows.keys()) and ("values" in inflows.keys()):
            input_dict["inflows"] = inflows
        # otherwise, assume a solute dictionary with {solute_name: concentration_value}
        else:
            unit = self.input_unit_set["molecularConcentration"]["oli_unit"]
            inflows_oli_names = {get_oli_name(k): v for k, v in inflows.items()}
            input_dict["inflows"] = {"unit": unit, "values": inflows_oli_names}

        if flash_method == "corrosion-rates":
            _logger.info(
                f"Ensure DBS file uses 'AQ' thermodynamic framework to use Corrosion Analyzer"
            )
            input_dict["corrosionParameters"] = self._configure_corrosion(
                contact_surface,
                flow_type,
                diameter,
                liq_velocity,
                gas_velocity,
                rot_velocity,
                shear_stress,
                roughness,
                nonaqueous_visc,
                water_cut_inversion,
                relative_visc_inversion,
            )

        json_input = self._add_to_json(
            flash_method,
            input_dict,
            included_solids,
            excluded_solids,
            use_scaling_rigorous,
            file_name,
        )
        return json_input

    def _configure_corrosion(
        self,
        contact_surface,
        flow_type,
        diameter,
        liq_velocity,
        gas_velocity,
        rot_velocity,
        shear_stress,
        roughness,
        nonaqueous_visc,
        water_cut_inversion,
        relative_visc_inversion,
    ):
        """
        Configure input dict for Corrosion Rates flash.

        :param contact_surface: string for contact surface metal name
        :param flow_type: string for flow configuration
        :param diameter: float for diameter of surface (i.e., pipe or rotor)
        :param liq_velocity: float for velocity of liquid flow
        :param gas_velocity: float for velocity of vapor flow, used in 'approximateMultiPhaseFlow'
        :param rot_velocity: float for rotational velocity
        :param shear_stress: float for defined shear stress, used in 'definedShearStress'
        :param roughness: float for pipe roughness, used in 'approximateMultiPhaseFlow'
        :param nonaqueous_visc: float for absolute viscosity of nonaqueous phase, used in 'approximateMultiPhaseFlow'
        :param water_cut_inversion: float for water cut at point of dispersion inversion, used in 'approximateMultiPhaseFlow'
        :param relative_visc_inversion: float for maximum relative viscosity of dispersion at inversion, used in 'approximateMultiPhaseFlow'

        :return config: dictionary for corrosion analysis parameters
        """

        valid_flow_types = [
            "static",
            "pipeFlow",
            "rotatingDisk",
            "rotatingCylinder",
            "completeAgitation",
            "definedShearStress",
            "approximateMultiPhaseFlow",
        ]
        if flow_type not in valid_flow_types:
            raise RuntimeError(
                f"Invalid flow_type: {flow_type}."
                f"Expected one of {', '.join(t for t in valid_flow_types)}"
            )

        config = {
            "calculationType": "isothermal",
            "corrosionParameters": {
                "contactSurface": contact_surface,
                "flowType": flow_type,
            },
        }

        def _try_float(v):
            try:
                val = float(v)
            except:
                val = None
            return val

        _check_args = lambda args: [arg for arg in args if arg is None]
        _check_floats = lambda args: [arg for arg in args if _try_float(arg) is None]
        if flow_type == "pipeFlow":
            args = [diameter, liq_velocity]
            missing_args = _check_args(args)
            not_floats = _check_floats(args)
        elif flow_type in ["rotatingDisk", "rotatingCylinder"]:
            args = [diameter, rot_velocity]
            missing_args = _check_args(args)
            not_floats = _check_floats(args)
        elif flow_type == "definedShearStress":
            args = [shear_stress]
            missing_args = _check_args(args)
            not_floats = _check_floats(args)
        elif flow_type == "approximateMultiPhaseFlow":
            args = [
                diameter,
                liq_velocity,
                gas_velocity,
                roughness,
                nonaqueous_visc,
                water_cut_inversion,
                relative_visc_inversion,
            ]
            missing_args = _check_args(args)
            not_floats = _check_floats(args)
        if missing_args:
            raise RuntimeError(
                f"Missing argument(s) for {flash_method}: {missing_args}"
            )
        if not_floats:
            raise RuntimeError(
                f"Invalid values for argument(s): {not_floats}. Expected value"
            )

        if flow_type == "pipeFlow":
            config["corrosionParameters"].update(
                {
                    "pipeDiameter": {
                        "value": diameter,
                        "unit": self.input_unit_set["pipeDiameter"]["oli_unit"],
                    },
                    "pipeFlowVelocity": {
                        "value": liq_velocity,
                        "unit": self.input_unit_set["pipeFlowVelocity"]["oli_unit"],
                    },
                }
            )
        elif flow_type == "rotatingDisk":
            config["corrosionParameters"].update(
                {
                    "diskDiameter": {
                        "value": diameter,
                        "unit": self.input_unit_set["diskDiameter"]["oli_unit"],
                    },
                    "diskRotationSpeed": {
                        "value": rot_velocity,
                        "unit": self.input_unit_set["diskRotationSpeed"]["oli_unit"],
                    },
                },
            )

        elif flow_type == "rotatingCylinder":
            config["corrosionParameters"].update(
                {
                    "rotorDiameter": {
                        "value": diameter,
                        "unit": self.input_unit_set["rotorDiameter"]["oli_unit"],
                    },
                    "rotorRotation": {
                        "value": rot_velocity,
                        "unit": self.input_unit_set["rotorRotation"]["oli_unit"],
                    },
                },
            )

        elif flow_type == "definedShearStress":
            config["corrosionParameters"].update(
                {
                    "shearStress": {
                        "value": shear_stress,
                        "unit": self.input_unit_set["shearStress"]["oli_unit"],
                    },
                },
            )

        elif flow_type == "approximateMultiPhaseFlow":
            config["corrosionParameters"].update(
                {
                    "pipeDiameter": {
                        "value": diameter,
                        "unit": self.input_unit_set["pipeDiameter"]["oli_unit"],
                    },
                    "liquidFlowInPipe": {
                        "value": liq_velocity,
                        "unit": self.input_unit_set["liquidFlowInPipe"]["oli_unit"],
                    },
                    "gasFlowInPipe": {
                        "value": gas_velocity,
                        "unit": self.input_unit_set["gasFlowInPipe"]["oli_unit"],
                    },
                    "pipeRoughness": {
                        "value": roughness,
                        "unit": self.input_unit_set["pipeRoughness"]["oli_unit"],
                    },
                    "viscAbs2ndLiq": {
                        "value": nonaqueous_visc,
                        "unit": self.input_unit_set["viscAbs2ndLiq"]["oli_unit"],
                    },
                    "waterCutAtPointOfDispersionInversion": water_cut_inversion,
                    "maxRelViscosityOfDispersionAtInversion": relative_visc_inversion,
                }
            )
        return config

    def _add_to_json(
        self,
        flash_method,
        input_data,
        included_solids,
        excluded_solids,
        use_scaling_rigorous,
        file_name,
    ):
        """
        Add input data to JSON.

        :param flash_method: string for flash calculation name
        :param input_data: data object from flash configuration function
        :param included_solids: list of solids to include in analysis
        :param excluded_solids: list of solids to exclude from analysis
        :param use_scaling_rigorous: bool to switch between Rigorous (default) and Estimated scaling computations
        :param file_name: string for file to write, if any
        """

        self._set_prescaling_calculation_mode(use_scaling_rigorous)
        json_input = {"params": {}}

        if flash_method == "wateranalysis":
            json_input["params"].update({"waterAnalysisInputs": input_data})
        else:
            json_input["params"] = input_data

        output_unit_set_info = {}
        for k, v in self.output_unit_set.items():
            output_unit_set_info[k] = v["oli_unit"]
        additional_params = {
            "optionalProperties": dict(self.optional_properties),
            "unitSetInfo": dict(output_unit_set_info),
        }
        if included_solids and excluded_solids:
            raise RuntimeError(
                "Invalid argument combination. "
                "Only one of included_solids and excluded_solids "
                "may be specified at once."
            )
        else:
            if included_solids:
                additional_params.update({"included_solids": list(included_solids)})
            if excluded_solids:
                additional_params.update({"excluded_solids": list(excluded_solids)})
        json_input["params"].update(additional_params)

        if file_name is not None:
            write_output(json_input, file_name)
        return json_input

    def _set_prescaling_calculation_mode(self, use_scaling_rigorous):
        """
        Set prescaling computation method based on argument.

        :param use_scaling_rigorous: boolean indicating desired state of 'rigorous' and 'estimated' optional properties
        """

        props = self.optional_properties
        if bool(use_scaling_rigorous) == bool(props["prescalingTendenciesRigorous"]):
            return
        new_values = {k: not v for k, v in props.items() if "prescaling" in k}
        props.update(new_values)

    def run_flash(
        self,
        flash_method,
        oliapi_instance,
        dbs_file_id,
        json_input,
        survey=None,
        file_name=None,
        # max_concurrent_processes=1000,
        # burst_job_tag=None,
        # batch_size=None,
    ):
        """
        Conduct single point analysis with initial JSON input, or conduct a survey on that input.

        :param flash_method: string for flash calculation name
        :param oliapi_instance: instance of OLI Cloud API
        :param dbs_file_id: string ID of DBS file
        :param json_input: JSON input for flash calculation
        :param survey: dictionary containing names and input values to modify in JSON
        :param file_name: string for file to write, if any

        :return processed_requests: results from processed OLI flash requests
        """

        if flash_method == "corrosion-rates":
            # check if DBS file is using AQ thermodynamic framework
            oliapi_instance.get_corrosion_contact_surfaces(dbs_file_id)

        if self.relative_inflows:
            _logger.info(
                f"relative_inflows={self.relative_inflows},"
                + " surveys will add values to initial state"
            )

        if survey is None:
            survey = {}
        num_samples = None
        for k, v in survey.items():
            if num_samples is None:
                num_samples = len(v)
            elif num_samples != len(v):
                raise RuntimeError(f"Length of list for key {k} differs from prior key")
        if num_samples is None:
            num_samples = 1
        requests_to_process = []
        for idx in range(num_samples):
            _logger.info(f"Flash sample #{idx+1} of {num_samples}")
            requests_to_process.append(
                {
                    "flash_method": flash_method,
                    "dbs_file_id": dbs_file_id,
                    "input_params": self.get_clone(
                        flash_method, json_input, idx, survey
                    ),
                }
            )
        processed_requests = oliapi_instance.process_request_list(
            requests_to_process,
            # burst_job_tag=burst_job_tag,
            # max_concurrent_processes=max_concurrent_processes,
            # batch_size=batch_size,
        )
        _logger.info("Completed running flash calculations")
        result = flatten_results(processed_requests)
        if file_name:
            write_output(result, file_name)
        return result

    def get_clone(self, flash_method, json_input, index, survey=None):
        """
        Iterate over a survey to create a modified clone from JSON input.

        :param flash_method: string for flash calculation name
        :param json_input: JSON input for flash calculation
        :param index: integer for index of incoming data
        :param survey: dictionary containing names and input values to modify in JSON

        :return clone: dictionary containing modified state variables and survey index
        """

        if survey is None:
            return json_input

        valid_flashes = [
            "wateranalysis",
            "isothermal",
            "isenthalpic",
            "bubblepoint",
            "dewpoint",
            "vapor-amount",
            "vapor-fraction",
            "isochoric",
            "setph",
            "precipitation-point",
            "corrosion-rates",
        ]
        if flash_method not in valid_flashes:
            raise RuntimeError(
                "Invalid flash_method: {flash_method}. Use one of {', '.join(valid_flashes)}"
            )

        clone = deepcopy(json_input)
        for k, v in survey.items():
            d = clone["params"]
            if flash_method == "wateranalysis":
                d = d["waterAnalysisInputs"]
                for param in d:
                    if param["name"].lower() == k.lower():
                        if self.relative_inflows:
                            param["value"] += v[index]
                        else:
                            param["value"] = v[index]
                        _logger.info(
                            f"Updating {k} for sample #{index} clone: new value = {param['value']}"
                        )
            else:
                if k in d:
                    pass
                elif k in d["inflows"]["values"]:
                    d = d["inflows"]["values"]
                elif hasattr(d, "corrosionParameters") and (
                    k in d["corrosionParameters"]
                ):
                    d = d["corrosionParameters"]
                else:
                    _logger.warning(f"Survey key {k} not found in JSON input.")
                if self.relative_inflows:
                    if isinstance(d[k], dict):
                        d[k]["value"] += v[index]
                        val = d[k]["value"]
                    else:
                        d[k] += v[index]
                        val = d[k]
                else:
                    if isinstance(d[k], dict):
                        d[k]["value"] = v[index]
                        val = d[k]["value"]
                    else:
                        d[k] = v[index]
                        val = d[k]
                _logger.info(
                    f"Updating {k} for sample #{index} clone: new value = {val}"
                )
        return clone

    def get_apparent_species_from_true(
        self,
        true_species_json,
        oliapi_instance,
        dbs_file_id,
        phase=None,
        file_name=None,
    ):
        """
        Run Water Analysis to get apparent species.

        :param true_species_json: JSON generated from true species
        :param oliapi_instance: instance of OLI Cloud API to call
        :param dbs_file_id: string ID of DBS file
        :param phase: string for inflows phase
        :param file_name: string for file to write, if any

        :return apparent_species: dictionary for molecular concentrations
        """

        stream_output = self.run_flash(
            "wateranalysis",
            oliapi_instance,
            dbs_file_id,
            true_species_json,
        )
        if phase is None:
            phase = "total"

        extracted_result = stream_output["result"][f"molecularConcentration_{phase}"]

        unit = self.input_unit_set["molecularConcentration"]["oli_unit"]
        concentrations = {k: v["values"][0] for k, v in extracted_result.items()}
        inflows = {"unit": unit, "values": concentrations}
        if file_name:
            write_output(inflows, file_name)
        return inflows


def flatten_results(processed_requests):

    _logger.info("Flattening OLI stream output ... ")

    props = []
    terminal_keys = ["unit", "value", "found", "fullVersion", "values"]

    def _find_props(data, path=None):
        """
        Get the path to all nested items in input data (recursive search).

        :param data: dictionary containing OLI flash output
        :param path: list of paths to endpoint

        :return props: list of nested path lists
        """
        path = path if path is not None else []
        if isinstance(data, dict):
            for k, v in data.items():
                if isinstance(v, (str, bool)):
                    props.append([*path, k])
                elif isinstance(v, list):
                    if all(k not in terminal_keys for k in v):
                        _find_props(v, [*path, k])
                elif isinstance(v, dict):
                    if all(k not in terminal_keys for k in v):
                        _find_props(v, [*path, k])
                    else:
                        props.append([*path, k])
        elif isinstance(data, list):
            for idx, v in enumerate(data):
                if isinstance(v, (dict, list)):
                    if all(k not in terminal_keys for k in v):
                        _find_props(v, [*path, idx])
                    else:
                        props.append([*path, idx])
        else:
            raise RuntimeError(f"Unexpected type for data: {type(data)}")

    def _get_nested_data(data, keys):

        for key in keys:
            data = data[key]
        return data

    def _extract_values(data, keys):
        values = _get_nested_data(data, keys)
        extracted_values = {}
        if isinstance(values, str):
            extracted_values = values
        elif isinstance(values, bool):
            extracted_values = bool(values)
        elif isinstance(values, dict):
            if any(k in values for k in ["group", "name", "fullVersion"]):
                if "value" in values:
                    extracted_values.update({"values": values["value"]})
                if "unit" in values:
                    unit = values["unit"] if values["unit"] else "dimensionless"
                    extracted_values.update({"units": unit})

            elif all(k in values for k in ["found", "phase"]):
                extracted_values = values
            else:
                unit = values["unit"] if values["unit"] else "dimensionless"
                if "value" in values:
                    extracted_values = {
                        "units": unit,
                        "values": values["value"],
                    }
                elif "values" in values:
                    extracted_values = {
                        k: {
                            "units": unit,
                            "values": values["values"][k],
                        }
                        for k, v in values["values"].items()
                    }
                elif "data" in values:
                    # intended for vaporDiffusivityMatrix
                    mat_dim = int(sqrt(len(values["data"])))
                    diffmat = reshape(values["data"], newshape=(mat_dim, mat_dim))

                    extracted_values = {
                        f'({values["speciesNames"][i]},{values["speciesNames"][j]})': {
                            "units": values["unit"],
                            "values": diffmat[i][j],
                        }
                        for i in range(len(diffmat))
                        for j in range(i, len(diffmat))
                    }
                else:
                    raise NotImplementedError(
                        f"results structure not accounted for. results:\n{values}"
                    )
        else:
            raise RuntimeError(f"Unexpected type for data: {type(values)}")
        return extracted_values

    def _create_input_dict(props, result):
        input_dict = {k: {} for k in set([prop[0] for prop in props])}
        for prop in props:
            k = prop[0]
            phase_tag = ""
            if "metaData" in prop:
                prop_tag = prop[-1]
            elif "result" in prop:
                # get property tag
                if isinstance(prop[-1], int):
                    prop_tag = prop[-2]
                else:
                    prop_tag = prop[-1]
                # get phase tag
                if any(k in prop for k in ["phases", "total"]):
                    if "total" in prop:
                        phase_tag = "total"
                    else:
                        phase_tag = prop[prop.index("phases") + 1]
            elif "submitted_requests" in prop:
                prop_tag = prop[-1]
                if "params" in prop:
                    if isinstance(prop[-1], int):
                        prop_tag = _get_nested_data(result, prop)["name"]
            else:
                _logger.warning(
                    f"Unexpected result:\n{result}\n\ninput_dict:\n{input_dict} from prop {prop}"
                )
                continue
            label = f"{prop_tag}_{phase_tag}" if phase_tag else prop_tag
            input_dict[k][label] = _extract_values(result, prop)
        return input_dict

    float_nan = float("nan")

    def _add_to_output(input_dict, output_dict, index, number_samples):
        """
        Add incoming flash results to output data.

        :param input_dict: dictionary for incoming data
        :param output_dict: dictionary for output data
        :param index: integer for index of incoming data
        :param number_samples: integer for total number of incoming data samples
        """

        for k, v in input_dict.items():
            try:
                val = float(v)
            except:
                val = None
            if val is not None:
                if k not in output_dict:
                    output_dict[k] = [float_nan] * number_samples
                output_dict[k][index] = val
            elif isinstance(v, str):
                if k not in output_dict:
                    if k in ["fullVersion", "units"]:
                        output_dict[k] = v
                    else:
                        output_dict[k] = [float_nan] * number_samples
                if k in ["fullVersion", "units"]:
                    if input_dict[k] != output_dict[k]:
                        raise Exception(f"Input and output do not agree for key {k}")
                else:
                    output_dict[k][index] = v
            elif isinstance(v, dict):
                if k not in output_dict:
                    output_dict[k] = {}
                _add_to_output(input_dict[k], output_dict[k], index, number_samples)
            else:
                raise Exception(f"Unexpected value: {v}")

    output_dict = {}
    num_samples = len(processed_requests)
    for idx, result in enumerate(processed_requests):
        props = []
        _find_props(result)
        input_dict = _create_input_dict(props, result)
        _add_to_output(input_dict, output_dict, idx, num_samples)
    return output_dict


def write_output(content, file_name):
    """
    Write dictionary-based content to file.

    :param content: dictionary of content to write
    :param file_name: string for name of file to write

    :param file_path: string for full path of written file
    """

    _logger.info(f"Saving content to {file_name}")
    with open(file_name, "w", encoding="utf-8") as f:
        json.dump(content, f)
    _logger.info("Save complete")
    file_path = Path(file_name).resolve()
    return file_path


def build_survey(survey_arrays, get_oli_names=False, file_name=None, mesh_grid=True):
    """
    Build a dictionary for modifying flash calculation parameters.

    :param survey_arrays: dictionary for variables and values to survey
    :param get_oli_names: bool switch to convert name into OLI name
    :param file_name: string for file to write, if any
    :param mesh_grid: if True (default) the input array will be combined to generate combination of all possible samples
        if False, the direct values in survey_arrays will be used

    :return survey: dictionary for product of survey variables and values
    """
    _name = lambda k: get_oli_name(k) if get_oli_names else k
    if mesh_grid:
        keys = [get_oli_name(k) if get_oli_names else k for k in survey_arrays]
        values = list(product(*(survey_arrays.values())))
        survey = {_name(keys[i]): [val[i] for val in values] for i in range(len(keys))}
    else:
        survey = {}
        values = None
        for key, arr in survey_arrays.items():
            survey[_name(key)] = arr
            if values is not None and len(values) != len(arr):
                raise ValueError(f"Length of list for key {key} differs from prior key")
            values = arr
    _logger.info(f"Survey contains {len(values)} items.")
    if file_name:
        write_output(survey, file_name)
    return survey


def get_survey_sample_conditions(survey, sample_points):
    """
    Return survey parameter values for one or more sample points.

    :param survey: dictionary for product of survey conditions and values
    :param sample_points: list of indices to get parameter values from

    :return sample_conditions: dictionary for parameter values for given samples
    """

    sample_conditions = {}
    for point in sample_points:
        sample_conditions[point] = {}
        for k, v in survey.items():
            sample_conditions[point][k] = v[point]
    _logger.debug(sample_conditions)
    return sample_conditions
