How to use the multicomponent aqueous solution (MCAS) property model
--------------------------------------------------------------------

The example below shows how to use the MCAS property model and display outputs for a state block. Property models allow
users to model the chemical and physical properties of simple systems without the use of unit models.

.. testsetup::

   # quiet idaes logs
   import idaes.logger as idaeslogger
   idaeslogger.getLogger('ideas.core').setLevel('CRITICAL')
   idaeslogger.getLogger('idaes.init').setLevel('CRITICAL')

.. testcode::

    # Import ConcreteModel from Pyomo
    from pyomo.environ import ConcreteModel, assert_optimal_termination
    # Import flowsheet block from IDAES core
    from idaes.core import FlowsheetBlock
    # Import solver from IDAES core
    from idaes.core.solvers import get_solver
    # Import MCAS property model
    import watertap.property_models.multicomp_aq_sol_prop_pack as props
    # Import utility tool for calculating scaling factors
    import idaes.core.util.scaling as iscale

    # Create a concrete model, flowsheet, and NaCl property parameter block.
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.MCASParameterBlock(solute_list=["Na_+", "Cl_-"],
                                               mw_data={"Na_+": 23e-3, 
                                                        "Cl_-": 35e-3},
                                               charge={"Na_+": 1, 
                                                       "Cl_-": -1})

    # Build the state block and specify a time (0 = steady state).
    m.fs.state_block = m.fs.properties.build_state_block([0])
    
    # Specify the state variables of the stream.
    m.fs.state_block[0].flow_mol_phase_comp["Liq", "Cl_-"].fix(0.483)
    m.fs.state_block[0].flow_mol_phase_comp["Liq", "Na_+"].fix(0.483)
    m.fs.state_block[0].flow_mol_phase_comp["Liq", "H2O"].fix(53.8)
    m.fs.state_block[0].pressure.fix(50e5)
    m.fs.state_block[0].temperature.fix(298.15)

    # Set scaling factors for component molar flowrates (variable * scaling factor should be between 0.01 and 100).
    m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1, index=("Liq", "Cl_-")
    )
    
    iscale.calculate_scaling_factors(m)

    # "Touch" build-on-demand variables so that they are created. If these are not touched before running the solver, the output would only display their initial values, not their actual values.
    m.fs.state_block[0].flow_mass_phase_comp
    m.fs.state_block[0].conc_mass_phase_comp
    m.fs.state_block[0].flow_vol_phase
    m.fs.state_block[0].molality_phase_comp
    m.fs.state_block[0].conc_mass_phase_comp
    m.fs.state_block[0].total_hardness
    m.fs.state_block[0].ionic_strength_molal
    m.fs.state_block[0].pressure_osm_phase
    
    # Create the solver object.
    solver = get_solver()

    # Solve the model.
    results = solver.solve(m)

    # Assert that the solution is "optimal"
    assert_optimal_termination(results)

    # Display the results.
    m.fs.state_block[0].display()

A portion of the displayed output is shown below.

.. testoutput::

   Block fs.state_block[0]

     Variables:
       temperature : State temperature
           Size=1, Index=None, Units=K
           Key  : Lower  : Value  : Upper  : Fixed : Stale : Domain
           None : 273.15 : 298.15 : 373.15 :  True :  True : NonNegativeReals
       pressure : State pressure
           Size=1, Index=None, Units=Pa
           Key  : Lower    : Value     : Upper : Fixed : Stale : Domain
           None : 100000.0 : 5000000.0 :  None :  True :  True : NonNegativeReals
       flow_mol_phase_comp : Component molar flow rate
           Size=3, Index=fs.state_block[0].flow_mol_phase_comp_index, Units=mol/s
           Key             : Lower : Value             : Upper : Fixed : Stale : Domain
           ('Liq', 'Cl_-') :     0 : 0.482999971880201 :  None :  True :  True : NonNegativeReals
            ('Liq', 'H2O') :     0 :              53.8 :  None :  True :  True : NonNegativeReals
           ('Liq', 'Na_+') :     0 :             0.483 :  None :  True :  True : NonNegativeReals
       conc_mol_phase_comp : Molar concentration
           Size=3, Index=fs.state_block[0].conc_mol_phase_comp_index, Units=mol/m**3
           Key             : Lower : Value              : Upper : Fixed : Stale : Domain
           ('Liq', 'Cl_-') :     0 :   484.738243699182 :  None : False : False :  Reals
            ('Liq', 'H2O') :     0 :  53993.62117868692 :  None : False : False :  Reals
           ('Liq', 'Na_+') :     0 : 484.73827192018166 :  None : False : False :  Reals
       conc_mass_phase_comp : Mass concentration
           Size=3, Index=fs.state_block[0].conc_mass_phase_comp_index, Units=kg/m**3
           Key             : Lower : Value              : Upper  : Fixed : Stale : Domain
           ('Liq', 'Cl_-') :     0 : 16.965838529471377 : 2000.0 : False : False :  Reals
            ('Liq', 'H2O') :     0 :  971.8851812163646 : 2000.0 : False : False :  Reals
           ('Liq', 'Na_+') :     0 :  11.14898025416418 : 2000.0 : False : False :  Reals
       dens_mass_phase : Mass density
           Size=1, Index=fs.state_block[0].dens_mass_phase_index, Units=kg/m**3
           Key : Lower : Value  : Upper  : Fixed : Stale : Domain
           Liq : 500.0 : 1000.0 : 2000.0 : False : False :  Reals
       mass_frac_phase_comp : Mass fraction
           Size=3, Index=fs.state_block[0].mass_frac_phase_comp_index, Units=dimensionless
           Key             : Lower : Value                : Upper : Fixed : Stale : Domain
           ('Liq', 'Cl_-') :     0 : 0.016965838529471377 : 1.001 : False : False :  Reals
            ('Liq', 'H2O') :     0 :   0.9718851812163645 : 1.001 : False : False :  Reals
           ('Liq', 'Na_+') :     0 :  0.01114898025416418 : 1.001 : False : False :  Reals
       flow_mass_phase_comp : Component Mass flowrate
           Size=3, Index=fs.state_block[0].flow_mass_phase_comp_index, Units=kg/s
           Key             : Lower : Value                : Upper : Fixed : Stale : Domain
           ('Liq', 'Cl_-') :     0 : 0.016904999015807037 :  None : False : False :  Reals
            ('Liq', 'H2O') :     0 :   0.9683999999999999 :  None : False : False :  Reals
           ('Liq', 'Na_+') :     0 : 0.011108999999999999 :  None : False : False :  Reals
       flow_vol_phase : Volumetric flow rate
           Size=1, Index=fs.properties.phase_list, Units=m**3/s
           Key : Lower : Value                 : Upper : Fixed : Stale : Domain
           Liq :     0 : 0.0009964139990158074 :  None : False : False :  Reals
       molality_phase_comp : Molality
           Size=2, Index=fs.state_block[0].molality_phase_comp_index, Units=mol/kg
           Key             : Lower : Value               : Upper : Fixed : Stale : Domain
           ('Liq', 'Cl_-') :     0 : 0.49876081358963353 :  None : False : False :  Reals
           ('Liq', 'Na_+') :     0 :  0.4987608426270136 :  None : False : False :  Reals
       total_hardness : total hardness as CaCO3
           Size=1, Index=None, Units=mg/l
           Key  : Lower : Value : Upper : Fixed : Stale : Domain
           None :     0 :     0 :  None :  True :  True : NonNegativeReals
       ionic_strength_molal : Molal ionic strength
           Size=1, Index=None, Units=mol/kg
           Key  : Lower : Value               : Upper : Fixed : Stale : Domain
           None :     0 : 0.49876082810832356 :  None : False : False : NonNegativeReals
       pressure_osm_phase : van't Hoff Osmotic pressure
           Size=1, Index=fs.properties.phase_list, Units=Pa
           Key : Lower : Value             : Upper : Fixed : Stale : Domain
           Liq :     0 : 2403290.623384758 :  None : False : False :  Reals
       ...


