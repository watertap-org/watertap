.. _mcas_how_to:

How to use the multicomponent aqueous solution (MCAS) property model
--------------------------------------------------------------------

The example below shows how to use the MCAS property model and display outputs for a state block. Property models allow
users to model the chemical and physical properties of simple systems without the use of unit models. For more technical details on the MCAS property model, see the associated :ref:`technical reference<mcas_tech_ref>`.

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

    # Create a concrete model and flowsheet.
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Create an instance of the MCAS property model.
    m.fs.properties = props.MCASParameterBlock(solute_list=["Na_+", "Cl_-"],
                                               mw_data={"Na_+": 23e-3, 
                                                        "Cl_-": 35e-3},
                                               charge={"Na_+": 1, 
                                                       "Cl_-": -1})

    # Build the state block and specify a time (0 = steady state).
    m.fs.state_block = m.fs.properties.build_state_block([0])
    
    # Specify the state variables of the stream (i.e., component molar flowrates, pressure, and temperature).
    m.fs.state_block[0].flow_mol_phase_comp["Liq", "Cl_-"].fix(0.483)
    m.fs.state_block[0].flow_mol_phase_comp["Liq", "Na_+"].fix(0.483)
    m.fs.state_block[0].flow_mol_phase_comp["Liq", "H2O"].fix(53.8)
    m.fs.state_block[0].pressure.fix(50e5)
    m.fs.state_block[0].temperature.fix(298.15)

    # Set scaling factors for component molar flowrates (variable * scaling factor should be between 0.01 and 100).
    m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "Na_+"))
    m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "Cl_-"))
    
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
           Key             : Lower : Value : Upper : Fixed : Stale : Domain
           ('Liq', 'Cl_-') :     0 : 0.483 :  None :  True :  True : NonNegativeReals
            ('Liq', 'H2O') :     0 :  53.8 :  None :  True :  True : NonNegativeReals
           ('Liq', 'Na_+') :     0 : 0.483 :  None :  True :  True : NonNegativeReals
       flow_mass_phase_comp : Component Mass flowrate
           Size=3, Index=fs.state_block[0].flow_mass_phase_comp_index, Units=kg/s
           Key             : Lower : Value                : Upper : Fixed : Stale : Domain
           ('Liq', 'Cl_-') :     0 :             0.016905 :  None : False : False :  Reals
            ('Liq', 'H2O') :     0 :   0.9683999999999999 :  None : False : False :  Reals
           ('Liq', 'Na_+') :     0 : 0.011108999999999999 :  None : False : False :  Reals
       ...

The default material flow basis (i.e., state variable) for the MCAS property model is component molar flowrate. 
However, the user can select component mass flowrate as the flow basis instead as follows.

.. testcode::

    # Import MaterialFlowBasis from the MCAS property model
    from watertap.property_models.multicomp_aq_sol_prop_pack import MaterialFlowBasis

    # Create a concrete model and flowsheet.
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Create an instance of the MCAS property model and use `material_flow_basis` argument to specify mass flowrate as the desired flow basis.
    m.fs.properties = props.MCASParameterBlock(solute_list=["Na_+", "Cl_-"],
                                               mw_data={"Na_+": 23e-3, 
                                                        "Cl_-": 35e-3},
                                               charge={"Na_+": 1, 
                                                       "Cl_-": -1},
                                               material_flow_basis=MaterialFlowBasis.mass)

    # Build the state block and specify a time (0 = steady state).
    m.fs.state_block = m.fs.properties.build_state_block([0])
    
    # Specify the state variables of the stream. Note, now we specify mass flowrate (`flow_mass_phase_comp`) instead of molar flowrate (`flow_mol_phase_comp`).
    m.fs.state_block[0].flow_mass_phase_comp["Liq", "Cl_-"].fix(0.0169)
    m.fs.state_block[0].flow_mass_phase_comp["Liq", "Na_+"].fix(0.0111)
    m.fs.state_block[0].flow_mass_phase_comp["Liq", "H2O"].fix(0.9684)
    m.fs.state_block[0].pressure.fix(50e5)
    m.fs.state_block[0].temperature.fix(298.15)

    # Set scaling factors for component mass flowrates (variable * scaling factor should be between 0.01 and 100).
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 10, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "Cl_-")
    )

    iscale.calculate_scaling_factors(m)
    
``assert_electroneutrality`` is an available method in MCAS. This can be used to assert and optionally adjust composition to enforce electroneutrality.
For a defined composition, i.e., the inlet composition, which is assumed to be known, set ``defined_state`` to True. To adjust composition to enforce electroneutrality, select the ion to adjust with the ``adjust_by_ion`` argument.
 
.. testcode::
  
    m.fs.state_block[0].assert_electroneutrality(defined_state=True, adjust_by_ion="Cl_-")

Output similar to what is shown below will appear to notify the user whether an adjustment in ion composition was made and by how much:

.. testoutput::

   Cl_- adjusted: fs.state_block[0].flow_mass_phase_comp['Liq',Cl_-] was adjusted from 0.0169 and fixed to 0.01689130427193779. Electroneutrality satisfied for fs.state_block[0]. Balance Result = 0.0