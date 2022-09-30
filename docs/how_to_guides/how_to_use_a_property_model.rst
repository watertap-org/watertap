How to use a property model
------------------------------------------------

The example below shows how to use a property model and display outputs for a state block. Property models allow
users to model the chemical and physical properties of simple systems without the use of unit models.

.. testsetup::

   # quiet idaes logs
   import idaes.logger as idaeslogger
   idaeslogger.getLogger('ideas.core').setLevel('CRITICAL')
   idaeslogger.getLogger('idaes.init').setLevel('CRITICAL')

.. testcode::

    # Import concrete model from Pyomo
    from pyomo.environ import ConcreteModel
    # Import flowsheet block from IDAES core
    from idaes.core import FlowsheetBlock
    # Import solver from IDAES core
    from idaes.core.solvers import get_solver
    # Import NaCl property model
    import watertap.property_models.NaCl_prop_pack as props
    # Import utility tool for calculating scaling factors
    import idaes.core.util.scaling as iscale

    # Create a concrete model, flowsheet, and NaCl property parameter block.
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()

    # Build the state block and specify a time (0 = steady state).
    m.fs.state_block = m.fs.properties.build_state_block([0], default={})

    # Fully specify the system.
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.035
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    feed_pressure = 50e5
    feed_temperature = 298.15

    m.fs.state_block[0].flow_mass_phase_comp['Liq', 'NaCl'].fix(feed_flow_mass * feed_mass_frac_NaCl)
    m.fs.state_block[0].flow_mass_phase_comp['Liq', 'H2O'].fix(feed_flow_mass * feed_mass_frac_H2O)
    m.fs.state_block[0].pressure.fix(feed_pressure)
    m.fs.state_block[0].temperature.fix(feed_temperature)

    # Set scaling factors for component mass flowrates (variable * scaling factor should be between 0.01 and 100).
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
    iscale.calculate_scaling_factors(m.fs)

    # "Touch" build-on-demand variables so that they are created. If these are not touched before running the solver, the output would only display their initial values, not their actual values.
    m.fs.state_block[0].dens_mass_phase['Liq']
    m.fs.state_block[0].conc_mass_phase_comp['Liq', 'NaCl']
    m.fs.state_block[0].flow_vol_phase['Liq']
    m.fs.state_block[0].molality_phase_comp['Liq', 'NaCl']
    m.fs.state_block[0].visc_d_phase['Liq']
    m.fs.state_block[0].diffus_phase_comp['Liq', 'NaCl']
    m.fs.state_block[0].enth_mass_phase['Liq']
    m.fs.state_block[0].pressure_osm_phase['Liq']

    # Create the solver object.
    solver = get_solver()

    # Solve the model and display the output.
    solver.solve(m, tee=False)
    m.fs.state_block[0].display()

The following output is displayed:

.. testoutput::

   Block fs.state_block[0]

     Variables:
       flow_mass_phase_comp : Mass flow rate
           Size=2, Index=fs.state_block[0].flow_mass_phase_comp_index, Units=kg/s
           Key             : Lower : Value : Upper : Fixed : Stale : Domain
            ('Liq', 'H2O') :   0.0 : 0.965 :  None :  True :  True : NonNegativeReals
           ('Liq', 'NaCl') :   0.0 : 0.035 :  None :  True :  True : NonNegativeReals
       temperature : State temperature
           Size=1, Index=None, Units=K
           Key  : Lower  : Value  : Upper  : Fixed : Stale : Domain
           None : 273.15 : 298.15 : 373.15 :  True :  True : NonNegativeReals
       pressure : State pressure
           Size=1, Index=None, Units=Pa
           Key  : Lower   : Value     : Upper      : Fixed : Stale : Domain
           None : 10000.0 : 5000000.0 : 50000000.0 :  True :  True : NonNegativeReals
       dens_mass_phase : Mass density
           Size=1, Index=fs.properties.phase_list, Units=kg/m**3
           Key : Lower : Value   : Upper  : Fixed : Stale : Domain
           Liq : 500.0 : 1021.46 : 2000.0 : False : False :  Reals
       mass_frac_phase_comp : Mass fraction
           Size=2, Index=fs.state_block[0].mass_frac_phase_comp_index
           Key             : Lower : Value : Upper : Fixed : Stale : Domain
            ('Liq', 'H2O') :   0.0 : 0.965 :  None : False : False :  Reals
           ('Liq', 'NaCl') :   0.0 : 0.035 :  None : False : False :  Reals
       conc_mass_phase_comp : Mass concentration
           Size=2, Index=fs.state_block[0].conc_mass_phase_comp_index, Units=kg/m**3
           Key             : Lower : Value             : Upper  : Fixed : Stale : Domain
            ('Liq', 'H2O') : 0.001 :          985.7089 : 2000.0 : False : False :  Reals
           ('Liq', 'NaCl') : 0.001 : 35.75110000000001 : 2000.0 : False : False :  Reals
       flow_vol_phase : Volumetric flow rate
           Size=1, Index=fs.properties.phase_list, Units=m**3/s
           Key : Lower : Value                 : Upper : Fixed : Stale : Domain
           Liq :   0.0 : 0.0009789908562254028 :  None : False : False :  Reals
       molality_phase_comp : Molality
           Size=1, Index=fs.state_block[0].molality_phase_comp_index, Units=mol/kg
           Key             : Lower  : Value              : Upper : Fixed : Stale : Domain
           ('Liq', 'NaCl') : 0.0001 : 0.6206267976011888 :    10 : False : False :  Reals
       visc_d_phase : Viscosity
           Size=1, Index=fs.properties.phase_list, Units=Pa*s
           Key : Lower  : Value      : Upper : Fixed : Stale : Domain
           Liq : 0.0001 : 0.00105525 :  0.01 : False : False :  Reals
       diffus_phase_comp : Diffusivity
           Size=1, Index=fs.state_block[0].diffus_phase_comp_index, Units=m**2/s
           Key             : Lower : Value              : Upper : Fixed : Stale : Domain
           ('Liq', 'NaCl') : 1e-10 : 1.471871345625e-09 : 1e-08 : False : False :  Reals
       enth_mass_phase : Specific enthalpy
           Size=1, Index=fs.properties.phase_list, Units=J/kg
           Key : Lower   : Value             : Upper     : Fixed : Stale : Domain
           Liq : 10000.0 : 99740.72571999999 : 1000000.0 : False : False :  Reals
       pressure_osm_phase : Osmotic pressure
           Size=1, Index=fs.properties.phase_list, Units=Pa
           Key : Lower : Value              : Upper      : Fixed : Stale : Domain
           Liq : 500.0 : 2852818.4460273827 : 50000000.0 : False : False :  Reals
       osm_coeff : Osmotic coefficient
           Size=1, Index=None
           Key  : Lower : Value              : Upper : Fixed : Stale : Domain
           None :   0.5 : 0.9271385000000001 :     2 : False : False :  Reals

     Objectives:
       None

     Constraints:
       eq_dens_mass_phase : Size=1
           Key  : Lower : Body : Upper
           None :   0.0 :  0.0 :   0.0
       eq_mass_frac_phase_comp : Size=2
           Key  : Lower : Body : Upper
            H2O :   0.0 :  0.0 :   0.0
           NaCl :   0.0 :  0.0 :   0.0
       eq_conc_mass_phase_comp : Size=2
           Key  : Lower : Body : Upper
            H2O :   0.0 :  0.0 :   0.0
           NaCl :   0.0 :  0.0 :   0.0
       eq_flow_vol_phase : Size=1
           Key  : Lower : Body : Upper
           None :   0.0 :  0.0 :   0.0
       eq_molality_phase_comp : Size=1
           Key  : Lower : Body                    : Upper
           NaCl :   0.0 : -1.1102230246251565e-16 :   0.0
       eq_visc_d_phase : Size=1
           Key  : Lower : Body : Upper
           None :   0.0 :  0.0 :   0.0
       eq_diffus_phase_comp : Size=1
           Key  : Lower : Body : Upper
           NaCl :   0.0 :  0.0 :   0.0
       eq_enth_mass_phase : Size=1
           Key  : Lower : Body : Upper
           None :   0.0 :  0.0 :   0.0
       eq_pressure_osm_phase : Size=1
           Key : Lower : Body : Upper
           Liq :   0.0 :  0.0 :   0.0
       eq_osm_coeff : Size=1
           Key  : Lower : Body : Upper
           None :   0.0 :  0.0 :   0.0
