from pyomo.environ import value, units as pyunits


class CCROConfiguration(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.default_configuration()

    def default_configuration(self):
        self.default_config = dict(
            raw_feed_conc=5 * pyunits.g / pyunits.liter,  # g/L
            raw_feed_flowrate=1 * pyunits.L / pyunits.s,  # L/min
            recycle_flowrate=10 * pyunits.L / pyunits.s,  # L/min
            flushing_conc="raw_feed_conc",  # g/L
            flushing_flowrate="recycle_flowrate",  # L/min
            flushing_efficiency=0.9 * pyunits.dimensionless,
            temperature=298 * pyunits.K,  # K
            p1_pressure_start="osmotic_pressure",  # psi
            osmotic_overpressure=3 * pyunits.dimensionless,
            p1_eff=0.8 * pyunits.dimensionless,
            p2_eff=0.8 * pyunits.dimensionless,
            A_comp=5.96e-12 * pyunits.meter / (pyunits.second * pyunits.Pa),
            B_comp=3.08e-08 * pyunits.meter / pyunits.second,
            membrane_area=100 * pyunits.meter**2,  # m2
            membrane_length="auto",  # m2
            channel_height=0.001 * pyunits.meter,  # m
            spacer_porosity=0.9 * pyunits.dimensionless,
            dead_volume="base_on_dead_volume_to_area_ratio",  # m3
            accumulation_time=5 * pyunits.second,
            dead_volume_to_area_ratio=value(1.2 * 0.001 * pyunits.meter * 0.9)
            * pyunits.m,
            pipe_to_module_ratio=0.2
            * pyunits.dimensionless,  # fraction of dead volume in pipes/pumps if we have hold up in RO
            # Total dead volume is module_area*dead_volume_to_area_ratio + module_area*dead_volume_to_area_ratio*pipe_to_module_ratio
        )
        self.update(self.default_config)
        self.display()

    def display(self):
        """Display current configuration."""
        print("--- CCRO Configuration ---")
        for key, item in self.items():
            print(f"{key}: {value(item)}, {pyunits.get_units(item)}")
        print("--------------------------")

    def __setitem__(self, key, value):
        """Override to prevent changing of existing keys."""
        if isinstance(value, (int, float)) and key not in self.default_config:
            raise ValueError(f"Cannot set {key} without defined pyomo units")
        else:
            if key not in self.default_config:
                raise KeyError(f"{key} is not a valid configuration key")
        super().__setitem__(key, value)
