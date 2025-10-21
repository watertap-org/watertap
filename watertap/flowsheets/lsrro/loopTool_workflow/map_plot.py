from psPlotKit.data_manager.ps_data_manager import PsDataManager
from psPlotKit.data_plotter.ps_map_plotter import MapPlotter
import numpy as np

if __name__ == "__main__":

    ##################################################################################
    # Data processing workflow follows these steps:
    # 1) Load data from h5 file using PsDataManager
    # 2) Register data keys for easier access and labeling
    # 2.a) Register keys you want to import, they will be imported from every simulation in h5 file, you
    #      can also limit to specific simulations using directories option
    # 3) Useing data manager its possible to reduce data from all simulations
    #   In this example, we simulated across 9 different stages, we cane use reduce_data option
    #  to find cost-optimal nubmer of stages, by specifing that we want to stack all data based on "number_of_stages"
    #  directory, and use LCOW as our minium objective, any data_key can be used.
    #  The reduction operation will create a new set of data called "stacked data" that can be plotted
    ###################################################################################
    data_manager = PsDataManager(
        [
            "output/lsrro_stage_sweep_analysisType_map_sweep.h5",
        ]
    )
    # register keys of interest
    # First input is the model key or file key you see in the h5 file,
    # return key is the key that will be returned to you
    # unit is desired final unit for conversion upon import
    # assign_units allows you to assign or force assignment of a unit to data during import, use it with
    # conversion_factor to define conversion from default unit to assigned unit
    #   common use case water flux, which might be saved with out default units, so you can assign
    #   unit of LMH, and conversion of 3600 to go from kg/m2/s to L/m2/h
    data_manager.register_data_key(
        "fs.water_recovery", return_key="Water recovery", units="%"
    )
    data_manager.register_data_key(
        "fs.feed.properties[0].conc_mass_phase_comp[Liq, NaCl]",
        return_key="Feed TDS",
        units="g/L",
    )
    data_manager.register_data_key("fs.costing.LCOW", return_key="LCOW")
    data_manager.load_data()
    data_manager.display()
    ### reduce data to find optimal number of stages for each condition
    data_manager.reduce_data(
        stack_keys="number_of_stages", data_key="LCOW", reduction_type="min"
    )
    data_manager.display()

    # We can select specific data sub sets of data to make it easier to work, and also
    # ensure we don't accidently plot wrond data
    data_manager.select_data("stacked_data", True)
    # generate working data
    wr = data_manager.get_selected_data()
    wr.display()

    # WE can sue MapPlotter tool to create our maps for cost and optimal number of stages
    # This can work really for any variable of interest.

    mp = MapPlotter(wr, save_folder="figures", save_name=f"LCOW_map")
    mp.plot_map(
        "stacked_data",
        ydata="Water recovery",
        xdata="Feed TDS",
        zdata="LCOW",
        zlevels=[0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20],
        axis_options={
            "ylabel": "Water recovery (%)",
            "xlabel": "Feed TDS (g/L)",
            "xticklabels": [5, 50, 100, 150, 200],
            "yticklabels": [30, 40, 50, 60, 70, 80, 90],
            "zlabel": f"LCOW ($\$$/m$^3$)",
            "zticks": [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20],
        },
    )
    mp = MapPlotter(data_manager, save_folder="figures", save_name=f"optimal_stages")
    mp.plot_map(
        "stacked_data",
        ydata="Water recovery",
        xdata="Feed TDS",
        zdata="number_of_stages",
        zlevels=np.arange(0.0, 5.5, 0.5),
        axis_options={
            "ylabel": "Water recovery (%)",
            "xlabel": "Feed TDS (g/L)",
            "xticklabels": [5, 50, 100, 150, 200, 250],
            "yticklabels": [30, 40, 50, 60, 70, 80, 90],
            "yticklabels": [30, 40, 60, 70, 80, 90],
            "zlabel": f"Number of Stages (#)",
            "zticks": [2, 3, 4, 5, 6, 7, 8, 9],
        },
    )
