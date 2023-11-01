import pandas as pd


def main():

    dens_mass = 1000.0
    visc_d = 0.001
    velocity_diluate = 0.01
    channel_height = 2.7e-4
    cell_width = 0.1
    spacer_porosity = 0.83
    spacer_specific_area = 10700
    # diffus_mass = 1.6e-9
    diffus_mass = 1.6e-16
    cell_length = 0.79
    pressure_drop = 1e5

    hydraulic_diameter_list = ["conventional", "specific", "fixed"]
    friction_factor_list = ["Guerri", "Kuroda", "fixed"]
    pressure_drop_list = ["None", "Experimental", "Darcy"]

    method = {}

    def cal_pressure_drop(method):
        friction_factor = cal_friction_factor(method)
        hydraulic_diameter, N_Re = cal_hydraulic_diamter(method)
        if method["pressure_drop_method"] == "None":
            pressure_drop_total = 0
        elif method["pressure_drop_method"] == "Experimental":
            pressure_drop_total = pressure_drop * cell_length
        elif method["pressure_drop_method"] == "Darcy":
            pressure_drop_total = (
                dens_mass
                * friction_factor
                * velocity_diluate**2
                * 0.5
                * hydraulic_diameter**-1
                * cell_length
            )
        return pressure_drop_total

    def cal_hydraulic_diamter(method):

        if method["hydraulic_diameter_method"] == "conventional":
            hydraulic_diameter = (
                2
                * channel_height
                * cell_width
                * (1 - spacer_porosity)
                * (channel_height + cell_width) ** -1
            )
        elif method["hydraulic_diameter_method"] == "specific":
            hydraulic_diameter = (
                4
                * spacer_porosity
                * (
                    2 * channel_height**-1
                    + (1 - spacer_porosity) * spacer_specific_area
                )
            ) ** -1
        elif method["hydraulic_diameter_method"] == "fixed":
            hydraulic_diameter = 1e-3
        else:
            print("wrong method")
        N_Re = dens_mass * velocity_diluate * hydraulic_diameter * visc_d**-1
        return (hydraulic_diameter, N_Re)

    def cal_friction_factor(method):
        hydraulic_diameter, N_Re = cal_hydraulic_diamter(method)
        if method["friction_factor_method"] == "Guerri":
            friction_factor = 4 * 50.6 * spacer_porosity**-7.06 * N_Re**-1
        elif method["friction_factor_method"] == "Kuroda":
            friction_factor = 4 * 9.6 * spacer_porosity**-1 * N_Re**-0.5
        elif method["friction_factor_method"] == "fixed":
            friction_factor = 10
        else:
            print("wrong")
        return friction_factor

    # result = {}
    # full_name = {}
    for pressure_drop_method in pressure_drop_list:
        for friction_factor_method in friction_factor_list:
            for hydraulic_diameter_method in hydraulic_diameter_list:
                method["hydraulic_diameter_method"] = hydraulic_diameter_method
                method["friction_factor_method"] = friction_factor_method
                method["pressure_drop_method"] = pressure_drop_method
                full_name = (
                    pressure_drop_method
                    + "_"
                    + friction_factor_method
                    + "_"
                    + hydraulic_diameter_method
                )
                pressure_drop_total = cal_pressure_drop(method=method)

                print(full_name, pressure_drop_total)

    # df_method = pd.DataFrame.from_dict( orient='index', columns=['Pressure Drop'])
    # outpath = r"C:\Users\kejia_hu\watertap-dev\watertap"
    # df_method.to_csv("pressure_result.csv")
    # print(df_method)


if __name__ == "__main__":
    m = main()
