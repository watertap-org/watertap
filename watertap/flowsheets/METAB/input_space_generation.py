#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

# input space generation

from idaes.core.surrogate.pysmo.sampling import LatinHypercubeSampling
import csv

__author__ = "Marcus Holly"


def create_samples(
    method="LHS",
    input_var_info=None,
    sample_numbers=10,
    csv_file="./results/input_data.csv",
):
    """
    Generate input samples over a defined variable space and write them to a CSV file

    Args:
        method (str, optional): Sampling strategy to use.
            This currently only supports the Latin Hypercube Sampling (LHS) method
        input_var_info (dict, optional): Ordered mapping of input variable names
        sample_numbers (int, optional): Number of sample points to generate
        csv_file (str, optional): Path to the output CSV file

    Returns:
        None. Results are written to ``csv_file`` and a confirmation message is
        printed on success.
    """
    if method == None:
        print("Please pick a sampling method")
        return
    elif method == "LHS":
        lower_bounds = [lb[0] for lb in input_var_info.values()]
        upper_bounds = [lb[1] for lb in input_var_info.values()]
        bounds = [lower_bounds, upper_bounds]
        samples = LatinHypercubeSampling(
            bounds, number_of_samples=sample_numbers, sampling_type="creation"
        ).sample_points()
        # Round each number in the matrix to 2 decimal places
        rounded_samples = [[round(sample, 2) for sample in row] for row in samples]
        print(type(rounded_samples))

    # Writing to CSV file
    with open(csv_file, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(input_var_info.keys())
        writer.writerows(rounded_samples)

    print(f'Sample file "{csv_file}" has been created successfully.')


if __name__ == "__main__":

    input_var_info = {
        "inf_fr": (5, 10),  # bounds
        "temp": (22, 35),
        "hrt": (1, 12),
    }

    create_samples(method="LHS", input_var_info=input_var_info, sample_numbers=20)
