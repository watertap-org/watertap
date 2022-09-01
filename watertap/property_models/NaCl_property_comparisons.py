import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():
    pitzer_h_file = "C:/Users/carso/Documents/MVC/watertap_results/Property comparisons/Pitzer Specific Enthalpy.csv"
    df_pitzer_h = pd.read_csv(pitzer_h_file)
    print(df_pitzer_h.columns)
    t = np.array(df_pitzer_h[:]['t'])
    P = np.array(df_pitzer_h[:]['P'])
    df_pitzer_h.drop('t', axis=1, inplace=True)
    df_pitzer_h.drop('P', axis=1, inplace=True)
    mass_fractions = []
    plt.figure()
    for col_name in df_pitzer_h.columns:
        w = molality_to_mass_fraction_NaCl(float(col_name))
        mass_fractions.append(w)
        plt.plot(t,df_pitzer_h[col_name],label='w = ' + str(round(w,3)))
    plt.xlabel('Temperature')
    plt.ylabel('Specific enthalpy [J/kg]')
    plt.legend()
    plt.show()


def molality_to_mass_fraction_NaCl(m):
    mol_weight = 1e-3*58.44 # kg/mol for NaCl
    return m*mol_weight/(1+m*mol_weight)

def seawater_properties(mass_fractions):

    return

if __name__ == "__main__":
    main()