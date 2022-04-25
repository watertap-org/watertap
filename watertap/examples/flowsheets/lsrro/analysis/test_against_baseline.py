import glob
import os.path as path
import numpy as np
import pandas as pd

from .lsrro_paper_analysis import main

def test_against_baseline():
    main()

    non_numeric_columns = ['A_case','B_case','AB_Tradeoff']
    for fn in glob.glob("./recheck_paper_cases_after_erd_fix_baseline/*.csv"):
        baseline = pd.read_csv(fn)
        test = pd.read_csv("./recheck_paper_cases_after_erd_fix/"+path.basename(fn))
    
        if not (baseline[non_numeric_columns] == test[non_numeric_columns]).all().all():
            raise ValueError("Difference in non-numeric values in "+path.basename(fn))
    
        test = test.T.drop(non_numeric_columns).T.astype(float)
        baseline = baseline.T.drop(non_numeric_columns).T.astype(float)
    
        if not np.isclose(baseline, test, equal_nan=True).all():
            raise ValueError("Difference in numeric values in "+path.basename(fn))
    
    print("All good")
