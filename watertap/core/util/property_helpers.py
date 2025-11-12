
import pandas as pd

def print_property_metadata(prop_pkg, return_df=False):
    """
    Print all supported properties from a WaterTAP/IDAES property package.
    
    Args:
        prop_pkg: The property model ParameterBlock (e.g., m.fs.properties)
        return_df: If True, returns a Pandas DataFrame instead of printing.
    """
    metadata = prop_pkg.get_metadata()
    vars, units, docs = [], [], []

    for v in metadata.properties:
        vars.append(v._name)
        units.append(str(v._units))
        docs.append(v._doc)

    if return_df:
        return pd.DataFrame({
            "Property Description": docs,
            "Model Attribute": vars,
            "Units": units
        })

    # Pretty-print
    name_col = "Model Attribute"
    desc_col = "Property Description"
    units_col = "Units"

    name_w = max(len(name_col), max(len(n) for n in vars)) + 2
    desc_w = max(len(desc_col), max(len(d) for d in docs)) + 2
    units_w = max(len(units_col), max(len(u) for u in units)) + 2

    print(f"{desc_col:<{desc_w}}{name_col:<{name_w}}{units_col:<{units_w}}")
    print("-" * (name_w + desc_w + units_w))

    for n, d, u in zip(vars, docs, units):
        print(f"{d:<{desc_w}}{n:<{name_w}}{u:<{units_w}}")
