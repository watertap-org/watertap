#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

__author__ = "Adam Atia"


def get_property_metadata(prop_pkg):
    """Return metadata for every property a package marks as supported.

    Args:
        prop_pkg: a property model ParameterBlock (e.g. m.fs.properties)

    Returns:
        list of dicts with keys "Description", "Name", "Units", sorted alphabetically by Description
    """
    pset = prop_pkg.get_metadata().properties
    rows = []
    for prop in pset:
        # TODO: switch prop._doc to prop.doc once doc property added in IDAES
        doc = getattr(prop, "doc", None) or prop._doc
        # indices can include None, phase_comp, etc.
        for idx in prop._indices:
            sub = prop[idx]
            if sub.supported:
                rows.append(
                    {"Description": doc, "Name": sub.name, "Units": str(sub.units)}
                )
        sorted_rows = sorted(rows, key=lambda r: r["Description"].lower())

    return sorted_rows


def print_property_metadata(prop_pkg):
    """Pretty-print supported properties as a fixed-width table."""
    rows = get_property_metadata(prop_pkg)
    if not rows:
        print("No supported properties found.")
        return
    cols = ["Description", "Name", "Units"]
    widths = {c: max(len(c), max(len(r[c]) for r in rows)) for c in cols}
    header = "  ".join(c.ljust(widths[c]) for c in cols)
    print(header)
    print("-" * len(header))
    for r in rows:
        print("  ".join(r[c].ljust(widths[c]) for c in cols))
