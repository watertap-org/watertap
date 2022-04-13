###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
"""
Utility script to document schema files.

Example usage and expected output::

    > python document_schemas.py my_output_dir/html/

    > ls html


    Directory: C:\\Users\\MyName\\my_output_dir\\html


    Mode                 LastWriteTime         Length Name
    ----                 -------------         ------ ----
    -a----         5/14/2021   9:00 AM          27345 component.html
    -a----         5/14/2021   9:00 AM           1324 component.json
    -a----         5/14/2021   9:00 AM          18781 reaction.html
    -a----         5/14/2021   9:00 AM           1034 reaction.json
    -a----         5/14/2021   9:00 AM           6391 schema_doc.css
    -a----         5/14/2021   9:00 AM            984 schema_doc.min.js

"""
# stdlib
import argparse
from json_schema_for_humans.generate import generate_from_file_object
import json
from pathlib import Path
import sys

# package
from watertap.edb.schemas import schemas

__author__ = "Dan Gunter (LBNL)"


def main():
    prs = argparse.ArgumentParser(description="Generate schema docs")
    prs.add_argument("directory", help="Directory to put generated schema docs")
    args = prs.parse_args()
    output_dir = Path(args.directory)
    for schema in "component", "reaction":
        schema_file = output_dir / f"{schema}.json"
        with schema_file.open("w") as f:
            json.dump(schemas[schema], f)
        output_file = (output_dir / f"{schema}.html").open("w")
        generate_from_file_object(schema_file.open("r"), output_file)
        print(f"Docs for {schema} at: {output_file.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
