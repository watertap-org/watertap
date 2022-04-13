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
    This file demonstrates how to connect to the public MongoDB Cloud server.

    NOTE: To use the public EDB database server, user's MUST have pymongo installed
    with the optional srv service. If not installed, please run the following from
    your [conda] environment.

    python -m pip install pymongo[srv]

    After you connect to the Cloud EDB, you can follow the other examples to build
    your IDAES config from that database. In the unit test, we use the cloud database
    to build a correct IDAES config for the same basic water example shown in file
    'the_basics'.

"""

# Import ElectrolyteDB object
from watertap.edb import ElectrolyteDB

__author__ = "Austin Ladshaw"


def connect_to_cloud_edb(url, db_name, check_connection=False):
    print("connecting to " + url)
    db = ElectrolyteDB(url=url, db=db_name, check_connection=check_connection)
    connected = db.can_connect(url=url, db=db_name)
    return (db, connected)
