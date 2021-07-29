###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################
'''
    This test is for the tutorials in the ../../tutorials/ directory.
    For now, it merely checks that the tutorials execute without raising
    an exception, it does not check for correctness
'''

import nbformat
import unittest
import glob
import os.path
from nbconvert.preprocessors import ExecutePreprocessor

_tutorials_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','tutorials')
_notebooks = glob.glob(os.path.join(_tutorials_dir, '**', '*.ipynb'), recursive=True)

class TestTutorials(unittest.TestCase):

    def test_tutorials(self):
        for notebook_filename in _notebooks:
            with self.subTest(notebook=os.path.basename(notebook_filename)):
                with open(notebook_filename) as f:
                    nb = nbformat.read(f, as_version=4)
                ep = ExecutePreprocessor()
                ep.preprocess(nb, 
                        {'metadata': {'path':f'{os.path.dirname(notebook_filename)}'}}) 

if __name__ == '__main__':
    unittest.main()
