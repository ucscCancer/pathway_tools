
import unittest
import os
import sys
import subprocess
import shutil
import time

from pathway_tools import dai_util

class TestCase(unittest.TestCase):    
    def test_multidim(self):
        #make sure that the internal table and multi_dim_iter go in the same order
        matrix = [
            0.999000, # [0, 0, 0]
            0.999000, # [1, 0, 0]
            0.999000, # [2, 0, 0]
            0.999000, # [0, 1, 0]
            0.000500, # [1, 1, 0]
            0.000500, # [2, 1, 0]
            0.999000, # [0, 2, 0]
            0.000500, # [1, 2, 0]
            0.000500, # [2, 2, 0]
            0.000500, # [0, 0, 1]
            0.000500, # [1, 0, 1]
            0.000500, # [2, 0, 1]
            0.000500, # [0, 1, 1]
            0.999000, # [1, 1, 1]
            0.000500, # [2, 1, 1]
            0.000500, # [0, 2, 1]
            0.000500, # [1, 2, 1]
            0.000500, # [2, 2, 1]
            0.000500, # [0, 0, 2]
            0.000500, # [1, 0, 2]
            0.000500, # [2, 0, 2]
            0.000500, # [0, 1, 2]
            0.000500, # [1, 1, 2]
            0.999000, # [2, 1, 2]
            0.000500, # [0, 2, 2]
            0.999000, # [1, 2, 2]
            0.999000, # [2, 2, 2]
        ]

        cpt = dai_util.CPT( [0,1,2], [3,3,3] )
        for i, j in enumerate(dai_util.multi_dim_iter([3,3,3])):
            cpt._table[j] = matrix[i]

        assert cpt._table._li == matrix

    def test_reorder(self):
        matrix_1 = [
            1, 2, 2,
            0, 1, 2,
            0, 0, 1
        ]

        matrix_2 = [
            1, 0, 0,
            2, 1, 0,
            2, 2, 1
        ]

        cpt = dai_util.CPT( [0,1], [3,3] )
        for i, j in enumerate(dai_util.multi_dim_iter([3,3])):
            cpt._table[j] = matrix_1[i]


        print "#", cpt.factors([1,0])
        print "#", cpt.factors([0,1])
        
        out = cpt.factors([1,0])

        print out
        print matrix_2

        assert out == matrix_2
