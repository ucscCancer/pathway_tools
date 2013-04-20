
import unittest
import os
import sys
import subprocess
import shutil
import time

from pathway_tools.infer import pathway_inference, InferenceConfig


class TestCase(unittest.TestCase):    
    def test_infer_1(self):

        for s in ["sample_1","sample_2","sample_3","sample_4","sample_5","sample_6"]:
            args = InferenceConfig(dogma='active_only.yaml', 
                pathway='toy_pathway_1.tab', 
                evidence=[('active', 'toy_data_1.tab')],
                sample=[s]
            )
            pathway_inference(args)


    def test_infer_2(self):

        for s in ["sample_1","sample_2","sample_3","sample_4"]:
            args = InferenceConfig(dogma='active_only.yaml', 
                pathway='toy_pathway_2.tab', 
                evidence=[('active', 'toy_data_2.tab')],
                sample=[s]
            )
            pathway_inference(args)
