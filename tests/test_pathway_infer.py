
import unittest
import os
import sys
import subprocess
import shutil
import time

from pathway_tools.infer import pathway_inference, InferenceConfig


class TestCase(unittest.TestCase):    
    def test_infer(self):

        for s in ["sample_1","sample_2","sample_3","sample_4","sample_5","sample_6"]:
            args = InferenceConfig(dogma='active_only.yaml', 
                pathway='toy_pathway.tab', 
                evidence=[('active', 'toy_data.tab')],
                sample=s
            )
            pathway_inference(args)
