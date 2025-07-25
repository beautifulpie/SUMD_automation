import argparse
import os
import sys
import numpy as np
import hashlib
import time
import json
import shutil
from datetime import datetime
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Atom import Atom
import logging

class NumpyEncoder(json.JSONEncoder):
    """Numpy 타입을 JSON 직렬화 가능하게 변환"""
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NumpyEncoder, self).default(obj)
