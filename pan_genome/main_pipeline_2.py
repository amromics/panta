import os
import shutil
import re
import json
import gzip
import csv
import logging
from Bio import SeqIO
import pandas as pd
from pan_genome.utils import *

logger = logging.getLogger(__name__)

def add_mcl_cluster(cluster, mcl_file):
    with open(mcl_file, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            genes = line.split('\t')
            cluster.append(genes)
    return 0
