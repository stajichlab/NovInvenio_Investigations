#!/usr/bin/env python3.12
"""AMI for the #132 priority-4 Leiden sweep, from the saved partitions.

Needs scikit-learn (user-site install for /usr/bin/python3.12, 2026-09-28).
Per dataset and resolution: mean and min pairwise AMI between seeds, and AMI
of seed 0 against seed 0 at the next lower resolution.

Usage: p4_ami.py <dataset>.labels.npz ... > p4_ami.tsv
"""
import itertools
import os
import sys

import numpy as np
from sklearn.metrics import adjusted_mutual_info_score as ami

print('dataset\tresolution\tseed_AMI_mean\tseed_AMI_min\tAMI_vs_prev_res')
for path in sys.argv[1:]:
    d = np.load(path)
    name = os.path.basename(path).replace('.labels.npz', '')
    res, lab = d['resolutions'], d['labels']          # lab[res, seed, family]
    for i, r in enumerate(res):
        pairs = [ami(a, b) for a, b in itertools.combinations(lab[i], 2)]
        prev = f'{ami(lab[i - 1][0], lab[i][0]):.3f}' if i else ''
        print(f'{name}\t{r:g}\t{np.mean(pairs):.3f}\t{min(pairs):.3f}\t{prev}', flush=True)
