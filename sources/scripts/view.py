#! /usr/bin/env python3
import sys
import pandas as pd


FILE=str(sys.argv[1])
pt = pd.read_table(FILE, delimiter=" ")
print(pt)
