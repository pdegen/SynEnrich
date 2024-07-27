# scripts/run_gseapy.py

import sys
import pandas as pd
import gseapy as gp

input_file = sys.argv[1]
output_file = sys.argv[2]

data = pd.read_csv(input_file)

data = pd.DataFrame([1,2,3,4,5])

# Save results
data.to_csv(output_file, index=False)
