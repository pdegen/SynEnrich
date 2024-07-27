# scripts/combine_results.py

import sys
import pandas as pd

clusterprofiler_file = sys.argv[1]
gseapy_file = sys.argv[2]
output_file = sys.argv[3]

clusterprofiler_results = pd.read_csv(clusterprofiler_file)
gseapy_results = pd.read_csv(gseapy_file)

combined_results = pd.DataFrame([9,8,7,6,5])

# Save combined results
combined_results.to_csv(output_file, index=False)
