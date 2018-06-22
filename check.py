import pandas as pd
import numpy as np

py = pd.read_csv("ratio_test.csv")
r = pd.read_csv("R_analysis.csv")

r.sort_values(by=["Guide.Sequence"], ascending=True, inplace=True)
py.sort_values(by=["sgRNA"], ascending=True, inplace=True)
r.reset_index(inplace=True)
py.reset_index(inplace=True)


compare = pd.DataFrame()
print(np.array_equal(py['sgRNA'], r['Guide.Sequence']))
compare['Gene'] = py['Gene']
compare['Python.Top5.Geom.Mean'] = py['Top5.Geom.Mean']
compare['R.Top5.Geom.Mean'] = r['Top5.Geom.Mean']
compare['Top5.Geom.Mean.Ratio'] = py['Top5.Geom.Mean'] / r['Top5.Geom.Mean']
compare['Python.Bot5.Geom.Mean'] = py['Bot5.Geom.Mean']
compare['R.Bot5.Geom.Mean'] = r['Bot5.Geom.Mean']
compare['Bot5.Geom.Mean.Ratio'] = py['Bot5.Geom.Mean'] / r['Bot5.Geom.Mean']

print(compare)
