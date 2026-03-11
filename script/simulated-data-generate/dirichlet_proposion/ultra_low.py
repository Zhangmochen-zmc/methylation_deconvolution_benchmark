import numpy as np
import pandas as pd

num_samples = 100
num_categories = 6
alpha_vector = [1, 1, 1, 1, 1, 1]

cell_types = [
    "B cells", 
    "CD4+ T cells", 
    "CD8+ T cells", 
    "Monocytes", 
    "Neutrophils",   
    "Natural Killer cells"
]
rare_idx = cell_types.index("Natural Killer cells")
rare_min = 0.005   
rare_max = 0.01   
data = np.random.dirichlet(alpha_vector, size=num_samples)
for i in range(num_samples):
    rare_prop = np.random.uniform(rare_min, rare_max)
    data[i, rare_idx] = rare_prop
  
    data[i] = data[i] / data[i].sum()

data = np.round(data, 4)

for i in range(num_samples):
    current_sum = np.sum(data[i])
    diff = 1.0 - current_sum
    
    if diff != 0:
        random_index = np.random.randint(0, num_categories)
        data[i, random_index] = round(data[i, random_index] + diff, 4)

df = pd.DataFrame(data, columns=cell_types)

df.to_csv("wgbslow_nk.csv", index=False, float_format='%.4f')
