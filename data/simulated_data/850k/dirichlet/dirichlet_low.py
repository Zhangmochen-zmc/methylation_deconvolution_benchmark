import numpy as np
import pandas as pd

num_samples = 100
num_categories = 6

alpha_vector = [2, 2, 2, 2, 0.5, 2]
data = np.random.dirichlet(alpha_vector, size=num_samples)

data = np.round(data, 4)

for i in range(num_samples):
    current_sum = np.sum(data[i])
    diff = 1.0 - current_sum

    if diff != 0:
        random_index = np.random.randint(0, num_categories)
        data[i, random_index] = round(data[i, random_index] + diff, 4)

cell_types = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Monocytes",
    "Neutrophils",
    "Natural Killer cells"
]
df = pd.DataFrame(data, columns=cell_types)
df.to_csv("850klow_neutrophil.csv", index=False, float_format='%.4f')
