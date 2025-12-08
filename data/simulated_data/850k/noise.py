import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def add_gaussian_noise_via_m_value(beta_matrix, noise_std=0.1, epsilon=1e-6):
    # Clipping
    beta_clipped = np.clip(beta_matrix, epsilon, 1 - epsilon)

    # Logit transform base 2
    #M = log2( Beta / (1 - Beta) )
    m_values = np.log2(beta_clipped / (1 - beta_clipped))

    #Gaussian noise
    noise = np.random.normal(loc=0, scale=noise_std, size=m_values.shape)
    m_noisy = m_values + noise

    # Inverse Logit
    # Beta = 2^M / (1 + 2^M)
    beta_noisy = (2**m_noisy) / (1 + 2**m_noisy)

    return beta_noisy

file_path = '../mix_data/simulated_matrix_850krandom.csv'
df = pd.read_csv(file_path, index_col=0)
print(f" {df.shape} (Row: CpG site, Column: sample)")


#0.2，0.5，1，2
target_noise_std = 6

noisy_data_array = add_gaussian_noise_via_m_value(df, noise_std=target_noise_std)

df_noisy = pd.DataFrame(noisy_data_array, index=df.index, columns=df.columns)


output_filename = '../mix_data/850krandom_noisy6.csv'
df_noisy.to_csv(output_filename)
print(f"Save to: {output_filename}")

#Visualization
original_flat = df.values.flatten()
noisy_flat = df_noisy.values.flatten()

if len(original_flat) > 5000:
    indices = np.random.choice(len(original_flat), 5000, replace=False)
    sample_orig = original_flat[indices]
    sample_noisy = noisy_flat[indices]
else:
    sample_orig = original_flat
    sample_noisy = noisy_noisy

plt.figure(figsize=(8, 8))
plt.scatter(sample_orig, sample_noisy, alpha=0.3, s=5)
plt.plot([0, 1], [0, 1], 'r--', label='y=x')
plt.title(f"Original vs Noisy Data (Noise STD={target_noise_std})")
plt.xlabel("Original Beta Values")
plt.ylabel("Noisy Beta Values")
plt.legend()
plt.show()
