
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson

np.random.seed(38)  # risultati riproducibili

# Simulazione
N = 1000
mu_teorica=300 #inutile quando sostituisco i dati
dati = np.random.poisson(mu_teorica, N)# qua inserisco i dati veri
mu= np.mean(dati)
std = np.std(dati, ddof=1)


# Statistiche campionarie

print(f"Media dei dati: {mu:.2f}")
print(f"Deviazione standard poissoniana: {np.sqrt(mu):.2f}")
print(f"Deviazione standard dei dati: {std:.2f}")

# Istogramma delle frequenze relative
bins = np.arange(dati.min() - 0.5, dati.max() + 1.5, 1)
k = np.arange(dati.min(), dati.max() + 1)

counts,a,b=plt.hist(dati, bins=bins, density=True, alpha=0.6,
         label="Dati simulati")
print(len(bins), len(k), len(dati), len(counts))
# Poisson teorica

plt.plot(k, poisson.pmf(k, mu), "r-",
         label="Poisson teorica")

chi2=np.sum((counts/N- poisson.pmf(k, mu))**2/poisson.pmf(k, mu_teorica))
print(chi2)

plt.xlabel("Conteggio")
plt.ylabel("Probabilità")
plt.legend()
plt.show()