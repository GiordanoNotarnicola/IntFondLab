import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson, chi2


np.random.seed(100)  # risultati riproducibili

# Simulazione

mu_teorica=0 #inutile quando sostituisco i dati
dati = np.array([10, 16, 12, 15, 15, 18, 10, 17, 24, 19, 13, 12, 19, 17, 12, 17, 20, 20, 20, 26, 25, 16, 19, 19, 24, 20, 15, 13, 19, 20, 18, 12, 11, 16, 19, 26, 15, 16, 15, 18, 22, 10, 16, 24, 13, 18, 17,
                14, 20, 20, 22, 18, 20, 17, 12, 16, 19, 17, 21, 11, 19, 15, 17, 18, 21, 18, 16, 18, 14, 31, 27, 17, 22, 16, 12, 11, 14, 24, 15, 21, 14, 15, 11, 11, 18, 19, 17, 16, 22, 17, 8, 13, 21, 16,
                15, 17, 13, 20, 15, 14, 16, 16, 17, 18, 16, 13, 28, 23, 17, 18, 10, 21, 19, 13, 16, 13, 24, 24, 17, 23, 16, 16, 21, 16, 13, 18, 20, 18, 20, 13, 25, 16, 18, 13, 20, 16, 21, 18, 14, 16, 
                14, 24, 13, 19, 17, 22, 18, 18, 14, 13, 10, 17, 22, 13, 13, 18, 16, 20, 22, 20, 18, 10, 23, 19, 18, 21, 15, 21, 19, 21, 19, 17, 17, 22, 17, 15, 16, 9, 18, 21, 18, 21, 12, 15, 24, 18, 
                8, 26, 15, 15, 15, 14, 15, 17, 21, 17, 13, 15, 21, 10, 15, 14, 20, 16, 15, 21, 14, 19, 12, 18, 14, 18, 20, 13, 23, 17, 16, 25, 14, 25, 16, 10, 21, 20, 24, 10, 18, 11, 20, 10, 23, 21, 
                15, 20, 19, 30, 12, 13, 18, 19, 15, 22, 14, 11, 16, 21, 15, 15, 16, 24, 18, 16, 16, 24, 19, 17, 14, 20, 13, 14, 24, 17, 15, 22, 19, 13, 13, 18, 13, 20, 13, 15, 23, 11, 16, 9, 14, 19, 
                23, 21, 18, 17, 17, 18, 16, 15, 21, 19, 22, 16, 14, 15, 16, 15, 11, 17, 26, 17, 20, 17, 13, 17, 23, 21, 15, 17, 12, 19])
N=len(dati)
print(N)
mu= np.mean(dati)
std = np.std(dati, ddof=1)


# Statistiche campionarie

print(f"Media dei dati: {mu:.2f}")
print(f"Deviazione standard poissoniana: {np.sqrt(mu):.2f}")
print(f"Deviazione standard dei dati: {std:.2f}")

# Istogramma delle frequenze relative
bin_width=1
bins = np.arange(dati.min() - 0.5, dati.max()+1.5 , bin_width)#costruzione dei bin per l'istogramma
k = np.arange(int(bins[0]), int(bins[-1]), bin_width)
counts,a,b=plt.hist(dati, bins=bins, density=True, label="Dati")

# Poisson teorica
#plt.plot(k, poisson.pmf(k, mu), "r-", label="Poisson teorica")


attese = poisson.cdf(bins[1:] - 0.5, mu) - poisson.cdf(bins[:-1] - 0.5, mu)
attese *= 1 / bin_width#sto utilizzando il teorema della media integrale
chisq = N * np.sum((counts - attese)**2 / attese) * bin_width
print(chisq, len(k))
df=len(k)-1
p_value = chi2.sf(chisq, df)
print(f"p-value = {p_value:.4f}")

plt.xlabel("Conteggio")
plt.ylabel("Probabilità")
plt.legend()
plt.show()