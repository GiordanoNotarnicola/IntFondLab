import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import expon, chi2


np.random.seed(100)  # risultati riproducibili

dati = np.array([560, 100, 1230, 200, 550, 740, 1730, 530, 1020, 930, 470, 1330, 1040, 730, 380, 960, 1360, 260, 580, 320,
1000, 870, 1690, 1680, 450, 1220, 240, 1380, 130, 910, 340, 970, 50, 280, 730, 470, 1380, 1120, 420, 1070, 1000,
2200, 1170, 870, 540, 1470, 1950, 20, 290, 340, 290])  # qua inserisco i dati veri
N=len(dati)
mu = np.mean(dati)  # media campionaria (stima del parametro di scala beta)
std = np.std(dati, ddof=1)

# Statistiche campionarie
print(f"Media dei tempi: {mu:.4f}")
print(f"Deviazione standard teorica (esponenziale): {mu:.4f}")  # Per l'esponenziale dev.std = media = beta
print(f"Deviazione standard dei dati: {std:.4f}")

bin_width=200
bins = np.arange(dati.min() - 0.5, dati.max()+1.5 , bin_width)#costruzione dei bin per l'istogramma
k = np.arange(int(bins[0]), int(bins[-1]), bin_width)
counts, a, b = plt.hist(dati, bins=bins, density=True, label="Dati")


# Esponenziale teorica (usiamo scale=mu per stimarla dai dati)
plt.plot(k, expon.pdf(k, scale=mu), "r-", label="Esponenziale teorica")

attese = expon.cdf(bins[1:] - 0.5, mu) - expon.cdf(bins[:-1] - 0.5, mu)
attese *= 1 / bin_width#sto utilizzando il teorema della media integrale
chisq = N * np.sum((counts - attese)**2 / attese) * bin_width
# Calcolo del Chi quadro confrontando la densità osservata con quella teorica nei centri dei bin
print(chisq, len(k))
print(attese)
df=len(k)-1
p_value = chi2.sf(chisq, df)
print(f"p-value = {p_value:.4f}")

plt.xlabel("Delta_t")
plt.ylabel("Densità di probabilità")
plt.legend()
plt.show()