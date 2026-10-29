import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import expon, chi2


np.random.seed(100)  # risultati riproducibili

dati = np.array([560, 100, 1230, 200, 550, 740, 1730, 530, 1020, 930, 470, 1330, 1040, 730, 380, 960, 1360, 260, 580, 320,
1000, 870, 1690, 1680, 450, 1220, 240, 1380, 130, 910, 340, 970, 50, 280, 730, 470, 1380, 1120, 420, 1070, 1000,
2200, 1170, 870, 540, 1470, 1950, 20, 290, 340, 290])  # dati veri

N = len(dati)
mu = np.mean(dati)  # media campionaria (stima del parametro di scala beta)
std = np.std(dati, ddof=1)

# Statistiche campionarie
print(f"Media dei tempi: {mu:.4f}")
print(f"Deviazione standard teorica (esponenziale): {mu:.4f}")  
print(f"Deviazione standard dei dati: {std:.4f}")

bin_width = 200
bins = np.arange(dati.min() - 0.5, dati.max() + 1.5, bin_width)  # costruzione dei bin per l'istogramma

# Istogramma con grafica migliorata
counts, a, b = plt.hist(dati, bins=bins, density=True, alpha=0.6, color="skyblue", edgecolor="black", label=f"Dati (bin width = {bin_width}"+r"$\, \mu s$)")
bin_center = 0.5 * (a[:-1] + a[1:])

attese = expon.cdf(bins[1:] - 0.5, scale=mu) - expon.cdf(bins[:-1] - 0.5, scale=mu)
attese *= 1 / bin_width  # teorema della media integrale
chisq = N * np.sum((counts - attese)**2 / attese) * bin_width

# Calcolo del Chi quadro
print(chisq, len(a))
df = len(a) - 1
p_value = chi2.sf(chisq, df)
print(f"p-value = {p_value:.4f}")

# Esponenziale teorica con chi2 e p-value in legenda
k = np.linspace(bins[0], bins[-1], 500)
plt.plot(k, expon.pdf(k, scale=mu), color="crimson", linewidth=2.5, label=f"Esponenziale Teorica")
plt.plot([], [], ' ', label=f" $\chi^2$ = {chisq:.2f}  (df = {df})\np-value = {p_value:.4f}")

# Personalizzazione di layout, titoli ed etichette ingrandite
plt.xlabel(r"$\Delta t \, [\mu s]$", fontsize=12, fontweight="bold")
plt.ylabel("Densità di Probabilità", fontsize=12, fontweight="bold")

plt.xticks(fontsize=11)
plt.yticks(fontsize=11)

plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(fontsize=11, frameon=True, facecolor="white", edgecolor="lightgray", shadow=True)
plt.tight_layout()

plt.show()