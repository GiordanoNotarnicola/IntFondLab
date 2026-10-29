import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson, chi2
from scipy.optimize import curve_fit
from scipy.special import factorial


np.random.seed(100)  # risultati riproducibili

# Simulazione

mu_teorica = 0  # inutile quando sostituisco i dati
dati = np.array([10, 16, 12, 15, 15, 18, 10, 17, 24, 19, 13, 12, 19, 17, 12, 17, 20, 20, 20, 26, 25, 16, 19, 19, 24, 20, 15, 13, 19, 20, 18, 12, 11, 16, 19, 26, 15, 16, 15, 18, 22, 10, 16, 24, 13, 18, 17,
                14, 20, 20, 22, 18, 20, 17, 12, 16, 19, 17, 21, 11, 19, 15, 17, 18, 21, 18, 16, 18, 14, 27, 17, 22, 16, 12, 11, 14, 24, 15, 21, 14, 15, 11, 11, 18, 19, 17, 16, 22, 17, 8, 13, 21, 16,
                15, 17, 13, 20, 15, 14, 16, 16, 17, 18, 16, 13, 28, 23, 17, 18, 10, 21, 19, 13, 16, 13, 24, 24, 17, 23, 16, 16, 21, 16, 13, 18, 20, 18, 20, 13, 25, 16, 18, 13, 20, 16, 21, 18, 14, 16, 
                14, 24, 13, 19, 17, 22, 18, 18, 14, 13, 10, 17, 22, 13, 13, 18, 16, 20, 22, 20, 18, 10, 23, 19, 18, 21, 15, 21, 19, 21, 19, 17, 17, 22, 17, 15, 16, 9, 18, 21, 18, 21, 12, 15, 24, 18, 
                8, 26, 15, 15, 15, 14, 15, 17, 21, 17, 13, 15, 21, 10, 15, 14, 20, 16, 15, 21, 14, 19, 12, 18, 14, 18, 20, 13, 23, 17, 16, 25, 14, 25, 16, 10, 21, 20, 24, 10, 18, 11, 20, 10, 23, 21, 
                15, 20, 19, 30, 12, 13, 18, 19, 15, 22, 14, 11, 16, 21, 15, 15, 16, 24, 18, 16, 16, 24, 19, 17, 14, 20, 13, 14, 24, 17, 15, 22, 19, 13, 13, 18, 13, 20, 13, 15, 23, 11, 16, 9, 14, 19, 
                23, 21, 18, 17, 17, 18, 16, 15, 21, 19, 22, 16, 14, 15, 16, 15, 11, 17, 26, 17, 20, 17, 13, 17, 23, 21, 15, 17, 12, 19, 29])
N = len(dati)
print(N)
mu = np.mean(dati)
std = np.std(dati, ddof=1)


# Statistiche campionarie

print(f"Media dei dati: {mu:.2f}")
print(f"Deviazione standard poissoniana: {np.sqrt(mu):.2f}")
print(f"Deviazione standard dei dati: {std:.2f}")

# Istogramma delle frequenze relative
bin_width = 1
bins = np.arange(dati.min() - 0.5, dati.max() + 1.5, bin_width)  # costruzione dei bin per l'istogramma
counts, a, b = plt.hist(dati, bins=bins, density=True, alpha=0.6, color="skyblue", edgecolor="black", label=f"Dati  (bin width = {bin_width}")

bin_center = 0.5 * (a[:-1] + a[1:])
err = np.sqrt(counts * (1 - counts) / N)


def fit(k, a, mu):
    return a * poisson.pmf(np.round(k), mu)

popt, pcov = curve_fit(fit, bin_center, counts, p0=[1, mu], sigma=err)
print(popt, np.sqrt(pcov[1, 1]), np.sqrt(popt[1]))
res = counts - fit(bin_center, *popt)
chisq = np.sum((res / err)**2)
df=len(bin_center)-len(popt)
p_value = chi2.sf(chisq, df)
print(f"p-value = {p_value:.4f}")

# Plot con grafica migliorata e Chi-quadro in legenda
plt.plot(bin_center, fit(bin_center, *popt), color="crimson", linewidth=2.5, label=f"Fit di Poisson ($\mu$ = {popt[1]:.2f} $\pm$ {pcov[1,1]:.2f})")
plt.errorbar(bin_center, counts, yerr=err, fmt=".", color="midnightblue", ecolor="gray", capsize=3, label=f"Dati con Incertezza")
plt.plot([], [], ' ', label=f"  $\hat\mu$ = {mu:.2f} $\pm$ {np.sqrt(mu):.2f} \n $\chi^2$ = {chisq:.2f}  (df = {df})\np-value = {p_value:.4f} ")
# Personalizzazione di layout, titoli ed etichette ingrandite
plt.xlabel("Conteggio", fontsize=12, fontweight="bold")
plt.ylabel("Densità di Probabilità", fontsize=12, fontweight="bold")

plt.xticks(fontsize=11)
plt.yticks(fontsize=11)

plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(fontsize=11, frameon=True, facecolor="white", edgecolor="lightgray", shadow=True)
plt.tight_layout()

plt.show()