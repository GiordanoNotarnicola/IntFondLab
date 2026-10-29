import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from scipy.optimize import curve_fit

# Dati sperimentali (altezze picco-picco in mV)
y = np.array([592, 440, 320, 232, 176, 128, 112, 64])

# Intervallo di tempo tra un picco e l'altro (24 ns)
dt = 24
t = np.arange(len(y)) * dt  # 0, 24, 48, 72, 96, 120, 144, 168 ns

# Incertezza costante su ogni misura: t_division / sqrt(3) -> 16 / sqrt(3) mV
err = np.ones(len(y)) * 16 / np.sqrt(3)

# Modello di fit per il decadimento esponenziale smorzato: A(t) = A0 * exp(-t / tau)
def fit_func(t, A0, tau):
    return A0 * np.exp(-t / tau)

# Fit non lineare con curve_fit
popt, pcov = curve_fit(fit_func, t, y, p0=[600, 50], sigma=err, absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))  # Incertezze dei parametri (A0 e tau)

print(f"Parametri fit: A0 = {popt[0]:.2f} ± {perr[0]:.2f}, tau = {popt[1]:.2f} ± {perr[1]:.2f}")

# Calcolo dei residui e del Chi-quadro
res = y - fit_func(t, *popt)
chisq = np.sum((res / err)**2)
df = len(y) - len(popt)  # 8 - 2 = 6 gradi di libertà
p_value = chi2.sf(chisq, df)

print(f"Chi2 = {chisq:.4f}, df = {df}, p-value = {p_value:.4f}")

# Plot con grafica curata
plt.figure(figsize=(8, 6))
t_plot = np.linspace(t[0], t[-1], 500)
plt.plot(t_plot, fit_func(t_plot, *popt), color="crimson", linewidth=2.5, label="Fit Esponenziale Smorzato")
plt.errorbar(t, y, yerr=err, fmt=".", color="midnightblue", ecolor="gray", capsize=4, markersize=6, label="Picchi Sperimentali")

# Voce fittizia in legenda per mostrare i parametri e il test del chi-quadro in modo ordinato
plt.plot([], [], ' ', label=f"$A_0$ = {popt[0]:.1f} $\pm$ {perr[0]:.1f} mV\n "+r"$\tau$"+f" = {popt[1]:.1f} $\pm$ {perr[1]:.1f} ns\n$\chi^2$ = {chisq:.2f} (df = {df})\np-value = {p_value:.3f}")

# Personalizzazione di layout, titoli ed etichette ingrandite
plt.title("Decadimento dell'Ampiezza Picco-Picco", fontsize=14, fontweight="bold", pad=15)
plt.xlabel(r"Tempo $t$ [ns]", fontsize=12, fontweight="bold")
plt.ylabel(r"Ampiezza Picco-Picco [mV]", fontsize=12, fontweight="bold")

plt.xticks(fontsize=11)
plt.yticks(fontsize=11)

plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(fontsize=11, frameon=True, facecolor="white", edgecolor="lightgray", shadow=True, loc="upper right")
plt.tight_layout()

plt.show()