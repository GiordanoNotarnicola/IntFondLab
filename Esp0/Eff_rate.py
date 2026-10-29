import numpy as np
import matplotlib.pyplot as plt

# Dati sperimentali
v2 = np.array([1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900])
pmt1 = np.array([29048, 29255, 29445, 29373, 29567, 29546, 29850, 30477, 29926])
pmt2 = np.array([5, 39, 210, 577, 1504, 6528, 9052, 10169, 11229])
pmt3 = np.array([73023, 66213, 66153, 67920, 66762, 70906, 66639, 66554, 66374])
doppie = np.array([274, 270, 268, 251, 273, 258, 278, 250, 253])
triple = np.array([1, 17, 73, 133, 175, 171, 189, 177, 172])

# Calcolo dell'efficienza e della sua incertezza binomiale
efficienza = triple / doppie
eff_err = np.sqrt(efficienza * (1 - efficienza) / doppie)

# Stampa dei valori calcolati per controllo
for v, eff, err in zip(v2, efficienza, eff_err):
    print(f"V2 = {v} V | Efficienza = {eff:.4f} ± {err:.4f}")

# Plot con grafica curata
plt.figure(figsize=(8, 6))
plt.errorbar(v2, efficienza, yerr=eff_err, fmt=".", color="crimson", ecolor="gray", 
             elinewidth=1.5, capsize=4, markersize=6, label="Efficienza Sperimentale")

# Personalizzazione di layout, titoli ed etichette ingrandite
plt.title("Curva di Efficienza in Funzione del Voltaggio $V_2$", fontsize=14, fontweight="bold", pad=15)
plt.xlabel(r"Voltaggio $V_2$ [V]", fontsize=12, fontweight="bold")
plt.ylabel(r"Efficienza ($\text{Triple} / \text{Doppie}$)", fontsize=12, fontweight="bold")

plt.xticks(fontsize=11)
plt.yticks(fontsize=11)

plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(fontsize=11, frameon=True, facecolor="white", edgecolor="lightgray", shadow=True, loc="lower right")
plt.tight_layout()

plt.show()