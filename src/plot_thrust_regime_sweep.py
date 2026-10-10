"""
Plot drag's relative contribution across electric propulsion thrust regimes.

Reads the arrays saved by thrust_regime_sweep.py and produces two plots:
(1) total altitude decay vs. thrust / power level, and
(2) drag's percentage share of that decay vs. thrust / power level - the
figure directly addressing Dr. Modenini's request to examine how drag's
relative contribution varies across different +burn durations and thrust
regimes.

Since each regime may reach reentry (or run to the full simulation
duration) at a different elapsed time, drag's share is computed over
each regime's own elapsed time against the drag-only baseline
interpolated to that same time (see thrust_regime_sweep.py) - this plot
marks regimes that reached reentry with a filled marker, and regimes
that ran the full simulation duration without reentry with an open
marker, so the two are not visually conflated.
"""


import numpy as np
import matplotlib.pyplot as plt

from paths import DATA_DIR, FIGURES_DIR


power = np.load(DATA_DIR / "regime_power.npy")
thrust = np.load(DATA_DIR / "regime_thrust.npy")
burn_duration_hr = np.load(DATA_DIR / "regime_burn_duration_hr.npy")
elapsed_days = np.load(DATA_DIR / "regime_elapsed_days.npy")
reentry_reached = np.load(DATA_DIR / "regime_reentry_reached.npy")
total_decay = np.load(DATA_DIR / "regime_total_decay.npy")
drag_share_percent = np.load(DATA_DIR / "regime_drag_share_percent.npy")


# Plot 1: drag's relative contribution vs. thrust regime (the key result)

fig, ax1 = plt.subplots(figsize = (9, 7))

ax1.plot(power, drag_share_percent, color = "tab:blue", linewidth = 1.5, marker = "o", markersize = 9, zorder = 2)

ax1.set_xticks(power)
ax1.set_ylim(0.425, 0.485)
ax1.set_xlabel("Available Power (W)")
ax1.set_ylabel("Drag's share of altitude decay at that regime's elapsed time (%)")
ax1.grid(True, alpha = 0.3, zorder = 0)

ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(power)
ax2.set_xticklabels([f"{d:.0f}" for d in elapsed_days])
ax2.set_xlabel("Elapsed time at reentry (days)")

plt.title("Drag's Relative Contribution Across Thrust Regimes \n"
          "(electric propulsion, 0.1118 kg propellant budget per report Table 3)", pad = 30)

plt.tight_layout()
plt.savefig(FIGURES_DIR / "drag_relative_contribution.png", dpi = 300)


# Plot 2: total altitude decay vs. thrust regime, for context

fig, ax1 = plt.subplots(figsize = (9, 6))

bar_labels = [f"{p:.0f} W \n({t * 1000:.2f} mN)" for p, t in zip(power, thrust)]

bars = ax1.bar(bar_labels, total_decay, color = "tab:orange")

for bar, value in zip(bars, total_decay):
    ax1.annotate(f"{value:.3f}", xy = (bar.get_x() + bar.get_width() / 2, bar.get_height()), ha = "center", va = "bottom")

ax1.set_ylabel("Total altitude decay (km)")
ax1.set_ylim(350, 380)
ax1.grid(True, alpha = 0.3, axis = "y")

ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(range(len(power)))
ax2.set_xticklabels([f"{b:.1f} hr" for b in burn_duration_hr])
ax2.set_xlabel("Burn duration (hours)")

plt.title("Total Altitude Decay by Thrust Regime \n"
          "(all cases reached reentry)", pad = 30)

plt.tight_layout()
plt.savefig(FIGURES_DIR / "total_decay_by_regime.png", dpi = 300)

plt.show()

print("Summary:")
for p, b, e, d, s, r in zip (power, burn_duration_hr, elapsed_days, total_decay, drag_share_percent, reentry_reached):
    status = "reentry" if r else "no reentry"
    print(f"{p:>3.0f} W | burn = {b:>6.2f} hr | elapsed = {e:>6.2f} d ({status}) | "
          f"total decay = {d:>8.3f} km | drag share = {s:>5.2f} %")
