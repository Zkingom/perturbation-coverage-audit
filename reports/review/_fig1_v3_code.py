import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import json
from pathlib import Path

R = Path(r"F:\codex_in_this\claude_science")

def apply_figure_style(sizes=(8, 7, 6)):
    plt.rcParams.update({
        'font.size': sizes[2],
        'axes.titlesize': sizes[1],
        'axes.labelsize': sizes[2],
        'xtick.labelsize': sizes[2],
        'ytick.labelsize': sizes[2],
        'legend.fontsize': sizes[2],
    })

def panel_letter(ax, letter):
    ax.text(-0.18, 1.05, f'({letter})', transform=ax.transAxes, fontsize=8, fontweight='bold', va='top')

apply_figure_style(sizes=(8, 7, 6))

dev = json.loads((R / "results/d1_v0.7dev/confirm_v07_on_dev/confirm_v07_dev.json").read_text(encoding="utf8"))
rp = json.loads((R / "results/conf_v0.7/rpe1/confirm_v07_result.json").read_text(encoding="utf8"))
full = json.loads((R / "results/conf_v0.7/gwps/confirm_v07_result.json").read_text(encoding="utf8"))

RES = {"dev": dev, "rpe1": rp, "gwps": full}

DCOL = {"dev": "#9a9a9a", "rpe1": "#1f5fa8", "gwps": "#c26a1a"}
DLAB = {"dev": "K562 essential (dev.)", "rpe1": "RPE1 (essential)", "gwps": "K562 non-essential (GWPS)"}

tiers = np.arange(1, 6)

meths = [("Mdist-N", "Similarity (norm.)"), ("Mdist-M", "Similarity (Mondrian)"),
         ("MD-GO", "GO distance"), ("MD-STR", "STRING distance"),
         ("Munc-N", "Uncertainty (norm.)"), ("Munc-M", "Uncertainty (Mondrian)")]
yy = np.arange(len(meths))[::-1]


def panel_a(ax):
    for j, ds in enumerate(("dev", "rpe1", "gwps")):
        y = [RES[ds]["tier_tables_rand"]["M0"]["tier_cov"][str(t)] for t in tiers]
        off = (j - 1) * 0.07
        ax.plot(tiers + off, y, marker="o", ms=3.5, lw=1.4, color=DCOL[ds],
                mfc="none" if ds == "dev" else DCOL[ds], label=DLAB[ds])
    ax.axhline(0.9, color="k", lw=0.7, ls=":")
    ax.text(0.85, 0.905, "nominal 0.9", fontsize=6, va="bottom")
    ax.set_xticks(tiers)
    ax.set_xticklabels([f"T{t}" for t in tiers])
    ax.set_xlabel("Effect-strength quintile")
    ax.set_ylabel("Coverage")
    ax.set_ylim(0.4, 1.05)
    ax.set_title("Strongest 20% under-covered", loc="left")
    ax.legend(frameon=False, loc="lower left", fontsize=6)


def panel_b(ax):
    for ds in ("rpe1", "gwps"):
        for (m, _), y0 in zip(meths, yy):
            e_ = RES[ds]["descriptive_rand"][m]["E1_vs_M0"]
            off = 0.15 if ds == "rpe1" else -0.15
            ax.errorbar(100 * e_["diff"], y0 + off,
                        xerr=[[100 * (e_["diff"] - e_["lo"])], [100 * (e_["hi"] - e_["diff"])]],
                        fmt="o", ms=3.5, color=DCOL[ds], elinewidth=1, capsize=0)
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(yy)
    ax.set_yticklabels([l for _, l in meths])
    ax.set_xlabel("T5 coverage change (pp)")
    ax.set_title("Difficulty proxies: no T5 gain", loc="left")
    ax.margins(x=0.06)


def panel_c(ax):
    for ds in ("dev", "rpe1", "gwps"):
        P = RES[ds]["P"]
        e_ = RES[ds]["step1"]["E1"]
        ax.errorbar(P["rho"], 100 * e_["diff"],
                    xerr=[[P["rho"] - P["rho_ci"][0]], [P["rho_ci"][1] - P["rho"]]],
                    yerr=[[100 * (e_["diff"] - e_["lo"])], [100 * (e_["hi"] - e_["diff"])]],
                    fmt="o", ms=5, color=DCOL[ds], mfc="none" if ds == "dev" else DCOL[ds],
                    elinewidth=1, capsize=0, zorder=3)
    ax.axhline(5, color="k", lw=0.6, ls=":")
    ax.text(0.02, 5.4, "pre-registered 5 pp margin", fontsize=6)
    ax.axhline(0, color="k", lw=0.6)
    ax.text(full["P"]["rho"] + 0.02, 100 * full["step1"]["E1"]["diff"] - 2.6,
            "GWPS", fontsize=6, color=DCOL["gwps"])
    ax.text(rp["P"]["rho"] - 0.03, 100 * rp["step1"]["E1"]["diff"] + 2.4,
            "RPE1", fontsize=6, color=DCOL["rpe1"], ha="right")
    ax.text(dev["P"]["rho"] + 0.02, 100 * dev["step1"]["E1"]["diff"] + 0.8,
            "dev.", fontsize=6, color=DCOL["dev"])
    ax.text(0.02, 16.3, "oracle: +48 / 52 / 47 pp", fontsize=6, va="top")
    ax.set_xlabel("Proxy–effect-strength ρ")
    ax.set_ylabel("T5 gain of proxy method (pp)")
    ax.set_xlim(0, 0.6)
    ax.set_ylim(-3, 17)
    ax.set_title("Repair tracks predictability", loc="left")


plt.close("all")

fig = plt.figure(figsize=(7.4, 2.8))
axa = fig.add_axes([0.065, 0.17, 0.25, 0.70])
axb = fig.add_axes([0.505, 0.17, 0.19, 0.70])
axc = fig.add_axes([0.785, 0.17, 0.205, 0.70])

panel_a(axa)
panel_b(axb)
panel_c(axc)

axa.legend(frameon=False, loc="center left", bbox_to_anchor=(0.0, 0.42), fontsize=6)

for a_, l in zip((axa, axb, axc), "abc"):
    panel_letter(a_, l)

fig.savefig("fig1_main_v3.png", dpi=300, bbox_inches="tight")
fig.savefig("fig1_main_v3.pdf", bbox_inches="tight")