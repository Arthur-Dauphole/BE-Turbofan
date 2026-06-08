"""
BE 10 - Turbofan : diagramme (T,s) du flux primaire, effet du Mach de vol.
Script AUTONOME : le modele de cycle est inclus, aucun autre fichier requis.
Dependances : numpy, matplotlib   (pip install numpy matplotlib)
"""

import numpy as np
import matplotlib.pyplot as plt


# =====================================================================
#  PARAMETRES (conception + hypotheses)
# =====================================================================
class Parametres:
    def __init__(self):
        # Conditions de vol (croisiere)
        self.M0 = 0.8
        self.P0 = 227e2        # Pa
        self.T0 = 217.0        # K
        # Gaz
        self.r = 287.0         # J/kg/K (avant combustion)
        self.r_e = 291.6       # J/kg/K (apres combustion)
        self.gamma = 1.4
        self.gamma_e = 1.33
        self.Pk = 42.8e6       # J/kg (PCI kerosene)
        # Rendements & pertes
        self.xi_e = 0.98
        self.eta_c_BP = 0.90
        self.eta_c_HP = 0.90
        self.eta_f = 0.92
        self.eta_comb = 0.99
        self.xi_cc = 0.95
        self.eta_m = 0.98
        self.eta_t_HP = 0.89
        self.eta_t_BP = 0.90
        self.xi_tuy = 0.98
        # Parametres de conception
        self.Tt4 = 1600.0
        self.OPR = 40.0
        self.pi_CHP = 22.0
        self.pi_f = 1.45
        self.BPR = 11.0
        # Objectif
        self.F_objectif = 21000.0
        # Contraintes de taille
        self.Mach_fan = 0.6
        self.hub_ratio = 0.3


def cp(gamma, r):
    return gamma * r / (gamma - 1.0)


# =====================================================================
#  MODELE DE CYCLE (uniquement ce qui est utile au diagramme T,s)
# =====================================================================
class Turbofan:
    def __init__(self, p: Parametres):
        self.p = p
        self.cp_air = cp(p.gamma, p.r)
        self.cp_gaz = cp(p.gamma_e, p.r_e)
        self.run()

    def _compression(self, Tt_in, pi, eta):
        g = self.p.gamma
        return Tt_in * pi ** ((g - 1.0) / (g * eta))

    def run(self):
        p = self.p
        cpa, cpg = self.cp_air, self.cp_gaz

        # Conditions amont
        a0 = np.sqrt(p.gamma * p.r * p.T0)
        V0 = p.M0 * a0
        Tt0 = p.T0 * (1 + (p.gamma - 1) / 2 * p.M0**2)
        Pt0 = p.P0 * (1 + (p.gamma - 1) / 2 * p.M0**2) ** (p.gamma / (p.gamma - 1))

        # Entree d'air (0 -> 2)
        Tt2 = Tt0
        Pt2 = p.xi_e * Pt0

        # Fan (2 -> 21)
        Tt21 = self._compression(Tt2, p.pi_f, p.eta_f)
        Pt21 = Pt2 * p.pi_f

        # Compresseur BP (21 -> 25), adapte pour atteindre l'OPR
        pi_BP = p.OPR / (p.pi_f * p.pi_CHP)
        Tt25 = self._compression(Tt21, pi_BP, p.eta_c_BP)
        Pt25 = Pt21 * pi_BP

        # Compresseur HP (25 -> 3)
        Tt3 = self._compression(Tt25, p.pi_CHP, p.eta_c_HP)
        Pt3 = Pt25 * p.pi_CHP

        # Chambre (3 -> 4)
        Tt4 = p.Tt4
        Pt4 = p.xi_cc * Pt3
        f = (cpg * Tt4 - cpa * Tt3) / (p.eta_comb * p.Pk - cpg * Tt4)

        # Turbine HP (4 -> 45) : entraine le CHP
        dW_CHP = cpa * (Tt3 - Tt25)
        dTt_THP = dW_CHP / (p.eta_m * (1 + f) * cpg)
        Tt45 = Tt4 - dTt_THP
        pi_THP = (Tt45 / Tt4) ** (p.gamma_e / (p.eta_t_HP * (p.gamma_e - 1)))
        Pt45 = Pt4 * pi_THP

        # Turbine BP (45 -> 5) : entraine fan + compresseur BP
        dW_fan = (1 + p.BPR) * cpa * (Tt21 - Tt2)
        dW_BP = cpa * (Tt25 - Tt21)
        dTt_TBP = (dW_fan + dW_BP) / (p.eta_m * (1 + f) * cpg)
        Tt5 = Tt45 - dTt_TBP
        pi_TBP = (Tt5 / Tt45) ** (p.gamma_e / (p.eta_t_BP * (p.gamma_e - 1)))
        Pt5 = Pt45 * pi_TBP

        # Tuyere primaire (5 -> 9), adaptee P9 = P0
        Pt9 = p.xi_tuy * Pt5
        Tt9 = Tt5
        T9 = Tt9 * (p.P0 / Pt9) ** ((p.gamma_e - 1) / p.gamma_e)
        V9 = np.sqrt(max(2 * cpg * (Tt9 - T9), 0.0))

        # Rendement global (pour la legende)
        Tt19 = Tt21
        Pt19 = p.xi_tuy * Pt21
        T19 = Tt19 * (p.P0 / Pt19) ** ((p.gamma - 1) / p.gamma)
        V19 = np.sqrt(max(2 * cpa * (Tt19 - T19), 0.0))
        F_par_mp = ((1 + f) * V9 - V0) + p.BPR * (V19 - V0)
        mp = p.F_objectif / F_par_mp
        ms = p.BPR * mp
        mk = f * mp
        W_pr = p.F_objectif * V0
        W_cy = 0.5 * mp * ((1 + f) * V9**2 - V0**2) + 0.5 * ms * (V19**2 - V0**2)
        eta_th = W_cy / (mk * p.Pk)
        eta_pr = W_pr / W_cy
        self.eta_glob = eta_th * eta_pr

        # Stockage des stations utiles au flux primaire
        self.T0, self.P0 = p.T0, p.P0
        self.Tt2, self.Pt2 = Tt2, Pt2
        self.Tt25, self.Pt25 = Tt25, Pt25
        self.Tt3, self.Pt3 = Tt3, Pt3
        self.Tt4, self.Pt4 = Tt4, Pt4
        self.Tt45, self.Pt45 = Tt45, Pt45
        self.Tt5, self.Pt5 = Tt5, Pt5
        self.T9 = T9


# =====================================================================
#  CONSTRUCTION DES POINTS DU CYCLE EN (T, s)
# =====================================================================
def cycle_points(M0):
    p = Parametres()
    p.M0 = M0
    t = Turbofan(p)
    cpa, cpg = cp(p.gamma, p.r), cp(p.gamma_e, p.r_e)

    # Entropie relative s = cp ln(T/Tref) - r ln(P/Pref), ref = amont statique
    Tref, Pref = p.T0, p.P0

    def s_air(T, P):
        return cpa * np.log(T / Tref) - p.r * np.log(P / Pref)

    def s_gaz(T, P):
        return cpg * np.log(T / Tref) - p.r_e * np.log(P / Pref)

    # (nom, T totale, P totale, entropie)
    pts = [
        ("0",  p.T0,   p.P0,   s_air(p.T0,   p.P0)),
        ("2",  t.Tt2,  t.Pt2,  s_air(t.Tt2,  t.Pt2)),
        ("25", t.Tt25, t.Pt25, s_air(t.Tt25, t.Pt25)),
        ("3",  t.Tt3,  t.Pt3,  s_air(t.Tt3,  t.Pt3)),
        ("4",  t.Tt4,  t.Pt4,  s_gaz(t.Tt4,  t.Pt4)),
        ("45", t.Tt45, t.Pt45, s_gaz(t.Tt45, t.Pt45)),
        ("5",  t.Tt5,  t.Pt5,  s_gaz(t.Tt5,  t.Pt5)),
        ("9",  t.T9,   p.P0,   s_gaz(t.T9,   p.P0)),
    ]
    return pts, t


# =====================================================================
#  TRACE
# =====================================================================
def main():
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    colors = {0.8: "#1f77b4", 0.7: "#d62728"}

    for M0 in [0.8, 0.7]:
        pts, t = cycle_points(M0)
        s = [pt[3] for pt in pts]
        T = [pt[1] for pt in pts]
        ax.plot(s, T, "-o", color=colors[M0], lw=2, ms=5,
                label=f"$M_0={M0}$ ($\\eta_g$={t.eta_glob:.3f})")
        # Annoter les stations sur la courbe de reference uniquement
        if M0 == 0.8:
            for name, Tt, Pt, sv in pts:
                ax.annotate(name, (sv, Tt), textcoords="offset points",
                            xytext=(6, 5), fontsize=9, color=colors[M0])

    ax.set_xlabel("Entropie relative  $s$  [J/kg/K]")
    ax.set_ylabel("Température totale  $T_t$  [K]")
    ax.set_title("Diagramme (T,s) du flux primaire — effet du Mach de vol")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("Q5_diagramme_Ts.png", dpi=130)
    print("Figure enregistree : Q5_diagramme_Ts.png")
    # plt.show()   # decommenter pour affichage interactif


if __name__ == "__main__":
    main()
