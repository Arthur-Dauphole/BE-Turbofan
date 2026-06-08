"""
BE 10 - Modèle paramétrique de cycle intensif d'un turbofan double-corps double-flux (LEAP-1A)
ISAE-SUPAERO, 2e année - Mécanique et thermodynamique des fluides
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Propriétés des gaz
# =============================================================================
gamma = 1.4          # ratio des chaleurs spécifiques (air)
r = 287.0            # J/(kg·K), constante gaz parfait air
cp = gamma * r / (gamma - 1)  # J/(kg·K)

gamma_star = 1.33    # ratio des chaleurs spécifiques (gaz brûlés)
r_star = 291.6       # J/(kg·K)
cp_star = gamma_star * r_star / (gamma_star - 1)

Pk = 42_800_000      # J/kg, pouvoir calorifique inférieur du kérosène


def turbofan_cycle(M0=0.8, T0=217.0, P0=22700.0, F_target=21000.0,
                   BPR=11.0, pi_f=1.45, OPR=40.0, pi_CHP=22.0, Tt4=1600.0,
                   # Rendements et pertes de charge
                   xi_e=0.98, eta_c_BP=0.90, eta_c_HP=0.90, eta_f=0.92,
                   eta_comb=0.99, xi_cc=0.95, eta_m=0.98,
                   eta_t_HP=0.89, eta_t_BP=0.90, xi_tuy=0.98,
                   # Sizing
                   M_fan=0.6, hub_tip_ratio=0.3):
    """
    Calcule le cycle thermodynamique complet d'un turbofan double-corps.
    Retourne un dictionnaire avec tous les résultats.
    """
    lam = BPR  # λ = ṁs/ṁp

    # --- Station 0 : conditions de vol ---
    a0 = np.sqrt(gamma * r * T0)
    V0 = M0 * a0
    Tt0 = T0 * (1 + (gamma - 1) / 2 * M0**2)
    Pt0 = P0 * (Tt0 / T0) ** (gamma / (gamma - 1))

    # --- Station 2 : sortie entrée d'air ---
    Tt2 = Tt0
    Pt2 = xi_e * Pt0

    # --- Taux de compression du compresseur BP ---
    pi_cBP = OPR / (pi_f * pi_CHP)

    # --- Fan (station 2 → sortie fan) ---
    Tt_fan = Tt2 * pi_f ** ((gamma - 1) / (gamma * eta_f))
    Pt_fan = pi_f * Pt2

    # --- Flux secondaire (station 21) : après le fan ---
    Tt21 = Tt_fan
    Pt21 = Pt_fan

    # --- Compresseur BP (station fan → 2.5) : flux primaire seulement ---
    Tt25 = Tt_fan * pi_cBP ** ((gamma - 1) / (gamma * eta_c_BP))
    Pt25 = pi_cBP * Pt_fan

    # --- Compresseur HP (station 2.5 → 3) ---
    Tt3 = Tt25 * pi_CHP ** ((gamma - 1) / (gamma * eta_c_HP))
    Pt3 = pi_CHP * Pt25

    # --- Chambre de combustion (station 3 → 4) ---
    # Bilan : ṁp·cp·Tt3 + ṁk·ηcomb·Pk = (ṁp + ṁk)·cp*·Tt4
    # f = ṁk/ṁp
    f = (cp_star * Tt4 - cp * Tt3) / (eta_comb * Pk - cp_star * Tt4)
    Pt4 = xi_cc * Pt3
    # Tt4 donné

    # --- Turbine HP (station 4 → 4.5) ---
    # Bilan arbre HP : ηm · (1+f)·cp*·(Tt4 - Tt45) = cp·(Tt3 - Tt25)
    delta_Tt_HP = cp * (Tt3 - Tt25) / (eta_m * (1 + f) * cp_star)
    Tt45 = Tt4 - delta_Tt_HP
    # Pression : relation polytropique (détente)
    Pt45 = Pt4 * (Tt45 / Tt4) ** (gamma_star / (eta_t_HP * (gamma_star - 1)))

    # --- Turbine BP (station 4.5 → 5) ---
    # Bilan arbre BP : ηm · (1+f)·cp*·(Tt45 - Tt5) = cp·[(1+λ)(Tt_fan - Tt2) + (Tt25 - Tt_fan)]
    W_fan = (1 + lam) * cp * (Tt_fan - Tt2)       # travail fan (primaire + secondaire)
    W_cBP = cp * (Tt25 - Tt_fan)                    # travail compresseur BP (primaire seul)
    W_LP_total = W_fan + W_cBP                       # travail total arbre BP par kg de flux primaire

    delta_Tt_BP = W_LP_total / (eta_m * (1 + f) * cp_star)
    Tt5 = Tt45 - delta_Tt_BP
    Pt5 = Pt45 * (Tt5 / Tt45) ** (gamma_star / (eta_t_BP * (gamma_star - 1)))

    # --- Tuyère primaire (station 9) : adaptée → P9 = P0 ---
    Pt9 = xi_tuy * Pt5
    Tt9 = Tt5
    ratio9 = 1 - (P0 / Pt9) ** ((gamma_star - 1) / gamma_star)
    if ratio9 <= 0:
        return None  # cycle non physique
    V9 = np.sqrt(2 * cp_star * Tt9 * ratio9)

    # --- Tuyère secondaire (station 19) : adaptée → P19 = P0 ---
    Pt19 = xi_tuy * Pt21
    Tt19 = Tt21
    V19 = np.sqrt(2 * cp * Tt19 * (1 - (P0 / Pt19) ** ((gamma - 1) / gamma)))

    # --- Poussée et débit ---
    # F = ṁp·[(1+f)·V9 - V0] + λ·ṁp·(V19 - V0)
    specific_thrust_primary = (1 + f) * V9 - V0 + lam * (V19 - V0)  # F/ṁp
    m_dot_p = F_target / specific_thrust_primary  # kg/s flux primaire
    m_dot_s = lam * m_dot_p                        # kg/s flux secondaire
    m_dot_total = m_dot_p + m_dot_s                # kg/s total air
    m_dot_k = f * m_dot_p                          # kg/s carburant

    # --- Dimensionnement (rayon fan) ---
    # Fan face : Mach M_fan, conditions Pt2, Tt2
    T2 = Tt2 / (1 + (gamma - 1) / 2 * M_fan**2)
    P2 = Pt2 * (T2 / Tt2) ** (gamma / (gamma - 1))
    rho2 = P2 / (r * T2)
    V2 = M_fan * np.sqrt(gamma * r * T2)
    A_fan = m_dot_total / (rho2 * V2)
    r_max = np.sqrt(A_fan / (np.pi * (1 - hub_tip_ratio**2)))
    d_max = 2 * r_max

    # --- Performances ---
    W_pr = F_target * V0                                                         # puissance propulsive
    W_cy = 0.5 * m_dot_p * ((1 + f) * V9**2 - V0**2) + 0.5 * m_dot_s * (V19**2 - V0**2)  # puissance du cycle
    W_chim = m_dot_k * Pk                                                        # puissance chimique

    eta_th = W_cy / W_chim        # rendement thermique
    eta_pr = W_pr / W_cy          # rendement propulsif
    eta_global = eta_th * eta_pr  # rendement global

    Cs = m_dot_k / F_target * 3600 / 10  # kg/(h·daN)
    f_sp = F_target / m_dot_total         # m/s, poussée spécifique

    results = {
        # Conditions de vol
        'M0': M0, 'V0': V0, 'T0': T0, 'P0': P0,
        'Tt0': Tt0, 'Pt0': Pt0,
        # Stations
        'Tt2': Tt2, 'Pt2': Pt2,
        'pi_cBP': pi_cBP,
        'Tt_fan': Tt_fan, 'Pt_fan': Pt_fan,
        'Tt21': Tt21, 'Pt21': Pt21,
        'Tt25': Tt25, 'Pt25': Pt25,
        'Tt3': Tt3, 'Pt3': Pt3,
        'Tt4': Tt4, 'Pt4': Pt4,
        'Tt45': Tt45, 'Pt45': Pt45,
        'Tt5': Tt5, 'Pt5': Pt5,
        'Tt9': Tt9, 'Pt9': Pt9, 'V9': V9,
        'Tt19': Tt19, 'Pt19': Pt19, 'V19': V19,
        # Débits
        'f': f, 'm_dot_p': m_dot_p, 'm_dot_s': m_dot_s,
        'm_dot_total': m_dot_total, 'm_dot_k': m_dot_k,
        # Dimensionnement
        'A_fan': A_fan, 'r_max': r_max, 'd_max': d_max,
        # Performances
        'Cs': Cs, 'f_sp': f_sp,
        'eta_th': eta_th, 'eta_pr': eta_pr, 'eta_global': eta_global,
        'W_pr': W_pr, 'W_cy': W_cy, 'W_chim': W_chim,
        # Paramètres d'entrée
        'BPR': BPR, 'pi_f': pi_f, 'OPR': OPR, 'pi_CHP': pi_CHP,
    }
    return results


def print_results(res):
    """Affiche les résultats principaux du cycle."""
    print("=" * 60)
    print("  RÉSULTATS DU CYCLE TURBOFAN")
    print("=" * 60)

    print(f"\n--- Conditions de vol ---")
    print(f"  M0 = {res['M0']:.2f},  V0 = {res['V0']:.1f} m/s")
    print(f"  Tt0 = {res['Tt0']:.1f} K,  Pt0 = {res['Pt0']:.0f} Pa")

    print(f"\n--- Taux de compression ---")
    print(f"  π_f = {res['pi_f']:.2f},  π_cBP = {res['pi_cBP']:.3f},  π_CHP = {res['pi_CHP']:.1f}")
    print(f"  OPR = {res['OPR']:.1f}")

    print(f"\n--- Températures totales (K) ---")
    print(f"  Tt2  = {res['Tt2']:.1f}")
    print(f"  Tt_fan = {res['Tt_fan']:.1f} (sortie fan)")
    print(f"  Tt25 = {res['Tt25']:.1f} (sortie compresseur BP)")
    print(f"  Tt3  = {res['Tt3']:.1f} (sortie compresseur HP)")
    print(f"  Tt4  = {res['Tt4']:.1f} (sortie chambre combustion)")
    print(f"  Tt45 = {res['Tt45']:.1f} (sortie turbine HP)")
    print(f"  Tt5  = {res['Tt5']:.1f} (sortie turbine BP)")

    print(f"\n--- Vitesses d'éjection ---")
    print(f"  V9  = {res['V9']:.1f} m/s (tuyère primaire)")
    print(f"  V19 = {res['V19']:.1f} m/s (tuyère secondaire)")

    print(f"\n--- Débits ---")
    print(f"  ṁp = {res['m_dot_p']:.2f} kg/s (primaire)")
    print(f"  ṁs = {res['m_dot_s']:.2f} kg/s (secondaire)")
    print(f"  ṁ_total = {res['m_dot_total']:.2f} kg/s")
    print(f"  ṁk = {res['m_dot_k']:.4f} kg/s (carburant)")
    print(f"  f = ṁk/ṁp = {res['f']:.5f}")

    print(f"\n--- Dimensionnement ---")
    print(f"  Diamètre fan = {res['d_max']:.3f} m")
    print(f"  Contrainte < 2 m : {'OK' if res['d_max'] < 2.0 else 'NON RESPECTÉE'}")

    print(f"\n--- Performances ---")
    print(f"  Cs   = {res['Cs']:.4f} kg/(h·daN)")
    print(f"  f_sp = {res['f_sp']:.1f} m/s")
    print(f"  η_th = {res['eta_th']:.4f} ({res['eta_th']*100:.2f}%)")
    print(f"  η_pr = {res['eta_pr']:.4f} ({res['eta_pr']*100:.2f}%)")
    print(f"  η    = {res['eta_global']:.4f} ({res['eta_global']*100:.2f}%)")
    print("=" * 60)


# =============================================================================
# Études paramétriques
# =============================================================================

def study_eta_th_vs_OPR(Tt4_values=[1400, 1600, 1800, 2000]):
    """Q2 : rendement thermique en fonction de l'OPR pour plusieurs Tt4."""
    OPR_range = np.linspace(10, 80, 100)
    fig, ax = plt.subplots(figsize=(8, 5))
    for Tt4 in Tt4_values:
        eta_th_list = []
        for opr in OPR_range:
            try:
                res = turbofan_cycle(OPR=opr, Tt4=Tt4)
                eta_th_list.append(res['eta_th'] if res else np.nan)
            except Exception:
                eta_th_list.append(np.nan)
        ax.plot(OPR_range, eta_th_list, label=f'$T_{{t4}}$ = {Tt4} K')
    ax.set_xlabel('OPR ($\\pi_c$)')
    ax.set_ylabel('Rendement thermique $\\eta_{th}$')
    ax.set_title('Rendement thermique vs OPR')
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig('figures/eta_th_vs_OPR.png', dpi=200)
    return fig


def study_eta_global_vs_BPR(pi_f_values=[1.2, 1.45, 1.7, 2.0]):
    """Q3 : rendement global en fonction du BPR pour plusieurs π_f."""
    BPR_range = np.linspace(2, 20, 100)
    fig, ax = plt.subplots(figsize=(8, 5))
    for pf in pi_f_values:
        eta_list = []
        for bpr in BPR_range:
            try:
                res = turbofan_cycle(BPR=bpr, pi_f=pf)
                eta_list.append(res['eta_global'] if res else np.nan)
            except Exception:
                eta_list.append(np.nan)
        ax.plot(BPR_range, eta_list, label=f'$\\pi_f$ = {pf}')
    ax.set_xlabel('BPR ($\\lambda$)')
    ax.set_ylabel('Rendement global $\\eta$')
    ax.set_title('Rendement global vs BPR')
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig('figures/eta_global_vs_BPR.png', dpi=200)
    return fig


def study_mach_effect(M0_values=[0.7, 0.8]):
    """Q5 : effet du Mach de vol sur la taille et les performances."""
    print("\n" + "=" * 60)
    print("  COMPARAISON M0 = 0.7 vs M0 = 0.8")
    print("=" * 60)
    for M0 in M0_values:
        res = turbofan_cycle(M0=M0)
        print(f"\n  M0 = {M0}")
        print(f"    Diamètre fan  = {res['d_max']:.3f} m")
        print(f"    ṁ_total       = {res['m_dot_total']:.2f} kg/s")
        print(f"    η_th          = {res['eta_th']*100:.2f}%")
        print(f"    η_pr          = {res['eta_pr']*100:.2f}%")
        print(f"    η_global      = {res['eta_global']*100:.2f}%")
        print(f"    Cs            = {res['Cs']:.4f} kg/(h·daN)")
        print(f"    f_sp          = {res['f_sp']:.1f} m/s")


# =============================================================================
# Exécution
# =============================================================================
if __name__ == "__main__":
    # --- Cas de référence LEAP-1A ---
    res = turbofan_cycle()
    print_results(res)

    # --- Études paramétriques ---
    fig1 = study_eta_th_vs_OPR()
    fig2 = study_eta_global_vs_BPR()
    study_mach_effect()

    plt.show()