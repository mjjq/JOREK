#!/usr/bin/env python3
"""Compare JxB_radial and gradP_radial Fourier outputs."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


dec_Jxb_gradp_ratio = 1e-2
gamma = 0.0085

n = 1
m = 1

step = 600
xlim = 0.6

def load_four_file(path: Path) -> np.ndarray:
    """Load four_results file, skipping the first row header."""
    return np.loadtxt(path, skiprows=1)

def main() -> None:
    jxb_path   = Path(f"four_results/JxB_radial_n{n:03d}_m{m:03d}_step{step:06d}")
    gradp_path = Path(f"four_results/gradP_radial_n{n:03d}_m{m:03d}_step{step:06d}")
    Ur_path    = Path(f"four_results/U_radial_n{n:03d}_m{m:03d}_step{step:06d}")
    rho_path   = Path(f"four_results/rho_n{n:03d}_m{m:03d}_step{step:06d}")
    psi_path   = Path(f"four_results/Psi_n{n:03d}_m{m:03d}_step{step:06d}")
    # UR_path    = Path(f"four_results/UR_n{n:03d}_m{m:03d}_step{step:06d}")
    # UZ_path    = Path(f"four_results/UZ_n{n:03d}_m{m:03d}_step{step:06d}")
    # Uphi_path  = Path(f"four_results/Uphi_n{n:03d}_m{m:03d}_step{step:06d}")
    RUxB_phi_path = Path(f"four_results/RUxB_phi_n{n:03d}_m{m:03d}_step{step:06d}")
    rho_divU_path = Path(f"four_results/rho_divU_n{n:03d}_m{m:03d}_step{step:06d}")
    U_gradRho_path = Path(f"four_results/U_gradRho_n{n:03d}_m{m:03d}_step{step:06d}")
    
    
    rhoeq_path = Path(f"four_results/rho_n{0:03d}_m{0:03d}_step{step:06d}")
    jxbeq_path =  Path(f"four_results/JxB_radial_n{0:03d}_m{0:03d}_step{step:06d}")
    gradpeq_path = Path(f"four_results/gradP_radial_n{0:03d}_m{0:03d}_step{step:06d}")

    # BReq_path = Path(f"four_results/BR_n{0:03d}_m{0:03d}_step{step:06d}")
    # BZeq_path = Path(f"four_results/BZ_n{0:03d}_m{0:03d}_step{step:06d}")
    # Bphi_eq_path = Path(f"four_results/Bphi_n{0:03d}_m{0:03d}_step{step:06d}")


    jxb = load_four_file(jxb_path)
    gradp = load_four_file(gradp_path)
    Ur = load_four_file(Ur_path)
    rho = load_four_file(rho_path)
    psi = load_four_file(psi_path)
    # UR = load_four_file(UR_path)
    # UZ = load_four_file(UZ_path)
    # Uphi = load_four_file(Uphi_path)
    RUxB_phi = load_four_file(RUxB_phi_path)
    rho_divU = load_four_file(rho_divU_path)
    U_gradRho = load_four_file(U_gradRho_path)

    rho_eq = load_four_file(rhoeq_path)
    jxb_eq = load_four_file(jxbeq_path)
    gradp_eq = load_four_file(gradpeq_path)
    # BReq = load_four_file(BReq_path)
    # BZeq = load_four_file(BZeq_path)
    # Bphi_eq = load_four_file(Bphi_eq_path)

    if jxb.shape[1] < 5 or gradp.shape[1] < 5:
        raise ValueError("Expected at least 5 columns: Psi_N, abs, real, imag, phase.")
    if jxb.shape[0] != gradp.shape[0]:
        raise ValueError("Input files have different row counts.")

    psi_n = jxb[:, 0]
    rhop = np.sqrt(np.clip(psi_n, a_min=0.0, a_max=None))
    mask = rhop <= xlim

    real_diff = jxb[:, 2] - gradp[:, 2]
    imag_diff = jxb[:, 3] - gradp[:, 3]

    diff_eq = jxb_eq[:, 1] - gradp_eq[:, 1]

    # real_UxB_phi = UR[:, 2] * BZeq[:, 1] - UZ[:, 2] * BReq[:, 1]
    # imag_UxB_phi = UR[:, 3] * BZeq[:, 1] - UZ[:, 3] * BReq[:, 1]
    real_RUxB_phi = RUxB_phi[:, 2]
    imag_RUxB_phi = RUxB_phi[:, 3]

    fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=False)
    axes = axes.ravel()

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    diff_color = colors[0]
    ur_color = colors[1]
    jxb_color = colors[2]
    gradp_color = colors[3]
    diff_rho_color = colors[4]

    axes[0].plot(rhop[mask], real_diff[mask] / rho_eq[mask, 1] * dec_Jxb_gradp_ratio, linestyle="-", label=rf"{dec_Jxb_gradp_ratio} $Re(\delta(J\times B) - \nabla\delta P)/\rho_0$", color=diff_color)
    axes[0].plot(rhop[mask], imag_diff[mask] / rho_eq[mask, 1] * dec_Jxb_gradp_ratio, linestyle="--", label=rf"{dec_Jxb_gradp_ratio} $Im(\delta(J\times B) - \nabla\delta P)/\rho_0$", color=diff_color)
    axes[0].plot(rhop[mask], - diff_eq[mask] * rho[mask,2] / rho_eq[mask, 1] ** 2, linestyle="-", label=r"$-(J_0xB_0 - \nabla P_0)*Re(\delta \rho)/\rho_0^2$", color=diff_rho_color)
    axes[0].plot(rhop[mask], - diff_eq[mask] * rho[mask,3] / rho_eq[mask, 1] ** 2, linestyle="--", label=r"$-(J_0xB_0 - \nabla P_0)*Im(\delta \rho)/\rho_0^2$", color=diff_rho_color)

    # axes[0].plot(rhop[mask], real_diff[mask] / rho_eq[mask, 1] - diff_eq[mask] * rho[mask,2] / rho_eq[mask, 1] ** 2, linestyle="-", label=r"Sum with RHS terms", color=diff_color)
    # axes[0].plot(rhop[mask], imag_diff[mask] / rho_eq[mask, 1] - diff_eq[mask] * rho[mask,3] / rho_eq[mask, 1] ** 2, linestyle="--", label=r"Sum with RHS terms", color=diff_color)

    axes[0].plot(rhop[mask], gamma * Ur[mask, 2], label=fr"$\gamma$ * Re(Ur)", linestyle="-", alpha=0.7, color=ur_color)
    axes[0].plot(rhop[mask], gamma * Ur[mask, 3], label=fr"$\gamma$ * Im(Ur)", linestyle="--", alpha=0.7, color=ur_color)
    axes[0].set_ylabel("Difference")
    axes[0].set_title(f"Momentum equation terms (m/n = {m}/{n})")
    axes[0].grid(True, alpha=0.3)
    legend0 = axes[0].legend()
    legend0.set_draggable(True)
    axes[1].set_xlim(0, xlim)

    axes[1].plot(rhop[mask], jxb[mask, 2] / rho_eq[mask, 1], linestyle="-", label=r"$Re\delta(J\times B)/\rho_0$", color=jxb_color)
    axes[1].plot(rhop[mask], jxb[mask, 3] / rho_eq[mask, 1], linestyle="--", label=r"$Im\delta(J\times B)/\rho_0$", color=jxb_color)
    axes[1].plot(rhop[mask], gradp[mask, 2] / rho_eq[mask, 1], linestyle="-", label=r"$Re(\nabla \delta P)/\rho_0$", alpha=0.7, color=gradp_color)
    axes[1].plot(rhop[mask], gradp[mask, 3] / rho_eq[mask, 1], linestyle="--", label=r"$Im(\nabla \delta P)/\rho_0$", alpha=0.7, color=gradp_color)
    axes[1].set_ylabel("Value")
    axes[1].set_title(f"JxB_radial and gradP_radial (m/n = {m}/{n})")
    axes[1].grid(True, alpha=0.3)
    legend1 = axes[1].legend()
    legend1.set_draggable(True)
    axes[1].set_xlim(0, xlim)

    axes[2].plot(rhop[mask], rho[mask, 2], linestyle="-", label=r"$Re(\delta \rho)$", color=colors[0])
    axes[2].plot(rhop[mask], rho[mask, 3], linestyle="--", label=r"$Im(\delta \rho)$", color=colors[0])
    axes[2].set_ylabel("Value")
    axes[2].set_title(rf"$\delta \rho$ (m/n = {m}/{n})")
    axes[2].grid(True, alpha=0.3)
    legend2 = axes[2].legend()
    legend2.set_draggable(True)
    axes[2].set_xlim(0, xlim)

    axes[3].plot(rhop[mask], diff_eq[mask], linestyle="-", label=r"$J_0xB_0 - \nabla P_0$", color=colors[1])
    axes[3].set_ylabel("Difference")
    axes[3].set_title(f"JxB_radial - gradP_radial at equilibrium (m/n = {0}/{0})")
    axes[3].grid(True, alpha=0.3)
    legend3 = axes[3].legend()
    legend3.set_draggable(True)
    axes[3].set_xlim(0, xlim)

    # axes[4].plot(rhop[mask], real_UxB_phi[mask], linestyle="-", label=r"$Re(\delta U x B_0)_\phi$", color=colors[0])
    # axes[4].plot(rhop[mask], imag_UxB_phi[mask], linestyle="--", label=r"$Im(\delta U x B_0)_\phi$", color=colors[0])
    axes[4].plot(rhop[mask], real_RUxB_phi[mask], linestyle="-", label=r"$Re(R \delta U x B_0)_\phi$", alpha=0.7, color=colors[0])
    axes[4].plot(rhop[mask], imag_RUxB_phi[mask], linestyle="--", label=r"$Im(R \delta U x B_0)_\phi$", alpha=0.7, color=colors[0])
    axes[4].plot(rhop[mask], gamma * psi[mask, 2], linestyle="-", label=fr"$\gamma$ * Re(Psi)", alpha=0.7, color=colors[1])
    axes[4].plot(rhop[mask], gamma * psi[mask, 3], linestyle="--", label=fr"$\gamma$ * Im(Psi)", alpha=0.7, color=colors[1])
    axes[4].set_ylabel("Value")
    axes[4].set_title(rf"Induction equation terms (m/n = {m}/{n})")
    axes[4].grid(True, alpha=0.3)
    legend4 = axes[4].legend()
    legend4.set_draggable(True)
    axes[4].set_xlim(0, xlim)

    axes[5].plot(rhop[mask], -rho_divU[mask, 2] - U_gradRho[mask, 2], linestyle="-", label=r"$Re(-\rho_0\nabla\cdot\delta \mathbf{U} - \delta \mathbf{U} \cdot \nabla \rho_0)$", color=colors[0])
    axes[5].plot(rhop[mask], -rho_divU[mask, 3] - U_gradRho[mask, 3], linestyle="--", label=r"$Im(-\rho_0\nabla\cdot\delta \mathbf{U} - \delta \mathbf{U} \cdot \nabla \rho_0)$", color=colors[0])
    axes[5].plot(rhop[mask], gamma * rho[mask, 2], linestyle="-", label=fr"$\gamma$ * Re($\delta \rho$)", color=colors[1])
    axes[5].plot(rhop[mask], gamma * rho[mask, 3], linestyle="--", label=fr"$\gamma$ * Im($\delta \rho$)", color=colors[1])
    axes[5].set_xlabel("rhop = sqrt(Psi_N)")
    axes[5].set_ylabel("Value")
    axes[5].set_title(rf"Continuity equation terms (m/n = {m}/{n})")
    axes[5].grid(True, alpha=0.3)
    legend5 = axes[5].legend()
    legend5.set_draggable(True)
    axes[5].set_xlim(0, xlim)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
