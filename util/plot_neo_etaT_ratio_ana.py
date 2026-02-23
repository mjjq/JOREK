#!/usr/bin/env python3
"""Compute and plot the analytical neo eta(T) ratio used in JOREK.

Implements the Fortran routine:
  get_neo_etaT_ratio_ana1(psin)

Formula:
  ratio = 1 / (1 - sqrt(sqrt(psin) * neo_minor_radius / neo_major_radius))**2

For sqrt(psin) <= XL1NCETA and XL1NCETA > 0, a smooth core blend is applied:
  aux_weight = 0.5 * (1 - cos(pi * sqrt(psin) / XL1NCETA))
  ratio = (1 - aux_weight) * 1.0 + aux_weight * ratio

Example usage:
  plot_neo_etaT_ratio_ana.py --neo-minor-radius 0.22 --neo-major-radius 1.7 --xl1nceta 0.4 --in-data jorek_eta_neo_ratio_ana_unif_psi
"""


import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def get_neo_etaT_ratio_ana1(psin: np.ndarray, neo_minor_radius: float, neo_major_radius: float, xl1nceta: float) -> np.ndarray:
    psin = np.asarray(psin, dtype=float)
    sqrt_psin = np.sqrt(np.clip(psin, 0.0, None))

    inside = np.sqrt(np.clip(sqrt_psin * neo_minor_radius / neo_major_radius, 0.0, None))
    denom = 1.0 - inside
    denom = np.where(np.abs(denom) < 1e-14, np.sign(denom) * 1e-14 + (denom == 0.0) * 1e-14, denom)

    ratio = 1.0 / (denom**2)

    if xl1nceta > 0.0:
        mask = sqrt_psin <= xl1nceta
        aux_weight = np.zeros_like(psin)
        aux_weight[mask] = 0.5 * (1.0 - np.cos(np.pi * sqrt_psin[mask] / xl1nceta))
        ratio[mask] = (1.0 - aux_weight[mask]) * 1.0 + aux_weight[mask] * ratio[mask]

    return ratio


def smooth_transition_c2(x: np.ndarray, y1: np.ndarray, y2: np.ndarray, xl1nceta: float) -> np.ndarray:
    """Smoothly blend y1->y2 with C2 continuity in [0, XL1NCETA]."""
    y2_smooth = np.array(y2, dtype=float, copy=True)
    if xl1nceta <= 0.0:
      return y2_smooth

    for i in range(len(x)):
      if x[i] <= xl1nceta:
        t = x[i] / xl1nceta
        weight = t**3 * (6.0 * t**2 - 15.0 * t + 10.0)
        y2_smooth[i] = (1.0 - weight) * y1[i] + weight * y2[i]

    return y2_smooth


def get_neo_etaT_ratio_ana2(psin: np.ndarray, neo_minor_radius: float, neo_major_radius: float, xl1nceta: float) -> np.ndarray:
    """Alternative analytical neo eta(T) ratio with C2 smooth transition."""
    psin = np.asarray(psin, dtype=float)
    sqrt_psin = np.sqrt(np.clip(psin, 0.0, None))

    inside = np.sqrt(np.clip(sqrt_psin * neo_minor_radius / neo_major_radius, 0.0, None))
    denom = 1.0 - inside
    denom = np.where(np.abs(denom) < 1e-14, np.sign(denom) * 1e-14 + (denom == 0.0) * 1e-14, denom)

    y1 = np.ones_like(psin)
    y2 = 1.0 / (denom**2)
    return smooth_transition_c2(sqrt_psin, y1, y2, xl1nceta)


def get_neo_etaT_ratio_ana2_d2_analytic(psin: np.ndarray, neo_minor_radius: float, neo_major_radius: float, xl1nceta: float) -> np.ndarray:
    psin = np.asarray(psin, dtype=float)
    d2 = np.full_like(psin, np.nan, dtype=float)

    positive = psin > 0.0
    if not np.any(positive):
      return d2

    x = psin[positive]
    sqrt_x = np.sqrt(x)

    c = np.sqrt(neo_minor_radius / neo_major_radius)
    u = c * np.power(x, 0.25)
    one_minus_u = 1.0 - u

    f = 1.0 / (one_minus_u**2)
    u_prime = u / (4.0 * x)
    u_second = -3.0 * u / (16.0 * x**2)
    f_prime = 2.0 * u_prime / (one_minus_u**3)
    f_second = 2.0 * u_second / (one_minus_u**3) + 6.0 * (u_prime**2) / (one_minus_u**4)

    y_second = f_second.copy()

    if xl1nceta > 0.0:
      mask_core = sqrt_x <= xl1nceta
      if np.any(mask_core):
        x_core = x[mask_core]
        sqrt_x_core = sqrt_x[mask_core]
        f_core = f[mask_core]
        f_prime_core = f_prime[mask_core]
        f_second_core = f_second[mask_core]

        t = sqrt_x_core / xl1nceta
        w = t**3 * (6.0 * t**2 - 15.0 * t + 10.0)
        w_t = 30.0 * t**2 * (t - 1.0)**2
        w_tt = 120.0 * t**3 - 180.0 * t**2 + 60.0 * t

        dt_dx = 1.0 / (2.0 * xl1nceta * sqrt_x_core)
        d2t_dx2 = -1.0 / (4.0 * xl1nceta * x_core**1.5)

        w_prime = w_t * dt_dx
        w_second = w_tt * dt_dx**2 + w_t * d2t_dx2

        y_second[mask_core] = w_second * (f_core - 1.0) + 2.0 * w_prime * f_prime_core + w * f_second_core

    d2[positive] = y_second
    return d2


def get_neo_etaT_ratio_ana1_d2_analytic(psin: np.ndarray, neo_minor_radius: float, neo_major_radius: float, xl1nceta: float) -> np.ndarray:
    psin = np.asarray(psin, dtype=float)
    d2 = np.full_like(psin, np.nan, dtype=float)

    positive = psin > 0.0
    if not np.any(positive):
      return d2

    x = psin[positive]
    sqrt_x = np.sqrt(x)
    c = np.sqrt(neo_minor_radius / neo_major_radius)
    u = c * np.power(x, 0.25)
    one_minus_u = 1.0 - u

    f = 1.0 / (one_minus_u**2)
    u_prime = u / (4.0 * x)
    u_second = -3.0 * u / (16.0 * x**2)
    f_prime = 2.0 * u_prime / (one_minus_u**3)
    f_second = 2.0 * u_second / (one_minus_u**3) + 6.0 * (u_prime**2) / (one_minus_u**4)

    g_second = f_second.copy()

    if xl1nceta > 0.0:
      mask_core = sqrt_x <= xl1nceta
      if np.any(mask_core):
        x_core = x[mask_core]
        sqrt_x_core = sqrt_x[mask_core]
        f_core = f[mask_core]
        f_prime_core = f_prime[mask_core]
        f_second_core = f_second[mask_core]

        arg = np.pi * sqrt_x_core / xl1nceta
        w = 0.5 * (1.0 - np.cos(arg))
        w_prime = np.pi / (4.0 * xl1nceta * sqrt_x_core) * np.sin(arg)
        w_second = (
          -np.pi / (8.0 * xl1nceta) * np.sin(arg) / (x_core**1.5)
          + np.pi**2 / (8.0 * xl1nceta**2) * np.cos(arg) / x_core
        )

        g_second[mask_core] = w_second * (f_core - 1.0) + 2.0 * w_prime * f_prime_core + w * f_second_core

    d2[positive] = g_second
    return d2


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute and plot get_neo_etaT_ratio_ana1.")
    parser.add_argument("--neo-minor-radius", type=float, required=True, help="neo_minor_radius")
    parser.add_argument("--neo-major-radius", type=float, required=True, help="neo_major_radius")
    parser.add_argument("--xl1nceta", type=float, required=True, help="XL1NCETA")
    parser.add_argument("--psin-min", type=float, default=1e-5, help="Minimum Psi_N")
    parser.add_argument("--psin-max", type=float, default=1.1, help="Maximum Psi_N (default: 1.1)")
    parser.add_argument("--npoints", type=int, default=550, help="Number of points (default: 550)")
    parser.add_argument("--out-data", type=Path, default=Path("jorek_eta_neo_ratio_ana_py.dat"), help="Output data file")
    parser.add_argument("--out-fig", type=Path, default=Path("jorek_eta_neo_ratio_ana_py.png"), help="Output figure file")
    parser.add_argument(
      "--in-data",
      type=Path,
      default=None,
      help="Optional input file with two columns: Psi_N and neo_etaT_ratio_ana.",
    )
    args = parser.parse_args()

    psin = np.linspace(args.psin_min, args.psin_max, args.npoints)

    # ratio = get_neo_etaT_ratio_ana1(psin, args.neo_minor_radius, args.neo_major_radius, args.xl1nceta)
    # second_derivative_analytic = get_neo_etaT_ratio_ana1_d2_analytic(psin, args.neo_minor_radius, args.neo_major_radius, args.xl1nceta)

    ratio = get_neo_etaT_ratio_ana2(psin, args.neo_minor_radius, args.neo_major_radius, args.xl1nceta)
    second_derivative_analytic = get_neo_etaT_ratio_ana2_d2_analytic(psin, args.neo_minor_radius, args.neo_major_radius, args.xl1nceta)

    first_derivative = np.gradient(ratio, psin, edge_order=2)
    second_derivative = np.gradient(first_derivative, psin, edge_order=2)

    np.savetxt(
      args.out_data,
      np.column_stack([psin, ratio, second_derivative, second_derivative_analytic]),
      header="Psi_N neo_etaT_ratio_ana d2_neo_etaT_ratio_ana_dPsiN2_num d2_neo_etaT_ratio_ana_dPsiN2_analytic",
      fmt="%.16e",
    )

    fig, axes = plt.subplots(2, 1, figsize=(7, 8), sharex=True)

    psin_num = None
    ratio_num = None
    second_derivative_num = None
    if args.in_data is not None:
      num_data = np.loadtxt(args.in_data)
      if num_data.ndim == 1:
        num_data = num_data.reshape(1, -1)
      if num_data.shape[1] < 2:
        raise ValueError("Input file must contain at least two columns: Psi_N and neo_etaT_ratio_ana.")
      psin_num = num_data[:, 0]
      ratio_num = num_data[:, 1]
      first_derivative_num = np.gradient(ratio_num, psin_num, edge_order=2)
      second_derivative_num = np.gradient(first_derivative_num, psin_num, edge_order=2)

    axes[0].plot(psin, ratio, lw=2, label="Analytical")
    if psin_num is not None:
      axes[0].plot(psin_num, ratio_num, lw=2, linestyle="--", color = 'k', label="From file")
    axes[0].set_ylabel("neo_etaT_ratio_ana")
    axes[0].set_title("Analytical neo eta(T) ratio")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].plot(psin, second_derivative, lw=2, label="numerical ddpsi of analytical ratio")
    axes[1].plot(psin, second_derivative_analytic, lw=2, linestyle=":", label="analytical ddpsi of analytical ratio")
    if psin_num is not None:
      axes[1].plot(psin_num, second_derivative_num, lw=2, linestyle="--", color = 'k', label="From file")
    axes[1].set_xlabel("Psi_N")
    axes[1].set_ylabel("d2/dPsi_N2")
    axes[1].set_title("Second derivative of neo eta(T) ratio")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(args.out_fig, dpi=150)
    plt.show()

    print(f"Saved data: {args.out_data}")
    print(f"Saved figure: {args.out_fig}")


if __name__ == "__main__":
    main()
