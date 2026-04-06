#!/usr/bin/env python3
"""Fit invariant-mass spectrum background with Gaussian Process Regression (GPyTorch) and signal with two-sided Crystal Ball function.
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from typing import Iterable, List, Sequence, Tuple

import numpy as np
import ROOT

import matplotlib
matplotlib.use("Agg") # Don't open any GUI while running
import matplotlib.pyplot as plt

import ROOT
import torch
import gpytorch # For Gaussian Process Regression
from scipy.optimize import curve_fit # For the signal crystal ball fit

# Same band limits as in the preliminary analysis
BAND_LIMITS = np.array([
    0.030, 0.070,  # pi0 left side
    0.080, 0.199,  # pi0 peak
    0.209, 0.249,  # pi0 right side
    0.257, 0.371,  # eta left side
    0.399, 0.739,  # eta peak
    0.767, 0.880,  # eta right side
], dtype=float)


def build_training_ranges(pi0_gate=None, eta_gate=None):
    """
    Build sideband training ranges with optional signal-gate overrides.
    """
    if pi0_gate is None:
        pi0_lo, pi0_hi = float(BAND_LIMITS[2]), float(BAND_LIMITS[3])
    else:
        if not isinstance(pi0_gate, np.ndarray):
            raise TypeError("pi0_gate must be a numpy.ndarray or None");
        if pi0_gate.shape != (2,):
            raise ValueError("pi0_gate must have shape (2,)")
        pi0_lo = max(0.0, pi0_gate[0])
        pi0_hi = min(1.0, pi0_gate[1])

    if eta_gate is None:
        eta_lo, eta_hi = float(BAND_LIMITS[8]), float(BAND_LIMITS[9])
    else:
        if not isinstance(eta_gate, np.ndarray):
            raise TypeError("eta_gate must be a numpy.ndarray or None");
        if eta_gate.shape != (2,):
            raise ValueError("eta_gate must have shape (2,)")
        eta_lo = max(0.0, eta_gate[0])
        eta_hi = min(1.0, eta_gate[1])

    ranges_pi0 = [(0.0, pi0_lo), (pi0_hi, eta_lo), (eta_hi, 1.0)] # For the pi0, fit the full background
    eta_start = 0.30
    ranges_eta = [(eta_start, eta_lo), (eta_hi, 1.0)]
    return ranges_pi0, ranges_eta, (pi0_lo, pi0_hi), (eta_lo, eta_hi) # For the eta, only fit the right side of the background


@dataclass
class SelectedSet:
    x: np.ndarray
    y: np.ndarray
    yerr: np.ndarray
    edges_low: np.ndarray
    edges_high: np.ndarray
    width: np.ndarray
    mask: np.ndarray
    idx: np.ndarray


def extract_training_set(h: "ROOT.TH1", intervals: Sequence[Tuple[float, float]]) -> SelectedSet:
    axis = h.GetXaxis()
    n = h.GetNbinsX()
    idx = np.arange(1, n + 1, dtype=int)
    low = np.array([axis.GetBinLowEdge(int(i)) for i in idx], dtype=float)
    high = np.array([axis.GetBinUpEdge(int(i)) for i in idx], dtype=float)
    x = np.array([axis.GetBinCenter(int(i)) for i in idx], dtype=float)
    y = np.array([h.GetBinContent(int(i)) for i in idx], dtype=float)
    yerr = np.array([h.GetBinError(int(i)) for i in idx], dtype=float)
    width = high - low

    mask = np.zeros_like(x, dtype=bool)
    for xmin, xmax in intervals:
        mask |= (high > xmin) & (low < xmax)

    return SelectedSet(
        x=x[mask], y=y[mask], yerr=yerr[mask],
        edges_low=low[mask], edges_high=high[mask], width=width[mask],
        mask=mask, idx=idx[mask]
    )

class ExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, min_lengthscale: float | None = None):
        super().__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()

        base_kernel = gpytorch.kernels.RBFKernel()
        if min_lengthscale is not None:
            
            base_kernel.register_constraint(
                "raw_lengthscale",
                gpytorch.constraints.GreaterThan(float(min_lengthscale))
            )
            
        self.covar_module = gpytorch.kernels.ScaleKernel(base_kernel)

        if min_lengthscale is not None:
            self.covar_module.base_kernel.lengthscale = max(float(min_lengthscale), 0.5)

    def forward(self, x):
        return gpytorch.distributions.MultivariateNormal(self.mean_module(x), self.covar_module(x))


@dataclass
class GPFitResult:
    test_x: torch.Tensor
    mean: torch.Tensor
    lo: torch.Tensor
    hi: torch.Tensor
    change_point: float | None = None


def _prepare_tensors(sel: SelectedSet, device: torch.device, dtype: torch.dtype):
    x_t = torch.from_numpy(sel.x).to(device=device, dtype=dtype)
    y_t = torch.from_numpy(sel.y).to(device=device, dtype=dtype)
    yerr_np = sel.yerr
    yerr_t = torch.from_numpy(yerr_np).to(device=device, dtype=dtype)
    train_x = x_t.unsqueeze(-1)
    # Normalize and center the data before the fit
    y_mean = y_t.mean()
    y_std = y_t.std().clamp_min(1e-12)
    train_y = (y_t - y_mean) / y_std
    obs_var = ((yerr_t ** 2) / (y_std ** 2)).clamp_min((1e-6) ** 2)
    return train_x, y_t, yerr_t, train_y, obs_var, y_mean, y_std

def fit_gp(sel: SelectedSet, n_iter: int, lr: float, device: torch.device, dtype: torch.dtype, quiet: bool,
           min_lengthscale: float | None = 1.0) -> GPFitResult:
    train_x, y_t, yerr_t, train_y, obs_var, y_mean, y_std = _prepare_tensors(sel, device, dtype)
    likelihood = gpytorch.likelihoods.FixedNoiseGaussianLikelihood(noise=obs_var)
    model = ExactGPModel(train_x, train_y, likelihood, min_lengthscale=min_lengthscale).to(device)
    print("lengthscale constraint:",
      model.covar_module.base_kernel.raw_lengthscale_constraint)
    print("initial lengthscale:",
      model.covar_module.base_kernel.lengthscale.item())

    model.train(); likelihood.train()
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

    for i in range(n_iter):
        opt.zero_grad()
        out = model(train_x)
        loss = -mll(out, train_y)
        loss.backward()
        opt.step()
        if (i % 50 == 0) and (not quiet):
            length = float(model.covar_module.base_kernel.lengthscale.item())
            oscale = float(model.covar_module.outputscale.item())
            print(f"[GP] iter {i:4d} loss={loss.item():.4f} length={length:.4f} os={oscale:.4f}")

    model.eval(); likelihood.eval()
    x0 = float(train_x.squeeze()[0].item())
    test_x = torch.linspace(x0, 1, 201, dtype=dtype, device=device).unsqueeze(-1)
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        pred = likelihood(model(test_x))
    mean_n = pred.mean
    lo_n, hi_n = pred.confidence_region()
    mean = mean_n * y_std + y_mean
    lo = lo_n * y_std + y_mean
    hi = hi_n * y_std + y_mean
    return GPFitResult(test_x=test_x.detach().cpu(), mean=mean.detach().cpu(), lo=lo.detach().cpu(), hi=hi.detach().cpu(), change_point=None)


def crystalball_two_sides(x, A, mu, sigma, alphaL, nL, alphaR, nR):
    x = np.asarray(x, dtype=float)
    s = np.clip(sigma, 1e-6, np.inf)
    alphaL = np.clip(alphaL, 1e-3, 50.0)
    nL = np.clip(nL, 1.01, 50.0)
    alphaR = np.clip(alphaR, 1e-3, 50.0)
    nR = np.clip(nR, 1.01, 50.0)
    z = (x - mu) / s
    gauss = np.exp(-0.5 * z**2)
    cb = gauss.copy()

    mask_left = z < -alphaL
    if np.any(mask_left):
        zt = z[mask_left]
        BL = nL / alphaL - alphaL
        logAL = nL * (np.log(nL) - np.log(alphaL)) - 0.5 * alphaL**2
        log_tail = np.clip(logAL - nL * np.log(BL - zt), -700, 700)
        cb[mask_left] = np.exp(log_tail)

    mask_right = z > alphaR
    if np.any(mask_right):
        zt = z[mask_right]
        BR = nR / alphaR - alphaR
        logAR = nR * (np.log(nR) - np.log(alphaR)) - 0.5 * alphaR**2
        log_tail = np.clip(logAR - nR * np.log(BR + zt), -700, 700)
        cb[mask_right] = np.exp(log_tail)

    return A * cb / np.sqrt(2 * np.pi * s**2)


def two_cb(x, a1, mu1, sigma1, alphaL1, nL1, alphaR1, nR1,
              a2, mu2, sigma2, alphaL2, nL2, alphaR2, nR2):
    return (crystalball_two_sides(x, a1, mu1, sigma1, alphaL1, nL1, alphaR1, nR1) +
            crystalball_two_sides(x, a2, mu2, sigma2, alphaL2, nL2, alphaR2, nR2))


def fit_signal_cb(full_x: np.ndarray, residuals: np.ndarray):
    # Initial amplitudes from residual integrals in peak bands
    x_t = torch.from_numpy(full_x.astype(np.float64))
    r_t = torch.from_numpy(residuals.astype(np.float64))

    m1 = (x_t > BAND_LIMITS[2]) & (x_t < BAND_LIMITS[3])
    m2 = (x_t > BAND_LIMITS[8]) & (x_t < BAND_LIMITS[9])
    A_pi0 = float(torch.trapz(r_t[m1], x_t[m1])) if int(m1.sum()) > 1 else 1.0
    A_eta = float(torch.trapz(r_t[m2], x_t[m2])) if int(m2.sum()) > 1 else 1.0

    init = [A_pi0, 0.135, 0.015, 1, 10, 1, 10,
            A_eta, 0.55, 0.03, 1, 10, 1, 10]

    # Loose bounds to stabilize fit
    lb = [0,    0.10, 0.001, 0.05, 1.01, 0.05, 1.01,
          0,    0.45, 0.001, 0.05, 1.01, 0.05, 1.01]
    ub = [np.inf, 0.18, 0.08,  10.0, 80.0, 10.0, 80.0,
          np.inf, 0.65, 0.12,  10.0, 80.0, 10.0, 80.0]

    popt, _ = curve_fit(two_cb, full_x, residuals, p0=init, bounds=(lb, ub), maxfev=20000)
    return popt, two_cb(full_x, *popt)


def read_histogram(root_path: str, hist_name: str):
    f = ROOT.TFile.Open(root_path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {root_path}")
    h = f.Get(hist_name)
    if h is None:
        raise RuntimeError(f"Histogram '{hist_name}' not found in {root_path}")
    print(hist_name)
    # Detach from file so file can close safely
    h.SetDirectory(0)
    f.Close()
    return h



def parse_cb_peaks_from_popt(popt: np.ndarray):
    """Extract symmetric peak parameters from two-CB fit array."""
    p = np.asarray(popt, dtype=float)
    return {
        "pi0": {"yield": float(p[0]), "mu": float(p[1]), "sigma": abs(float(p[2]))},
        "eta": {"yield": float(p[7]), "mu": float(p[8]), "sigma": abs(float(p[9]))},
    }


def run_full_fit(h, *, device, dtype, args, pi0_gate=None, eta_gate=None):
    """Single-pass fit returning curves + peak fit result dict."""
    ranges_pi0, ranges_eta, pi0_gate, eta_gate = build_training_ranges(
        pi0_gate=pi0_gate,
        eta_gate=eta_gate,
    )
    sel_pi0 = extract_training_set(h, ranges_pi0)
    sel_eta = extract_training_set(h, ranges_eta)
    full_set = extract_training_set(h, [(0.0, 1.0)])

    res_pi0 = fit_gp(sel_pi0, n_iter=args.iters_pi0, lr=args.lr, device=device, dtype=dtype, quiet=not args.verbose,
                     min_lengthscale=args.min_pi0_lengthscale)
    eta_min_ls = None if (args.min_eta_lengthscale is not None and args.min_eta_lengthscale <= 0) else args.min_eta_lengthscale
    res_eta = fit_gp(sel_eta, n_iter=args.iters_eta, lr=args.lr, device=device, dtype=dtype, quiet=not args.verbose,
                     min_lengthscale=eta_min_ls)

    x_pi0 = res_pi0.test_x.squeeze().numpy()
    x_eta = res_eta.test_x.squeeze().numpy()
    m_pi0 = res_pi0.mean.squeeze().numpy()
    lo_pi0 = res_pi0.lo.squeeze().numpy()
    hi_pi0 = res_pi0.hi.squeeze().numpy()
    m_eta = res_eta.mean.squeeze().numpy()
    lo_eta = res_eta.lo.squeeze().numpy()
    hi_eta = res_eta.hi.squeeze().numpy()

    # Merge the two backgrounds
    mask_pi0 = x_pi0 <= 0.30
    mask_eta = x_eta > 0.30
    full_x_bkg = np.concatenate([x_pi0[mask_pi0], x_eta[mask_eta]])
    mean_bkg_curve = np.concatenate([m_pi0[mask_pi0], m_eta[mask_eta]])
    lo_bkg_curve = np.concatenate([lo_pi0[mask_pi0], lo_eta[mask_eta]])
    hi_bkg_curve = np.concatenate([hi_pi0[mask_pi0], hi_eta[mask_eta]])

    # Full curve
    full_x = full_set.x.astype(np.float64)
    full_y = full_set.y.astype(np.float64)
    full_yerr = full_set.yerr.astype(np.float64)

    bkg_mean = np.interp(full_x, full_x_bkg, mean_bkg_curve)
    residuals = full_y - bkg_mean # signal
    popt, fitted_signal = fit_signal_cb(full_x, residuals)
    peaks = parse_cb_peaks_from_popt(popt)

    return {
        "full_x": full_x, "full_y": full_y, "full_yerr": full_yerr,
        "bkg_mean": bkg_mean, "fitted_signal": fitted_signal,
        "full_x_bkg": full_x_bkg, "lo_bkg": lo_bkg_curve, "hi_bkg": hi_bkg_curve,
        "popt": popt, "peaks": peaks, "pi0_gate": pi0_gate, "eta_gate": eta_gate,
        "pi0_change_point": res_pi0.change_point,
    }


def build_informed_gate_from_previous_pass(pass1, k_pi0=4.0, k_eta=4.0):
    peaks = pass1["peaks"]
    pi0_gate = np.array([peaks["pi0"]["mu"] - float(k_pi0)*peaks["pi0"]["sigma"], \
                         peaks["pi0"]["mu"] + float(k_pi0)*peaks["pi0"]["sigma"]],
                        dtype=float)
    eta_gate = np.array([peaks["eta"]["mu"] - float(k_eta)*peaks["eta"]["sigma"], \
                         peaks["eta"]["mu"] + float(k_eta)*peaks["eta"]["sigma"]],
                        dtype=float)
    return pi0_gate, eta_gate

def save_root_graphs(full_x, bkg_mean, fitted_signal, graph_template: str, output_path: str):

    # Length check
    if len(full_x) != len(fitted_signal):
        raise ValueError(
            f"Length mismatch for signal graph: len(full_x)={len(full_x)} "
            f"!= len(fitted_signal)={len(fitted_signal)}"
        )
    if len(full_x) != len(bkg_mean):
        raise ValueError(
            f"Length mismatch for background graph: len(full_x)={len(full_x)} "
            f"!= len(bkg_mean)={len(bkg_mean)}"
        )

    # Create ROOT graphs
    signal_name=f"{graph_template}_signal";
    bkg_name=f"{graph_template}_bkg";

    graph_signal = ROOT.TGraph(len(full_x), full_x, fitted_signal)
    graph_signal.SetName(signal_name)
    graph_signal.SetTitle(";M_{#gamma#gamma} [GeV];Signal");

    graph_bkg = ROOT.TGraph(len(full_x), full_x, bkg_mean)
    graph_bkg.SetName(bkg_name)
    graph_bkg.SetTitle(";M_{#gamma#gamma} [GeV];Background");

    fout = ROOT.TFile(output_path, "RECREATE")
    if not fout or fout.IsZombie():
        raise OSError(f"Could not create ROOT file: {output_root_path}")

    fout.cd();
    graph_signal.Write();
    graph_bkg.Write();
    fout.Close()

def make_final_plot(full_x, full_y, full_yerr, bkg_mean, fitted_signal, full_x_bkg, lo_bkg, hi_bkg, output_path: str):
    fig = plt.figure(figsize=(8, 6))
    plt.errorbar(full_x, full_y, full_yerr, fmt='ko', ms=1.0, label='Observed (rate)')
    plt.plot(full_x, bkg_mean, 'g', label='GP Bkg mean')
    plt.plot(full_x, bkg_mean + fitted_signal, 'r', label='GP Bkg + Crystal Ball mean')
    plt.plot(full_x, fitted_signal, 'b', ms=1.0, label='Fitted Crystal Ball Signal')
    plt.fill_between(full_x_bkg, lo_bkg, hi_bkg, alpha=0.3, label='GP 95%')
    plt.xlabel(r"$M_{\gamma\gamma}$ [GeV]")
    plt.ylabel("Counts / [2 MeV]")
    plt.legend()
    plt.tight_layout()

    print("output_path={}".format(output_path))
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    print(f"save {output_path}.pdf")
    fig.savefig(f"{output_path}.pdf", dpi=200, bbox_inches='tight')
    fig.savefig(f"{output_path}.png", dpi=200, bbox_inches='tight')
    plt.close(fig)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Invariant-mass Gaussian Process Regression")
    p.add_argument("--input", required=True, help="Input ROOT file path")
    p.add_argument("--hist", default="h_pair_mass_pt_2", help="Histogram name")
    p.add_argument("--output", required=True, help="Output figure path (saved to png and pdf)")
    p.add_argument("--output-root", required=True, help="Output ROOT path")
    p.add_argument("--graph-template", default="graph_", help="Name of graph in ROOT file")
    p.add_argument("--iters-pi0", type=int, default=1000)
    p.add_argument("--iters-eta", type=int, default=1000)
    p.add_argument("--lr", type=float, default=0.01)
    p.add_argument("--min-pi0-lengthscale", type=float, default=None,
                   help="Minimum GP RBF lengthscale for pi0 GP (GeV).")
    p.add_argument("--min-eta-lengthscale", type=float, default=1.0,
                   help="Minimum GP RBF lengthscale for eta GP (GeV). Set <=0 to disable.")
    p.add_argument("--verbose", action="store_true", help="Print training logs")
    p.add_argument("--number-pass", type=int, default=1,
                   help="Optionally do several pass using previous pass peak fit to define signal gates.")
    p.add_argument("--k-pi0-mask", type=float, default=4.0,
                   help="Pass-2 informed pi0 gate half-width = max(k * sigma_pi0_pass1, min-halfwidth).")
    p.add_argument("--k-eta-mask", type=float, default=4.0,
                   help="Pass-2 informed eta gate half-width = max(k * sigma_eta_pass1, min-halfwidth).")
    p.add_argument("--save-pass1-overlay", default=None,
                   help="Optional path to save pass-1 overlay when using --two-pass-informed-gate.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    dtype = torch.float64

    h = read_histogram(args.input, args.hist)

    # Pass 1 (nominal signal range)
    pass1 = run_full_fit(
        h, device=device, dtype=dtype, args=args,
        pi0_gate=None,
        eta_gate=None
    )

    final_result = pass1
    pass2 = None

    # Optionally run new passes, informed with the signal window determined in the first pass
    print("args.number_pass={}".format(args.number_pass))
    for ipass in range(1, args.number_pass):
        pass1 = final_result

        if args.save_pass1_overlay:
            print("save pass1 overlay")
            make_final_plot(
                full_x=pass1["full_x"],
                full_y=pass1["full_y"],
                full_yerr=pass1["full_yerr"],
                bkg_mean=pass1["bkg_mean"],
                fitted_signal=pass1["fitted_signal"],
                full_x_bkg=pass1["full_x_bkg"],
                lo_bkg=pass1["lo_bkg"],
                hi_bkg=pass1["hi_bkg"],
                output_path=args.save_pass1_overlay
            )
            
        pi0_gate_2, eta_gate_2 = build_informed_gate_from_previous_pass(
            pass1,
            k_pi0=args.k_pi0_mask,
            k_eta=args.k_eta_mask
        )
        informed_gate = (pi0_gate_2, eta_gate_2)
        pass2 = run_full_fit(
            h, device=device, dtype=dtype, args=args,
            pi0_gate=pi0_gate_2,
            eta_gate=eta_gate_2
        )
        final_result = pass2

    make_final_plot(
        full_x=final_result["full_x"],
        full_y=final_result["full_y"],
        full_yerr=final_result["full_yerr"],
        bkg_mean=final_result["bkg_mean"],
        fitted_signal=final_result["fitted_signal"],
        full_x_bkg=final_result["full_x_bkg"],
        lo_bkg=final_result["lo_bkg"],
        hi_bkg=final_result["hi_bkg"],
        output_path=args.output
    )

    save_root_graphs(
        full_x=final_result["full_x"],
        bkg_mean=final_result["bkg_mean"],
        fitted_signal=final_result["fitted_signal"],
        graph_template=args.graph_template,
        output_path=args.output_root
    )
    
    if args.verbose:
        print("Saved:", args.output)
        print(f"final pi0 gate: [{final_result['pi0_gate'][0]:.6f}, {final_result['pi0_gate'][1]:.6f}] GeV")
        print(f"final eta gate: [{final_result['eta_gate'][0]:.6f}, {final_result['eta_gate'][1]:.6f}] GeV")
        print("Final CB parameters:", np.array2string(np.asarray(final_result["popt"]), precision=5))
        print("Final peaks:", final_result["peaks"])
        if pass2 is not None:
            y1 = pass1["peaks"]["pi0"]["yield"]; y2 = pass2["peaks"]["pi0"]["yield"]
            print(f"Pass1 pi0 yield={y1:.6g}, Pass2 pi0 yield={y2:.6g}, delta={y2-y1:.6g}, rel={(y2-y1)/y1 if y1!=0 else np.nan:.3%}")
            y1e = pass1["peaks"]["eta"]["yield"]; y2e = pass2["peaks"]["eta"]["yield"]
            print(f"Pass1 eta yield={y1e:.6g}, Pass2 eta yield={y2e:.6g}, delta={y2e-y1e:.6g}, rel={(y2e-y1e)/y1e if y1e!=0 else np.nan:.3%}")
            print(f"Informed gate used: pi0=[{informed_gate[0][1]:.6f},{informed_gate[0][1]:.6f}] GeV, eta=[{informed_gate[1][0]:.6f},{informed_gate[1][1]:.6f}] GeV")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
