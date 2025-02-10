import numpy as np
from scipy.optimize import curve_fit


def lifetime_decay_fit(detuning, x0, x1):
    """Modified ansatz for fitting infidelity curves."""
    return x0 / (x1 + detuning)


def fit_infidelity(detuning_list, infidelity_list):
    """Fits the infidelity data to the modified power-law model with better convergence."""
    p0 = [1, 1]  # Improved initial guess
    params, _ = curve_fit(lifetime_decay_fit, detuning_list, infidelity_list, p0=p0)
    # print(params)
    return params


def _compute_snail_aware_max_dBm(frequency_GHz, f_SNAIL):
    """Returns the maximum allowed pump power (in dBm) for a given pump frequency,
    using a linear (V-shaped) speed limit with a minimum at f0 = f_SNAIL/2.

    Assumptions:
      - At f_left (0.3 GHz): allowed pump power is d_left (25 dBm).
      - At f0 (f_SNAIL/2): allowed pump power is d_min (5 dBm).
    """  # noqa: D205
    d_min = 5  # Minimum allowed dBm at f0
    d_left = 25  # Allowed dBm at f_left
    f_left = f_SNAIL - 4.0  # 0.3  # Left endpoint frequency in GHz
    f0 = f_SNAIL / 2  # (e.g. ~2.138 GHz)
    slope = (d_left - d_min) / (f0 - f_left)
    return d_min + slope * np.abs(frequency_GHz - f0)


def _fit_epsilon(dBm_calib, t_f_calib, g3, lambda_factor, w_pump, w_snail):
    """Calibrates the conversion factor (X_factor) between pump power (in dBm)
    and the effective drive amplitude ε.

    Using:
      (1) (π/2) = 6 g3 |η| λ² t_f_calib
      (2) |η| = (ε * w_snail) / (w_pump² - w_snail²)

    Rearranging, one finds:
      |η|_calib = π / (12 g3 λ² t_f_calib)
      ε_calib = |η|_calib (w_pump² - w_snail²) / w_snail

    Then, assuming a linear relation ε = X_factor * (pump power in dBm):
      X_factor = ε_calib / dBm_calib
    """  # noqa: D205
    eta_calib = np.pi / (12 * t_f_calib * g3 * lambda_factor**2)
    epsilon_calib = eta_calib * (w_pump**2 - w_snail**2) / w_snail
    X_factor = epsilon_calib / dBm_calib
    return epsilon_calib, X_factor


def _compute_gate_duration(frequency_GHz, f_SNAIL, X_factor, g3, lambda_factor):
    """Computes the gate duration t_f (in seconds) for a given pump frequency (in GHz)
    using the calibrated conversion factor (X_factor). Uses the global f_SNAIL.

    The relations used are:
      (1) |η| = (ε * w_snail) / (w_pump² - w_snail²), with ε = X_factor * (pump power in dBm)
      (2) (π/2) = 6 g3 |η| λ² t_f   →   t_f = π / (12 g3 λ² |η|)

    Here, w_pump and w_snail are the pump and SNAIL angular frequencies (rad/s), respectively.
    """  # noqa: D205
    w_pump = 2 * np.pi * frequency_GHz * 1e9
    w_snail = 2 * np.pi * f_SNAIL * 1e9
    max_dBm = _compute_snail_aware_max_dBm(frequency_GHz, f_SNAIL)
    epsilon = X_factor * max_dBm
    eta_val = epsilon * w_snail / (w_pump**2 - w_snail**2)
    t_f = np.pi / (12 * eta_val * g3 * lambda_factor**2)
    return t_f


def speedlimit_infidelity_params(f_SNAIL, t_f_calib, T1, g3, lambdaq):  # noqa: D103
    f0 = f_SNAIL / 2  # (e.g. ~2.138 GHz)
    # f0 to f0+1 GHz => 0 to 2000 MHz detuning
    pump_freq_range = np.linspace(f0, f0 - 2, 100)
    detuning_mhz_list = np.abs((pump_freq_range - f0) * 1e3)

    test_ghz = f_SNAIL / 2 - 1.0  # Calibration pump frequency in GHz

    w_pump_calib = 2 * np.pi * test_ghz * 1e9
    w_snail_calib = 2 * np.pi * f_SNAIL * 1e9
    dBm_calib = _compute_snail_aware_max_dBm(test_ghz, f_SNAIL)

    epsilon_calib, X_factor = _fit_epsilon(
        dBm_calib, t_f_calib, g3, lambdaq, w_pump_calib, w_snail_calib
    )

    # Compute the gate durations for each pump frequency using the calibrated X_factor.
    detuned_durations = np.array(
        [
            _compute_gate_duration(f, f_SNAIL, X_factor, g3, lambdaq)
            for f in pump_freq_range
        ]
    )

    # Estimate the infidelity for each frequency using: infidelity = exp(-t_f/T1)
    fidelity_results = 1 - np.exp(-detuned_durations / T1)

    infidelity_params = fit_infidelity(detuning_mhz_list, fidelity_results)
    # fit_line = lifetime_decay_fit(detuning_mhz_list, *infidelity_params)

    return infidelity_params, fidelity_results
