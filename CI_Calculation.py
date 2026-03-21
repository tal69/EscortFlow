"""
Name:        EscortFlowSim_v5

Purpose:     Take the raw output pickle file of EscortFlowSim_v5.py and process it to create 95% confidence
             intervals for the steady state mean lead time, flow time, and waiting time
             It first apply MSER-5 procedure tp remove warmup period and than optimize batch (block) length
             to minimize variance of batch means while keeping acceptable lag-1 auto correlation

Author:      Tal Raviv,  talraviv@tau.ac.il, feel free to use this code but please contact me and let me know if you do
Create:      11/3/2026 (with some help from ChatGPT)
"""

import argparse
import glob
import pickle
import numpy as np
import mser5

def _lag1_autocorr(z: np.ndarray) -> float:
    """Sample lag-1 autocorrelation of a 1D array z."""
    z = np.asarray(z, dtype=float)
    k = z.size
    if k < 3:
        return np.nan
    x = z[:-1]
    y = z[1:]
    x = x - x.mean()
    y = y - y.mean()
    denom = np.sqrt(np.dot(x, x) * np.dot(y, y))
    if denom > 0:
        return float(np.dot(x, y) / denom)
    # Constant (or near-constant) series: treat as uncorrelated for batching checks.
    return 0.0


def choose_batches_min_var_with_lag1(
    x,
    d: int = 0,
    k_min: int = 20,
    b_min: int = 1,
    b_max: int | None = None,
):
    """
    Choose nonoverlapping batch size b (and implied k = floor(m/b)) that MINIMIZES
    the sample variance of the batch means, subject to:

        |rho1(B)| <= 1.96 / sqrt(k)

    where B are the batch means and k is the number of batches.

    Parameters
    ----------
    x : array-like
        Raw output observations (e.g., per-customer values).
    d : int
        Number of initial observations to drop (warm-up deletion).
    k_min : int
        Minimum number of batches allowed (for stability of CI, etc.). Default 20.
        If you truly don't care, you can set k_min=3.
    b_min : int
        Smallest batch size to consider.
    b_max : int | None
        Largest batch size to consider. If None, set to m // k_min.

    Returns
    -------
    result : dict
        {
          "b": chosen batch size,
          "k": number of batches,
          "rho1": lag-1 autocorr of batch means,
          "threshold": 1.96/sqrt(k),
          "s2": sample variance of batch means (ddof=1),
          "batch_means": array of batch means (length k),
          "used_m": number of post-warmup observations actually used (= k*b)
        }

        If no feasible solution exists, returns None.
    """
    x = np.asarray(x, dtype=float)
    if d < 0 or d >= x.size:
        raise ValueError("d must satisfy 0 <= d < len(x).")

    y = x[d:]
    m = y.size
    if m < 3:
        return None

    if k_min < 3:
        k_min = 3

    if b_max is None:
        b_max = m // k_min
    b_max = int(max(b_min, min(b_max, m // 3)))  # need k>=3 for rho1

    best = None

    for b in range(int(b_min), b_max + 1):
        k = m // b
        if k < k_min:
            continue

        used = k * b
        B = y[:used].reshape(k, b).mean(axis=1)

        rho1 = _lag1_autocorr(B)
        if not np.isfinite(rho1):
            continue

        thr = 1.96 / np.sqrt(k)
        if abs(rho1) > thr:
            continue

        s2 = float(np.var(B, ddof=1))  # sample variance of batch means

        cand = {
            "b": b,
            "k": k,
            "rho1": float(rho1),
            "threshold": float(thr),
            "s2": s2,
            "batch_means": B,
            "used_m": used,
        }

        if best is None:
            best = cand
        else:
            # primary objective: minimize s2
            # tie-breakers: prefer larger k (more batches), then smaller b
            if (cand["s2"] < best["s2"] - 1e-15 or
                (abs(cand["s2"] - best["s2"]) <= 1e-15 and cand["k"] > best["k"]) or
                (abs(cand["s2"] - best["s2"]) <= 1e-15 and cand["k"] == best["k"] and cand["b"] < best["b"])):
                best = cand

    return best


def summarize_metric(x):
    """
    Summarize one metric time series using the current pipeline:
    1) warmup deletion via MSER-5,
    2) batch-size selection with lag-1 constraint,
    3) mean and CI half-widths from batch means.

    Returns a dict with batching diagnostics and:
    - `half_width_95`: 95% CI half-width (normal approximation)
    - `half_width_99`: 99% CI half-width (normal approximation)
    """
    x = np.asarray(x, dtype=float)
    try:
        obs_del = int(mser5.mser5(x)["d_obs"])
        x_post_warmup = x[obs_del:]
        res = choose_batches_min_var_with_lag1(x_post_warmup, d=0, k_min=20)
    except Exception as exc:
        return {
            "obs_deleted": np.nan,
            "obs_used": np.nan,
            "batch_size": np.nan,
            "num_batches": np.nan,
            "lag1_autocorr": np.nan,
            "lag1_threshold": np.nan,
            "batch_means_variance": np.nan,
            "mean": np.nan,
            "half_width_95": np.nan,
            "half_width_99": np.nan,
            "res": None,
            "error": str(exc),
        }

    if res is None:
        return {
            "obs_deleted": obs_del,
            "obs_used": 0,
            "batch_size": np.nan,
            "num_batches": np.nan,
            "lag1_autocorr": np.nan,
            "lag1_threshold": np.nan,
            "batch_means_variance": np.nan,
            "mean": np.nan,
            "half_width_95": np.nan,
            "half_width_99": np.nan,
            "res": None,
            "error": "no feasible batch size found after warmup deletion",
        }

    mean_est = float(np.mean(res["batch_means"]))
    half_width_95 = 1.96 * np.sqrt(res["s2"] / res["k"])
    half_width_99 = 2.576 * np.sqrt(res["s2"] / res["k"])
    return {
        "obs_deleted": obs_del,
        "obs_used": int(res["used_m"]),
        "batch_size": int(res["b"]),
        "num_batches": int(res["k"]),
        "lag1_autocorr": float(res["rho1"]),
        "lag1_threshold": float(res["threshold"]),
        "batch_means_variance": float(res["s2"]),
        "mean": mean_est,
        "half_width_95": float(half_width_95),
        "half_width_99": float(half_width_99),
        "res": res,
        "error": "",
    }


def csv_cell(value):
    """Serialize NaN/None as an empty CSV cell."""
    if value is None:
        return ""
    if isinstance(value, (float, np.floating)) and np.isnan(value):
        return ""
    return str(value)


# ---- example usage ----
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--pickle-file",
        dest="pickle_files",
        nargs="+",
        help="One or more pickle file names to read raw simulation output from, wildcards such as *.p are ok",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--csv",
        dest="csv_file",
        help="CSV file name for summary output",
        default="block_res.csv",
    )
    parser.add_argument(
        "-H",
        "--header",
        action="store_true",
        help="Write CSV header row",
    )

    args = parser.parse_args()
    expanded_pickle_files = []
    for pattern in args.pickle_files:
        if glob.has_magic(pattern):
            matches = sorted(glob.glob(pattern))
            if not matches:
                parser.error(f"Pattern matched no files: {pattern}")
            expanded_pickle_files.extend(matches)
        else:
            expanded_pickle_files.append(pattern)

    # Keep first occurrence order while removing duplicates.
    pickle_files = list(dict.fromkeys(expanded_pickle_files))

    header_cols = [
        "Pickle File", "Algorithm Name", "Queue Management", "Seed", "Request Rate", "PBS Dimensions", "Output Cells",
        "Number of outputs", "Number of Escorts", "Number of Requests", "Simulation End Time",
        "Fractional Horizon", "Integer Horizon", "Execution Horizon", "Time Limit", "Attention Limit",
        "Max Opt Gap", "Actual Attention", "Non Optimal", "Max Gap", "Heuristic Solutions",
        "Lead Time Mean", "Lead Time CI Half Width 95%",
        "Waiting Time Mean", "Waiting Time CI Half Width 95%",
        "Flow Time Mean", "Flow Time CI Half Width 95%",
        "Lead Time Deleted", "Lead Time Used", "Lead Time Batch Size", "Lead Time Number of Batches",
        "Lead Time Lag1 Autocorr", "Lead Time Lag1 Threshold", "Lead Time Batch Means Variance",
        "Waiting Time Deleted", "Waiting Time Used", "Waiting Time Batch Size", "Waiting Time Number of Batches",
        "Waiting Time Lag1 Autocorr", "Waiting Time Lag1 Threshold", "Waiting Time Batch Means Variance",
        "Flow Time Deleted", "Flow Time Used", "Flow Time Batch Size", "Flow Time Number of Batches",
        "Flow Time Lag1 Autocorr", "Flow Time Lag1 Threshold", "Flow Time Batch Means Variance"
    ]
    with open(args.csv_file, "a") as f:
        if args.header:
            f.write(",".join(header_cols) + "\n")
        for pickle_file in pickle_files:
            with open(pickle_file, "rb") as pf:
                raw = pickle.load(pf)

            # Backward-compatible parsing:
            # v5+ raw format includes request_rate right after seed.
            if len(raw) == 22:
                (alg_name, queue_management, seed, request_rate, fractional_horizon, integer_horizon,
                 exec_horizon, time_limit, max_balls_in_air, max_opt_gap, Lx, Ly, O, E_orig, arrivals,
                 departures, start_move, moves, actual_max_balls, non_optimal, max_gap, heuristic_sol) = raw
            elif len(raw) == 21:
                (alg_name, queue_management, seed, fractional_horizon, integer_horizon,
                 exec_horizon, time_limit, max_balls_in_air, max_opt_gap, Lx, Ly, O, E_orig, arrivals,
                 departures, start_move, moves, actual_max_balls, non_optimal, max_gap, heuristic_sol) = raw
                request_rate = np.nan
            else:
                raise ValueError(
                    f"Unsupported raw tuple format in {pickle_file}. Expected length 21 or 22, got {len(raw)}"
                )

            simulation_end_time = int(np.max(departures)) if departures.size else 0
            late_departure_indices = np.flatnonzero(departures > arrivals[-1])
            if late_departure_indices.size > 0:
                cooldown_cutoff_arrival_time = arrivals[late_departure_indices[0]]
                cooldown_trim_mask = arrivals <= cooldown_cutoff_arrival_time
            else:
                cooldown_trim_mask = np.ones_like(arrivals, dtype=bool)

            lead_time = (departures - arrivals)[cooldown_trim_mask]
            waiting_time = (start_move - arrivals)[cooldown_trim_mask]
            flow_time = (departures - start_move)[cooldown_trim_mask]

            lead_stats = summarize_metric(lead_time)
            waiting_stats = summarize_metric(waiting_time)
            flow_stats = summarize_metric(flow_time)

            o_str = str(O).replace('"', '""')
            pbs_dimensions = f"{Lx}x{Ly}"

            row_vals = [
                pickle_file, alg_name, queue_management, seed, request_rate, pbs_dimensions, f'"{o_str}"', len(O), len(E_orig),
                len(arrivals), simulation_end_time,
                fractional_horizon, integer_horizon, exec_horizon, time_limit, max_balls_in_air,
                max_opt_gap, actual_max_balls, non_optimal, max_gap, heuristic_sol,
                lead_stats["mean"], lead_stats["half_width_95"],
                waiting_stats["mean"], waiting_stats["half_width_95"],
                flow_stats["mean"], flow_stats["half_width_95"],
                lead_stats["obs_deleted"], lead_stats["obs_used"], lead_stats["batch_size"], lead_stats["num_batches"],
                lead_stats["lag1_autocorr"], lead_stats["lag1_threshold"], lead_stats["batch_means_variance"],
                waiting_stats["obs_deleted"], waiting_stats["obs_used"], waiting_stats["batch_size"], waiting_stats["num_batches"],
                waiting_stats["lag1_autocorr"], waiting_stats["lag1_threshold"], waiting_stats["batch_means_variance"],
                flow_stats["obs_deleted"], flow_stats["obs_used"], flow_stats["batch_size"], flow_stats["num_batches"],
                flow_stats["lag1_autocorr"], flow_stats["lag1_threshold"], flow_stats["batch_means_variance"],
            ]
            f.write(",".join(csv_cell(v) for v in row_vals) + "\n")

            print(f"\nResults for {pickle_file}")
            for metric_name, stats in [("Lead Time", lead_stats), ("Waiting Time", waiting_stats), ("Flow Time", flow_stats)]:
                if stats["res"] is None:
                    print(f"{metric_name}: {stats['error']}.")
                    continue
                print(
                    f"{metric_name}: b={stats['batch_size']}, k={stats['num_batches']}, "
                    f"rho1={stats['lag1_autocorr']:.4f}, thr={stats['lag1_threshold']:.4f}, "
                    f"s^2={stats['batch_means_variance']:.6g}, mean={stats['mean']:.6g}, "
                    f"95% CI half-width={stats['half_width_95']:.6g}"
                )
