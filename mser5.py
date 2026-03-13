import numpy as np

def mser5(x, batch_size=5, min_remaining_batches=10):
    """
    MSER-5 warm-up estimator (Mean Squared Error Reduction with batch size 5).

    See:
    White, K. P., Jr. (1997). An effective truncation heuristic for bias reduction in simulation output. Simulation, 69(6), 323–334.
    Spratt, S. C. (1998). Heuristics for the startup problem. M.S. Thesis, Department of Systems Engineering, University of Virginia.
    White, K. P., Jr., Cobb, M. J., & Spratt, S. C. (2000). A comparison of five steady-state truncation heuristics for simulation. In Proceedings of the 2000 Winter Simulation Conference.
    Franklin, W. W., & White, K. P., Jr. (2008). Stationarity tests and MSER-5… In Proceedings of the 2008 Winter Simulation Conference.


    Steps:
      1) Partition x into nonoverlapping batches of size 5 (or batch_size).
      2) Replace each batch by its mean -> y_1,...,y_N (batch means).
      3) For each candidate deletion d (in batches), compute
            MSER(d) = s^2_d / n_d
         where s^2_d is the sample variance (ddof=1) of y_{d+1...N}
         and n_d is the number of remaining batches.
      4) Choose d minimizing MSER(d).
      5) Return deletion in *original observations* = d * batch_size.

    Parameters
    ----------
    x : array-like
        Sequential output observations (e.g., per-customer values).
    batch_size : int
        Use 5 for MSER-5; kept as a parameter for convenience.
    min_remaining_batches : int
        Small guardrail: require at least this many batches after deletion.

    Returns
    -------
    result : dict
        {
          "d_obs": deletion point in original observations,
          "d_batches": deletion point in batches,
          "batch_size": batch_size,
          "scores": MSER score array for each d_batches (NaN where invalid),
          "batch_means": y array of batch means,
          "N_batches": number of batch means used,
        }

        If not enough data to compute, raises ValueError.
    """
    x = np.asarray(x, dtype=float)
    if batch_size <= 0:
        raise ValueError("batch_size must be positive.")
    if x.size < 2 * batch_size:
        raise ValueError("Not enough data: need at least 2 batches to start.")

    # Form batch means
    N = x.size // batch_size
    if N < 2:
        raise ValueError("Not enough complete batches.")
    xb = x[:N * batch_size].reshape(N, batch_size)
    y = xb.mean(axis=1)  # batch means, length N

    if N - 0 < max(2, min_remaining_batches):
        raise ValueError(
            f"Not enough batches ({N}) to leave min_remaining_batches={min_remaining_batches}."
        )

    scores = np.full(N, np.nan, dtype=float)

    # Candidate deletions in batches: d = 0..N-min_remaining_batches
    max_d = N - min_remaining_batches
    for d in range(0, max_d + 1):
        z = y[d:]
        n = z.size
        if n < 2:
            continue
        s2 = np.var(z, ddof=1)
        scores[d] = s2 / n  # estimated MSE of the sample mean of remaining batches

    if not np.isfinite(scores).any():
        raise ValueError("No valid MSER scores computed (try smaller min_remaining_batches).")

    d_star = int(np.nanargmin(scores))
    return {
        "d_obs": d_star * batch_size,
        "d_batches": d_star,
        "batch_size": batch_size,
        "scores": scores,
        "batch_means": y,
        "N_batches": N,
    }


# ----------------- Example -----------------
if __name__ == "__main__":
    rng = np.random.default_rng(0)
    x = np.r_[np.linspace(5, 0, 3000) + rng.normal(0, 1, 3000),  # transient drift
              rng.normal(0, 1, 20000)]                           # steady-ish part

    res = mser5(x, batch_size=5, min_remaining_batches=30)
    print("Drop (observations):", res["d_obs"])
    print("Drop (batches):", res["d_batches"])
    print("Total batches:", res["N_batches"])