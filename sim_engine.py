"""
Standalone data generation for tempo_subtypes.
No dependency on the bebms_pkg package.

Incorporates:
  - Inlined Mallows-Kendall sampling (from bebms/mallows_kendall.py)
  - dist_type-based very_irregular_distribution (random non-Normal family
    per biomarker per dataset, to stress-test distribution-agnostic
    ordering recovery)
  - flip_directions are fixed at +1 for all biomarkers. Progression direction
    is encoded in the params via sign(theta_mean - phi_mean).
"""

from typing import List, Optional, Dict, Any
import json
import pandas as pd
import numpy as np
import os
from collections import defaultdict
from bisect import bisect_right


# ============================================================
# Mallows-Kendall helpers (inlined from bebms)
# ============================================================

def _theta_to_phi(theta: float) -> float:
    return float(np.exp(-theta))


def _phi_to_theta(phi: float) -> float:
    return float(-np.log(phi))


def _check_theta_phi(theta, phi):
    """Convert between theta/phi; exactly one must be non-None."""
    if not ((phi is None) ^ (theta is None)):
        raise ValueError("Exactly one of theta or phi must be provided.")
    if phi is None and not isinstance(theta, list):
        phi = _theta_to_phi(theta)
    if theta is None and not isinstance(phi, list):
        theta = _phi_to_theta(phi)
    if phi is None and isinstance(theta, list):
        phi = [_theta_to_phi(t) for t in theta]
    if theta is None and isinstance(phi, list):
        theta = [_phi_to_theta(p) for p in phi]
    return np.array(theta), np.array(phi)


def _merge(left, right):
    result, count = [], 0
    i, j, left_len = 0, 0, len(left)
    while i < left_len and j < len(right):
        if left[i] <= right[j]:
            result.append(left[i]); i += 1
        else:
            result.append(right[j]); count += left_len - i; j += 1
    result += left[i:]
    result += right[j:]
    return result, count


def _merge_sort_dist(lst):
    lst = list(lst)
    if len(lst) <= 1:
        return lst, 0
    mid = len(lst) // 2
    left, a  = _merge_sort_dist(lst[:mid])
    right, b = _merge_sort_dist(lst[mid:])
    sorted_, c = _merge(left, right)
    return sorted_, a + b + c


def _v_to_ranking(v, n):
    """Decomposition vector → permutation (0-indexed ranks)."""
    rem = list(range(n))
    rank = np.full(n, np.nan)
    for i in range(len(v)):
        rank[i] = rem[v[i]]
        rem.pop(v[i])
    return rank.astype(int)

"""Original authors: https://github.com/ekhiru/top-k-mallows
"""
def mallows_sample(
    m: int,
    n: int,
    *,
    theta: Optional[float] = None,
    phi: Optional[float] = None,
    s0: Optional[np.ndarray] = None,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Sample m complete rankings of n items from a Mallows model.

    Parameters
    ----------
    m    : number of rankings to generate
    n    : length of each ranking
    theta: dispersion (provide theta OR phi, not both)
    phi  : dispersion
    s0   : consensus ranking (default: identity 0..n-1)
    rng  : numpy Generator for reproducibility

    Returns
    -------
    ndarray of shape (m, n)  — each row is a rank array (sigma[i] = rank of item i)
    """
    theta, phi = _check_theta_phi(theta, phi)
    theta = np.full(n - 1, float(theta))

    if s0 is None:
        s0 = np.arange(n)

    rnge = np.arange(n - 1)
    # Normalisation constants psi[j] = sum_{r=0}^{n-j-1} exp(-theta[j]*r)
    psi = (1.0 - np.exp((-n + rnge) * theta)) / (1.0 - np.exp(-theta))

    # Probability tables for decomposition-vector components
    vprobs = np.zeros((n, n))
    for j in range(n - 1):
        vprobs[j, 0] = 1.0 / psi[j]
        for r in range(1, n - j):
            vprobs[j, r] = np.exp(-theta[j] * r) / psi[j]
        # Normalise row to guard against floating-point drift
        row_sum = vprobs[j, :].sum()
        if row_sum > 0:
            vprobs[j, :] /= row_sum

    samples = []
    for _ in range(m):
        v = [int(rng.choice(n, p=vprobs[i, :])) for i in range(n - 1)]
        v += [0]
        ranking = _v_to_ranking(v, n)
        samples.append(ranking)

    # Apply consensus permutation
    return np.array([s[s0] for s in samples])


# ============================================================
# Utility functions
# ============================================================

def dirichlet_multinomial(
    subtype_assignment_prior: int,
    total_participant: int,
    n_subtypes: int,
    rng: np.random.Generator,
) -> np.ndarray:
    alpha = np.ones(n_subtypes) * subtype_assignment_prior
    if n_subtypes > total_participant:
        raise ValueError(
            f"Cannot assign at least one item to {n_subtypes} subtypes "
            f"with only {total_participant} participants."
        )
    base_counts = np.ones(n_subtypes, dtype=int)
    remaining_items = total_participant - n_subtypes
    p = rng.dirichlet(alpha)
    additional_counts = rng.multinomial(remaining_items, p)
    return base_counts + additional_counts


def kendalls_w(rank_matrix: np.ndarray) -> float:
    """Kendall's W (coefficient of concordance) for complete rankings."""
    n_raters, n_items = rank_matrix.shape
    R = np.sum(rank_matrix, axis=0)
    R_bar = np.mean(R)
    SS = np.sum((R - R_bar) ** 2)
    return 12 * SS / (n_raters ** 2 * (n_items ** 3 - n_items))


# ============================================================
# Data-generation helpers
# ============================================================

def get_rank(sorted_et: np.ndarray, val: float) -> int:
    return int(bisect_right(sorted_et, val))


def very_irregular_distribution(
    dist_type: int,
    bm_params: Dict[str, float],
    state: str = "affected",
    size: int = 100_000,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """
    Non-normal samples for one biomarker/state.

    dist_type (0-5) is randomly assigned per biomarker per dataset at the call
    site, making the function agnostic to biomarker identity — suitable for
    arbitrarily many biomarkers (low-dim and high-dim alike).

        0: Triangular + Normal + Exponential
        1: Pareto + Uniform + Logistic
        2: Beta + Exponential(signed) + Normal(spike)
        3: Gamma + Weibull + Normal(±σ)
        4: Heavy-tailed Cauchy (clipped)
        5: Bimodal Normal-spike + Logistic
    """
    if rng is None:
        rng = np.random.default_rng()

    mean = bm_params["theta_mean"] if state == "affected" else bm_params["phi_mean"]
    std  = bm_params["theta_std"]  if state == "affected" else bm_params["phi_std"]

    base = np.zeros(size)
    seg1, seg2, seg3 = np.array_split(np.arange(size), 3)

    if dist_type == 0:
        base[seg1] = rng.triangular(mean - 2*std, mean - 1.5*std, mean, size=len(seg1))
        base[seg2] = rng.normal(mean + std, 0.3*std, size=len(seg2))
        base[seg3] = rng.exponential(scale=0.7*std, size=len(seg3)) + mean - 0.5*std
    elif dist_type == 1:
        base[seg1] = rng.pareto(1.5, size=len(seg1)) * std + mean - 2*std
        base[seg2] = rng.uniform(mean - 1.5*std, mean + 1.5*std, size=len(seg2))
        base[seg3] = rng.logistic(loc=mean, scale=std, size=len(seg3))
    elif dist_type == 2:
        base[seg1] = rng.beta(0.5, 0.5, size=len(seg1)) * 4*std + mean - 2*std
        base[seg2] = rng.exponential(scale=std*0.4, size=len(seg2)) * rng.choice([-1, 1], size=len(seg2)) + mean
        base[seg3] = rng.normal(mean, std*0.5, size=len(seg3)) + rng.choice([0, std*2], size=len(seg3))
    elif dist_type == 3:
        base[seg1] = rng.gamma(shape=2, scale=0.5*std, size=len(seg1)) + mean - std
        base[seg2] = rng.weibull(1.0, size=len(seg2)) * std + mean - std
        base[seg3] = rng.normal(mean, std*0.5, size=len(seg3)) + rng.choice([-1, 1], size=len(seg3)) * std
    elif dist_type == 4:
        raw = rng.standard_cauchy(size=size) * std + mean
        raw += rng.normal(0, 0.2*std, size=size)
        base = np.clip(raw, mean - 4*std, mean + 4*std)
    else:  # dist_type == 5
        spike = size // 10
        base[:spike] = rng.normal(mean, 0.2*std, size=spike)
        base[spike:]  = rng.logistic(loc=mean + std, scale=2*std, size=size - spike)

    base += rng.normal(0, 0.2*std, size=size)
    base  = np.clip(base, mean - 5*std, mean + 5*std)
    return base


def generate_measurements_ebm(
    params: Dict[str, Dict[str, float]],
    event_time_dict: Dict[str, float],
    shuffled_biomarkers: np.ndarray,
    experiment_name: str,
    all_kjs: np.ndarray,
    all_diseased: np.ndarray,
    keep_all_cols: bool,
    dist_type_dict: Dict[str, int],
    rng: Optional[np.random.Generator] = None,
) -> List[Dict[str, Any]]:
    """
    EBM binary-switch measurements.

    dist_type_dict is pre-sampled once per dataset (before the subtype loop)
    so that all subtypes share the same distribution family per biomarker.
    """
    if rng is None:
        rng = np.random.default_rng()

    data = []
    irreg_dict: Dict = defaultdict(dict)
    if "xnjNonNormal" in experiment_name:
        for biomarker in shuffled_biomarkers:
            bm_params = params[biomarker]
            for state in ("affected", "nonaffected"):
                irreg_dict[biomarker][state] = very_irregular_distribution(
                    dist_type_dict[biomarker], bm_params, state=state, rng=rng
                )

    for participant_id, (disease_stage, is_diseased) in enumerate(
        zip(all_kjs, all_diseased)
    ):
        for biomarker in shuffled_biomarkers:
            bm_params  = params[biomarker]
            event_time = event_time_dict[biomarker]

            if is_diseased and disease_stage >= event_time:
                state             = "affected"
                distribution_param = "theta"
            else:
                state             = "nonaffected"
                distribution_param = "phi"

            if "xnjNormal" in experiment_name:
                measurement = rng.normal(
                    bm_params[f"{distribution_param}_mean"],
                    bm_params[f"{distribution_param}_std"],
                )
            else:
                measurement = rng.choice(irreg_dict[biomarker][state], size=1)[0]

            record: Dict[str, Any] = {
                "participant": participant_id,
                "biomarker":   biomarker,
                "measurement": measurement,
                "diseased":    is_diseased,
            }
            if keep_all_cols:
                record.update({
                    "event_time": event_time,
                    "k_j":        disease_stage,
                    "affected":   disease_stage >= event_time,
                })
            data.append(record)
    return data


def generate_measurements_sigmoid(
    experiment_name: str,
    event_time_dict: Dict[str, float],
    all_kjs: np.ndarray,
    all_diseased: np.ndarray,
    shuffled_biomarkers: np.ndarray,
    params: Dict[str, Dict[str, float]],
    keep_all_cols: bool,
    noise_std_parameter: float,
    flip_directions: Dict[str, int],
    rng: Optional[np.random.Generator] = None,
) -> List[Dict[str, Any]]:
    """
    Sigmoid-progression measurements.

    flip_directions is pre-sampled once per dataset (before the subtype loop)
    so that all subtypes share the same direction per biomarker.  This ensures
    subtypes are separated by their orderings, not trivially by direction.
    """
    if rng is None:
        rng = np.random.default_rng()

    # Precompute R and rho
    for biomarker in params:
        theta_m = params[biomarker]["theta_mean"]
        phi_m   = params[biomarker]["phi_mean"]
        theta_v = params[biomarker]["theta_std"] ** 2
        phi_v   = params[biomarker]["phi_std"]   ** 2
        params[biomarker]["R"]   = theta_m - phi_m
        params[biomarker]["rho"] = max(1, abs(theta_m - phi_m) / np.sqrt(theta_v + phi_v))

    data = []
    max_stage = len(shuffled_biomarkers)

    for participant_id, (disease_stage, is_diseased) in enumerate(
        zip(all_kjs, all_diseased)
    ):
        if experiment_name.startswith("xiNearNormalWithNoise"):
            noise_std = max_stage * noise_std_parameter
            noises    = rng.normal(loc=0, scale=noise_std, size=max_stage)

        for biomarker_idx, biomarker in enumerate(shuffled_biomarkers):
            event_time = event_time_dict[biomarker]
            if experiment_name.startswith("xiNearNormalWithNoise"):
                event_time = float(np.clip(event_time + noises[biomarker_idx], 0, max_stage))

            bm_params           = params[biomarker]
            healthy_measurement = rng.normal(bm_params["phi_mean"], bm_params["phi_std"])

            if is_diseased:
                prog_magnitude = flip_directions[biomarker] * bm_params["R"]
                prog_rate      = bm_params["rho"]
                sigmoid_term   = prog_magnitude / (1 + np.exp(-prog_rate * (disease_stage - event_time)))
                measurement    = sigmoid_term + healthy_measurement
            else:
                measurement = healthy_measurement

            record: Dict[str, Any] = {
                "participant": participant_id,
                "biomarker":   biomarker,
                "measurement": measurement,
                "diseased":    is_diseased,
            }
            if keep_all_cols:
                record.update({
                    "event_time": event_time,
                    "k_j":        disease_stage,
                    "affected":   disease_stage >= event_time,
                })
            data.append(record)
    return data


def generate_data(
    filename: str,
    experiment_name: str,
    params: Dict[str, Dict[str, float]],
    n_participants: int,
    healthy_ratio: float,
    output_dir: str,
    m: int,
    dirichlet_alpha: Dict[str, List[float]],
    beta_params: Dict[str, Dict[str, float]],
    prefix: Optional[str],
    suffix: Optional[str],
    keep_all_cols: bool,
    fixed_biomarker_order: bool,
    noise_std_parameter: float,
    true_order_and_stages_dict: Dict,
    rng: np.random.Generator,
    save2file: bool,
    flip_directions: Dict[str, int],
    dist_type_dict: Dict[str, int],
) -> pd.DataFrame:
    """
    Generate one (sub)dataset for a given experiment configuration.

    flip_directions and dist_type_dict are passed in (pre-sampled per dataset)
    to ensure consistency across subtypes within the same dataset.
    """
    assert n_participants > 0
    assert 0 <= healthy_ratio <= 1

    if fixed_biomarker_order:
        shuffled_biomarkers = np.array(list(params.keys()))
    else:
        shuffled_biomarkers = rng.permutation(np.array(list(params.keys())))

    max_stage   = len(shuffled_biomarkers)
    event_times = np.arange(1, max_stage + 1)

    n_healthy  = int(n_participants * healthy_ratio)
    n_diseased = n_participants - n_healthy

    # ----------------------------------------------------------------
    if "kjOrdinal" in experiment_name:
        event_time_dict = dict(zip(shuffled_biomarkers, event_times))

        if "Uniform" in experiment_name:
            if len(dirichlet_alpha["uniform"]) != max_stage:
                dirichlet_alphas = [dirichlet_alpha["uniform"][0]] * max_stage
            else:
                dirichlet_alphas = dirichlet_alpha["uniform"]
            stage_probs = rng.dirichlet(dirichlet_alphas)
        else:
            stage_probs = rng.dirichlet(dirichlet_alpha["multinomial"][:max_stage])

        stage_counts   = rng.multinomial(n_diseased, stage_probs)
        disease_stages = np.repeat(np.arange(1, max_stage + 1), stage_counts)
        all_kjs        = np.concatenate([np.zeros(n_healthy), disease_stages])
        all_diseased   = all_kjs > 0

        shuffle_idx  = rng.permutation(n_participants)
        all_kjs      = all_kjs[shuffle_idx]
        all_diseased = all_diseased[shuffle_idx]

        data = generate_measurements_ebm(
            params, event_time_dict, shuffled_biomarkers, experiment_name,
            all_kjs, all_diseased, keep_all_cols,
            dist_type_dict=dist_type_dict, rng=rng,
        )
        true_stages = [int(x) for x in all_kjs]
        true_stages_continuous = []

    # ----------------------------------------------------------------
    else:
        epsilon = 1e-8

        if experiment_name.startswith("xi"):
            if fixed_biomarker_order:
                event_time_raw = np.sort(rng.beta(
                    a=beta_params["near_normal"]["alpha"],
                    b=beta_params["near_normal"]["beta"],
                    size=max_stage,
                ))
            else:
                event_time_raw = rng.beta(
                    a=beta_params["near_normal"]["alpha"],
                    b=beta_params["near_normal"]["beta"],
                    size=max_stage,
                )
            event_times = event_time_raw * max_stage + epsilon

        event_time_dict = dict(zip(shuffled_biomarkers, event_times))

        if "kjContinuousUniform" in experiment_name:
            disease_stages_raw = rng.beta(
                a=beta_params["uniform"]["alpha"],
                b=beta_params["uniform"]["beta"],
                size=n_diseased,
            ) + epsilon
        else:
            disease_stages_raw = rng.beta(
                a=beta_params["regular"]["alpha"],
                b=beta_params["regular"]["beta"],
                size=n_diseased,
            ) + epsilon

        disease_stages = disease_stages_raw * max_stage
        all_kjs        = np.concatenate([np.zeros(n_healthy), disease_stages])
        all_diseased   = all_kjs > 0

        shuffle_idx  = rng.permutation(n_participants)
        all_kjs      = all_kjs[shuffle_idx]
        all_diseased = all_diseased[shuffle_idx]

        if "sigmoid" in experiment_name:
            data = generate_measurements_sigmoid(
                experiment_name, event_time_dict, all_kjs, all_diseased,
                shuffled_biomarkers, params, keep_all_cols,
                noise_std_parameter=noise_std_parameter,
                flip_directions=flip_directions, rng=rng,
            )
        else:
            data = generate_measurements_ebm(
                params, event_time_dict, shuffled_biomarkers, experiment_name,
                all_kjs, all_diseased, keep_all_cols,
                dist_type_dict=dist_type_dict, rng=rng,
            )

        sorted_event_times = sorted(event_times)
        true_stages = [get_rank(sorted_event_times, x) for x in all_kjs]
        true_stages_continuous = (
            [float(x) for x in all_kjs]
            if "kjContinuous" in experiment_name
            else []
        )

    # ----------------------------------------------------------------
    df = pd.DataFrame(data)
    if save2file:
        df.to_csv(os.path.join(output_dir, f"{filename}.csv"), index=False)

    true_order_and_stages_dict[filename]["true_order"] = dict(zip(
        sorted(event_time_dict, key=lambda x: event_time_dict[x]),
        range(1, max_stage + 1),
    ))
    true_order_and_stages_dict[filename]["true_stages"] = true_stages.copy()
    true_order_and_stages_dict[filename]["true_stages_continuous"] = (
        true_stages_continuous.copy()
    )
    if experiment_name.startswith("xi"):
        true_order_and_stages_dict[filename]["true_order_continuous"] = (
            event_time_dict.copy()
        )
    else:
        true_order_and_stages_dict[filename]["true_order_continuous"] = {}
    return df


def dirichlet_near_normal(
    n_biomarkers: int, peak_height: float = 4.25, min_height: float = 0.35
) -> list:
    x      = np.arange(n_biomarkers)
    center = (n_biomarkers - 1) / 2
    sigma  = n_biomarkers / 6
    curve  = np.exp(-0.5 * ((x - center) / sigma) ** 2)
    curve  = (curve - curve.min()) / (curve.max() - curve.min())
    curve  = curve * (peak_height - min_height) + min_height
    return curve.tolist()


def generate(
    experiment_name: str = "sn_kjOrdinalDM_xnjNormal",
    params: Optional[Dict] = None,
    params_file: str = "params.json",
    js: List[int] = [50, 200, 500, 1000],
    rs: List[float] = [0.1, 0.25, 0.5, 0.75, 0.9],
    num_of_datasets_per_combination: int = 50,
    output_dir: str = "data",
    seed: int = 53,
    dirichlet_alpha: Optional[Dict] = None,
    beta_params: Optional[Dict] = None,
    prefix: Optional[str] = None,
    suffix: Optional[str] = None,
    keep_all_cols: bool = False,
    fixed_biomarker_order: bool = True,
    noise_std_parameter: float = 0.05,
    temperature_lo: float = 0.1,
    temperature_hi: float = 1.0,
    n_sub_lo: int = 2,
    n_sub_hi: int = 5,
    save2file: bool = False,
    provided_subtype_orders: Optional[np.ndarray] = None,
    subtype_dirichlet_priors: Optional[List[int]] = None,
    subtype_length_lo: int = 5,
) -> Dict:
    """
    Generate multiple datasets for all (j, r, variant) combinations.

    Returns a dict mapping filename → metadata dict (true orderings, stages, etc.).
    """
    if dirichlet_alpha is None:
        dirichlet_alpha = {
            "uniform":     [100],
            "multinomial": [0.35, 0.85, 1.55, 2.45, 3.45, 4.25, 4.25,
                            3.45, 2.45, 1.55, 0.85, 0.35],
        }
    if beta_params is None:
        beta_params = {
            "near_normal": {"alpha": 2.0, "beta": 2.0},
            "uniform":     {"alpha": 1,   "beta": 1},
            "regular":     {"alpha": 5,   "beta": 2},
        }

    os.makedirs(output_dir, exist_ok=True)

    def get_total(r: float, n_diseased: int) -> int:
        return round(n_diseased / (1 - r))

    if params is None:
        with open(params_file) as f:
            params = json.load(f)

    biomarker_names = sorted(params.keys())
    rng             = np.random.default_rng(seed)
    true_order_and_stages_dict = defaultdict(dict)

    # ------------------------------------------------------------
    # GLOBAL flip directions.
    # Progression direction is already encoded by sign(theta_mean - phi_mean),
    # so every biomarker uses flip=+1 in both lowdim and highdim. This keeps
    # highdim aligned with the lowdim/ADNI convention instead of adding an
    # extra random direction flip.
    #
    # NOTE: flip_directions is only used by the sigmoid measurement model
    # (generate_measurements_sigmoid), which is invoked for "sigmoid"
    # experiments. The EBM measurement model (kjOrdinal experiments)
    # does not use flip_directions — direction is implicit in which
    # distribution (theta vs phi) is sampled from.
    # ------------------------------------------------------------
    flip_directions_global = {b: 1 for b in biomarker_names}

    for participant_count in js:
        for healthy_ratio in rs:
            for variant in range(num_of_datasets_per_combination):
                sub_seed = int(rng.integers(0, 1_000_000))
                sub_rng  = np.random.default_rng(sub_seed)

                if num_of_datasets_per_combination >= 2:
                    filename = f"j{participant_count}_r{healthy_ratio}_E{experiment_name}_m{variant}"
                else:
                    filename = f"j{participant_count}_r{healthy_ratio}_E{experiment_name}"
                if prefix:
                    filename = f"{prefix}_{filename}"
                if suffix:
                    filename = f"{filename}_{suffix}"

                if len(params) != len(dirichlet_alpha["multinomial"]):
                    dirichlet_alpha["multinomial"] = dirichlet_near_normal(
                        n_biomarkers=len(params)
                    )

                n_healthy  = participant_count * healthy_ratio
                n_diseased = int(participant_count - n_healthy)

                # ----------------------------------------------------------
                # Pre-sample per-dataset constants shared across all subtypes
                # ----------------------------------------------------------
                # Flip direction: GLOBAL, sampled once before the dataset
                # loop. See block above for lowdim/highdim regimes.
                flip_directions = flip_directions_global

                # Non-normal family: one integer (0-5) per biomarker per dataset.
                dist_type_dict = {
                    b: int(sub_rng.integers(0, 6))
                    for b in biomarker_names
                }

                # ----------------------------------------------------------
                # Subtype orderings
                # ----------------------------------------------------------
                if provided_subtype_orders is None:
                    N_SUB       = int(sub_rng.integers(n_sub_lo, n_sub_hi + 1))
                    TEMPERATURE = float(sub_rng.uniform(temperature_lo, temperature_hi))

                    s0 = sub_rng.permutation(np.arange(len(params)))
                    len_unique_sampled = 0
                    while len_unique_sampled != N_SUB:
                        mk_sampled = mallows_sample(
                            m=N_SUB, n=len(params),
                            theta=TEMPERATURE, s0=s0, rng=sub_rng,
                        )
                        len_unique_sampled = len(set(tuple(arr) for arr in mk_sampled))
                    SUBTYPE_RANKINGS = mk_sampled
                else:
                    SUBTYPE_RANKINGS = provided_subtype_orders
                    N_SUB            = len(SUBTYPE_RANKINGS)
                    TEMPERATURE      = float("nan")

                W = kendalls_w(SUBTYPE_RANKINGS)

                # ----------------------------------------------------------
                # Subtype sizes (Dirichlet-Multinomial)
                # ----------------------------------------------------------
                SUBTYPE_LENGTHS = np.zeros(N_SUB, dtype=np.int64)
                max_attempt = 50
                success     = False
                for _ in range(max_attempt + 1):
                    alpha_val = int(rng.choice(subtype_dirichlet_priors))
                    SUBTYPE_LENGTHS = dirichlet_multinomial(
                        subtype_assignment_prior=alpha_val,
                        total_participant=n_diseased,
                        n_subtypes=N_SUB,
                        rng=sub_rng,
                    )
                    if np.all(SUBTYPE_LENGTHS >= subtype_length_lo):
                        success = True
                        break

                if not success:
                    raise ValueError(
                        f"Failed to ensure all subtypes have ≥ {subtype_length_lo} participants!"
                    )

                true_order_and_stages_dict[filename]["N_SUB"]        = int(N_SUB)
                true_order_and_stages_dict[filename]["TEMPERATURE"]  = TEMPERATURE
                true_order_and_stages_dict[filename]["TRUE_ORDERINGS"] = SUBTYPE_RANKINGS
                true_orderings_continuous = []
                true_order_and_stages_dict[filename]["CONCENTRATION"] = float(W)

                # ----------------------------------------------------------
                # Generate data per subtype then assemble
                # ----------------------------------------------------------
                FULL_DF               = []
                new_participant_start = 0

                for subtype_idx in range(N_SUB):
                    subtype_p_count  = get_total(
                        r=healthy_ratio, n_diseased=int(SUBTYPE_LENGTHS[subtype_idx])
                    )
                    subtype_ranking = SUBTYPE_RANKINGS[subtype_idx]

                    # Build params_use in event-order for this subtype
                    params_use = {}
                    for sorted_idx in np.argsort(subtype_ranking):
                        bm_str = biomarker_names[sorted_idx]
                        params_use[bm_str] = params[bm_str]

                    subtype_dict: Dict = defaultdict(dict)

                    df = generate_data(
                        filename=filename,
                        experiment_name=experiment_name,
                        params=params_use,
                        n_participants=subtype_p_count,
                        healthy_ratio=healthy_ratio,
                        output_dir=output_dir,
                        m=variant,
                        dirichlet_alpha=dirichlet_alpha,
                        beta_params=beta_params,
                        prefix=prefix,
                        suffix=suffix,
                        keep_all_cols=keep_all_cols,
                        fixed_biomarker_order=True,
                        noise_std_parameter=noise_std_parameter,
                        true_order_and_stages_dict=subtype_dict,
                        rng=sub_rng,
                        save2file=False,          # written once after assembly
                        flip_directions=flip_directions,
                        dist_type_dict=dist_type_dict,
                    )
                    order_continuous = subtype_dict[filename].get(
                        "true_order_continuous", {}
                    )
                    if order_continuous:
                        true_orderings_continuous.append([
                            float(order_continuous[b]) for b in biomarker_names
                        ])

                    if len(set(df["diseased"])) != 2:
                        raise ValueError("Zero-length subtype cluster!")

                    if not keep_all_cols:
                        diseased_dict_old = dict(zip(df.participant, df.diseased))
                        dff = df.pivot(
                            index="participant", columns="biomarker", values="measurement"
                        )
                        dff = dff.reindex(columns=biomarker_names, level=1)
                        dff.columns.name = None
                        dff.reset_index(inplace=True, drop=False)
                        dff["diseased"] = dff["participant"].map(diseased_dict_old)
                        dff.sort_values(by="participant", inplace=True)
                    else:
                        dff = df.copy()
                        dff.sort_values(by="participant", inplace=True)

                    old_unique  = pd.unique(df["participant"])
                    stage_map   = dict(zip(df["participant"].unique(),
                                          subtype_dict[filename]["true_stages"]))
                    stage_continuous = subtype_dict[filename].get(
                        "true_stages_continuous", []
                    )
                    stage_continuous_map = dict(
                        zip(df["participant"].unique(), stage_continuous)
                    ) if len(stage_continuous) > 0 else {}
                    dff["stage_assignments"]   = dff["participant"].map(stage_map)
                    if stage_continuous_map:
                        dff["stage_assignments_continuous"] = (
                            dff["participant"].map(stage_continuous_map)
                        )
                    dff["subtype_assignments"] = subtype_idx
                    new_ids    = np.arange(new_participant_start,
                                           new_participant_start + len(dff))
                    old_to_new = dict(zip(old_unique, new_ids))
                    dff["participant"] = dff["participant"].map(old_to_new)

                    FULL_DF.append(dff)
                    new_participant_start += subtype_p_count

                full_data = (
                    pd.concat(FULL_DF, ignore_index=True)
                    .sort_values(by="participant")
                )
                subtype_assignments = [
                    int(subtype) if bool(diseased) else None
                    for subtype, diseased in zip(
                        full_data["subtype_assignments"],
                        full_data["diseased"],
                    )
                ]
                true_order_and_stages_dict[filename]["TRUE_SUBTYPE_ASSIGNMENTS"] = (
                    subtype_assignments
                )
                if true_orderings_continuous:
                    true_order_and_stages_dict[filename][
                        "TRUE_ORDERINGS_CONTINUOUS"
                    ] = true_orderings_continuous
                else:
                    true_order_and_stages_dict[filename][
                        "TRUE_ORDERINGS_CONTINUOUS"
                    ] = []
                true_order_and_stages_dict[filename]["TRUE_STAGE_ASSIGNMENTS"] = list(
                    full_data["stage_assignments"]
                )
                if "stage_assignments_continuous" in full_data.columns:
                    true_order_and_stages_dict[filename][
                        "TRUE_STAGE_ASSIGNMENTS_CONTINUOUS"
                    ] = list(full_data["stage_assignments_continuous"])
                else:
                    true_order_and_stages_dict[filename][
                        "TRUE_STAGE_ASSIGNMENTS_CONTINUOUS"
                    ] = []
                true_order_and_stages_dict[filename]["DISEASED_ARR"] = [
                    int(bool(x)) for x in full_data["diseased"]
                ]
                drop_columns = ["subtype_assignments", "stage_assignments"]
                if "stage_assignments_continuous" in full_data.columns:
                    drop_columns.append("stage_assignments_continuous")
                full_data.drop(
                    columns=drop_columns, inplace=True
                )
                full_data.to_csv(f"{output_dir}/{filename}.csv", index=False)

    print(f"Data generation complete. Files saved in {output_dir}/")
    return true_order_and_stages_dict
