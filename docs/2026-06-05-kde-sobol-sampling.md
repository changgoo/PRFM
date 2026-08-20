# KDE-Sobol Sampling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the LHS candidate-pool approach in `prfm/phangs_sampling.py` with KDE-marginal-quantile Sobol sampling over the 5 design fields only (no f_mol, no SFR), and deliver two standalone executable scripts for running the sampling and generating diagnostic figures.

**Architecture:** Add `synthesize_kde_sobol()` and `synthesize_expanded_kde_sobol()` to `PHANGSSamplingDesigner`; the KDE is fit on `log10(design_fields)` only (5D), marginal quantile functions are extracted from a large auxiliary KDE draw, and a scrambled Sobol sequence is mapped through those quantile functions. A `qshear_max=1.5` hard bound rejects super-Keplerian points with sequential replacement that preserves progressive nesting. Two scripts in `project/scripts/` use this API end-to-end.

**Tech Stack:** `scipy.stats.qmc.Sobol`, `scipy.stats.gaussian_kde`, `numpy`, `astropy`, `matplotlib`, existing `prfm.phangs` and `prfm.phangs_sampling` modules.

---

## File Map

| Action | Path | Responsibility |
|--------|------|---------------|
| Modify | `prfm/phangs_sampling.py` | Add `SynthesisMethod` literals, `SamplingConfig` fields, `_compute_marginal_quantiles()`, `synthesize_kde_sobol()`, `synthesize_expanded_kde_sobol()`, update `sample()` dispatch |
| Modify | `tests/test_phangs.py` | Unit tests for new Sobol methods using synthetic data (no real PHANGS data) |
| Create | `project/scripts/run_sampling.py` | CLI: load data → select reference → KDE-Sobol at n=64/128/256 → fairness report → save CSVs |
| Create | `project/scripts/plot_sampling_diagnostics.py` | CLI: load data → sample → save correlation matrix, distribution overlay, pair overlay figures |
| Create | `project/output/.gitkeep` | Empty dir for script outputs (tracked in git) |

---

## Task 1: Update `SamplingConfig` and `SynthesisMethod`

**Files:**
- Modify: `prfm/phangs_sampling.py` (top of file, `SamplingConfig` dataclass)
- Test: `tests/test_phangs.py`

- [ ] **Step 1: Write failing tests for new config fields**

Add this class at the bottom of `tests/test_phangs.py`:

```python
class TestSamplingConfig:
    def test_default_synthesis_method_is_kde_sobol(self):
        from prfm.phangs_sampling import SamplingConfig
        cfg = SamplingConfig()
        assert cfg.synthesis_method == "kde_sobol"

    def test_design_fields_are_five(self):
        from prfm.phangs_sampling import SamplingConfig
        cfg = SamplingConfig()
        assert cfg.design_fields == [
            "Sigma_gas", "Sigma_star", "H_star", "Omega", "qshear"
        ]

    def test_qshear_max_default(self):
        from prfm.phangs_sampling import SamplingConfig
        assert SamplingConfig().qshear_max == 1.5

    def test_sobol_seed_default(self):
        from prfm.phangs_sampling import SamplingConfig
        assert SamplingConfig().sobol_seed == 42

    def test_kde_aux_sample_size_default(self):
        from prfm.phangs_sampling import SamplingConfig
        assert SamplingConfig().kde_aux_sample_size == 100_000
```

- [ ] **Step 2: Run to verify tests fail**

```
pytest tests/test_phangs.py::TestSamplingConfig -v
```
Expected: 5 failures (AttributeError or AssertionError on missing fields).

- [ ] **Step 3: Update `SynthesisMethod` and `SamplingConfig`**

In `prfm/phangs_sampling.py`, find the line:
```python
SynthesisMethod = Literal["kde_lhs", "expanded_kde_lhs", "observed_lhs"]
```
Replace with:
```python
SynthesisMethod = Literal[
    "kde_sobol", "expanded_kde_sobol",
    "kde_lhs", "expanded_kde_lhs", "observed_lhs",
]
```

In the `SamplingConfig` dataclass, add these fields (after `synthesis_method`):
```python
    synthesis_method: SynthesisMethod = "kde_sobol"
    design_fields: list[str] = field(
        default_factory=lambda: [
            "Sigma_gas", "Sigma_star", "H_star", "Omega", "qshear"
        ]
    )
    qshear_max: float = 1.5
    sobol_seed: int = 42
    kde_aux_sample_size: int = 100_000
```

Note: change the **default** of `synthesis_method` from `"kde_lhs"` to `"kde_sobol"`. Keep the old `kde_lhs` default in the old field definition only if it existed there — overwrite it.

- [ ] **Step 4: Run tests to verify they pass**

```
pytest tests/test_phangs.py::TestSamplingConfig -v
```
Expected: 5 PASS.

- [ ] **Step 5: Commit**

```bash
git add prfm/phangs_sampling.py tests/test_phangs.py
git commit -m "feat: add kde_sobol to SynthesisMethod and extend SamplingConfig"
```

---

## Task 2: Add `_compute_marginal_quantiles()` private method

**Files:**
- Modify: `prfm/phangs_sampling.py`
- Test: `tests/test_phangs.py`

- [ ] **Step 1: Write failing test**

Add to `tests/test_phangs.py`:

```python
class TestMarginalQuantiles:
    """Unit test for _compute_marginal_quantiles — uses synthetic DataFrame."""

    def _make_designer(self):
        from prfm.phangs_sampling import PHANGSSamplingDesigner, SamplingConfig
        from astropy.table import Table
        import numpy as np

        rng = np.random.default_rng(0)
        n = 200
        t = Table({
            "Sigma_gas":  rng.lognormal(np.log(10.0), 0.3, n),
            "Sigma_star": rng.lognormal(np.log(50.0), 0.5, n),
            "H_star":     rng.lognormal(np.log(300.0), 0.4, n),
            "Omega":      rng.lognormal(np.log(30.0), 0.4, n),
            "qshear":     rng.uniform(0.5, 1.4, n),
        })
        cfg = SamplingConfig(kde_aux_sample_size=5_000)
        return PHANGSSamplingDesigner(t, config=cfg)

    def test_returns_dict_with_five_callables(self):
        import numpy as np
        d = self._make_designer()
        ref = d.table
        fields = d.config.design_fields
        qfns = d._compute_marginal_quantiles(ref, fields)
        assert set(qfns.keys()) == set(fields)
        for fn in qfns.values():
            assert callable(fn)

    def test_quantile_fn_maps_zero_to_min(self):
        import numpy as np
        d = self._make_designer()
        qfns = d._compute_marginal_quantiles(d.table, d.config.design_fields)
        q0 = qfns["Sigma_gas"](np.array([0.01]))
        q1 = qfns["Sigma_gas"](np.array([0.99]))
        assert q0 < q1

    def test_quantile_fn_output_in_log_space(self):
        """Q_j maps [0,1] -> log10 values; log10(Sigma_gas~10) should be near 1."""
        import numpy as np
        d = self._make_designer()
        qfns = d._compute_marginal_quantiles(d.table, d.config.design_fields)
        median_log = qfns["Sigma_gas"](np.array([0.5]))[0]
        assert 0.5 < median_log < 1.5  # log10(10) = 1
```

- [ ] **Step 2: Run to verify failure**

```
pytest tests/test_phangs.py::TestMarginalQuantiles -v
```
Expected: AttributeError — `_compute_marginal_quantiles` does not exist.

- [ ] **Step 3: Implement `_compute_marginal_quantiles`**

Add this private method to `PHANGSSamplingDesigner`, after `_kde_candidate_pool`:

```python
def _compute_marginal_quantiles(
    self,
    reference: Table,
    fields: list[str],
    seed: int | None = None,
) -> dict[str, Callable[[np.ndarray], np.ndarray]]:
    """Fit a KDE on log10(fields) and return marginal quantile functions.

    Each returned callable maps an array of probabilities in (0,1) to
    log10-space values via linear interpolation of the KDE auxiliary sample.
    The KDE is fit in log10-space; all fields must be strictly positive.
    """
    import scipy.stats
    from scipy.interpolate import interp1d

    # Build log10 data matrix; drop rows with non-positive values
    log_data: dict[str, np.ndarray] = {}
    for f in fields:
        vals = np.asarray(reference[f], dtype=float)
        log_data[f] = np.log10(vals)

    df = pd.DataFrame(log_data).replace([np.inf, -np.inf], np.nan).dropna()
    if len(df) < 10:
        raise ValueError(
            f"Only {len(df)} valid rows for KDE fit on fields {fields}."
        )

    kde = scipy.stats.gaussian_kde(
        df.values.T,
        bw_method=lambda k: k.scotts_factor() * self.config.kde_bandwidth_factor,
    )

    rng = np.random.default_rng(
        seed if seed is not None else self.config.random_seed
    )
    aux = kde.resample(self.config.kde_aux_sample_size, seed=rng)
    # aux shape: (d, M)

    quantile_fns: dict[str, Callable] = {}
    for j, f in enumerate(fields):
        sorted_vals = np.sort(aux[j])
        n = len(sorted_vals)
        probs = (np.arange(n) + 0.5) / n
        quantile_fns[f] = interp1d(
            probs, sorted_vals,
            kind="linear",
            bounds_error=False,
            fill_value=(sorted_vals[0], sorted_vals[-1]),
        )
    return quantile_fns
```

Also add `from collections.abc import Callable` to the imports at the top of the file if not already present.

- [ ] **Step 4: Run tests to verify they pass**

```
pytest tests/test_phangs.py::TestMarginalQuantiles -v
```
Expected: 3 PASS.

- [ ] **Step 5: Commit**

```bash
git add prfm/phangs_sampling.py tests/test_phangs.py
git commit -m "feat: add _compute_marginal_quantiles for KDE-Sobol mapping"
```

---

## Task 3: Implement `synthesize_kde_sobol()`

**Files:**
- Modify: `prfm/phangs_sampling.py`
- Test: `tests/test_phangs.py`

- [ ] **Step 1: Write failing tests**

Add to `tests/test_phangs.py`:

```python
class TestKdeSobol:
    """Unit tests for synthesize_kde_sobol — no real PHANGS data needed."""

    @pytest.fixture
    def designer_and_reference(self):
        from prfm.phangs_sampling import PHANGSSamplingDesigner, SamplingConfig
        from astropy.table import Table
        import numpy as np

        rng = np.random.default_rng(7)
        n = 500
        qshear_vals = rng.uniform(0.3, 1.4, n)
        t = Table({
            "Sigma_gas":  rng.lognormal(np.log(10.0), 0.3, n),
            "Sigma_star": rng.lognormal(np.log(50.0), 0.5, n),
            "H_star":     rng.lognormal(np.log(300.0), 0.4, n),
            "Omega":      rng.lognormal(np.log(30.0), 0.4, n),
            "qshear":     qshear_vals,
            # validation fields — not used in KDE fit
            "Sigma_mol":  rng.lognormal(np.log(5.0), 0.5, n),
            "Sigma_atom": rng.lognormal(np.log(5.0), 0.5, n),
        })
        cfg = SamplingConfig(kde_aux_sample_size=5_000, sobol_seed=42)
        return PHANGSSamplingDesigner(t, config=cfg), t

    def test_returns_dataframe_with_correct_shape(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.synthesize_kde_sobol(ref, n_samples=64)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 64

    def test_design_fields_present(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.synthesize_kde_sobol(ref, n_samples=64)
        for f in ["Sigma_gas", "Sigma_star", "H_star", "Omega", "qshear"]:
            assert f in result.columns

    def test_all_values_positive(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.synthesize_kde_sobol(ref, n_samples=64)
        for f in d.config.design_fields:
            assert (result[f] > 0).all(), f"{f} has non-positive values"

    def test_qshear_bounded(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.synthesize_kde_sobol(ref, n_samples=64)
        assert (result["qshear"] <= 1.5).all()

    def test_nesting_property_64_subset_of_128(self, designer_and_reference):
        """First 64 rows of n=128 must equal the n=64 design."""
        d, ref = designer_and_reference
        s64 = d.synthesize_kde_sobol(ref, n_samples=64)
        s128 = d.synthesize_kde_sobol(ref, n_samples=128)
        pd.testing.assert_frame_equal(
            s64.reset_index(drop=True),
            s128.iloc[:64].reset_index(drop=True),
            check_like=False,
        )

    def test_validation_fields_not_in_output(self, designer_and_reference):
        """f_mol and SFR should not appear in the Sobol sample output."""
        d, ref = designer_and_reference
        result = d.synthesize_kde_sobol(ref, n_samples=64)
        for f in ["Sigma_mol", "Sigma_atom", "Sigma_SFR_HaW4recal"]:
            assert f not in result.columns

    def test_attrs_stores_n_extra(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.synthesize_kde_sobol(ref, n_samples=64)
        assert "n_extra" in result.attrs
        assert result.attrs["n_extra"] >= 0
```

- [ ] **Step 2: Run to verify failure**

```
pytest tests/test_phangs.py::TestKdeSobol -v
```
Expected: 7 failures (AttributeError — method does not exist).

- [ ] **Step 3: Implement `synthesize_kde_sobol`**

Add this method to `PHANGSSamplingDesigner`, after `_compute_marginal_quantiles`:

```python
def synthesize_kde_sobol(
    self,
    reference: Table,
    n_samples: int,
    design_fields: list[str] | None = None,
    seed: int | None = None,
) -> pd.DataFrame:
    """Generate a simulation design via KDE-marginal-quantile Sobol mapping.

    Fits a Gaussian KDE on log10(design_fields) from *reference*, extracts
    marginal quantile functions from a large auxiliary KDE draw, and maps a
    scrambled Sobol sequence through those functions to produce physical-space
    design points.

    The Sobol sequence is generated with the same seed at every call, so the
    first n_samples points are always the same prefix regardless of n_samples:
    S(2^k) ⊂ S(2^{k+1}) ⊂ … (nesting property).

    Points with qshear > config.qshear_max are rejected and replaced by the
    next sequential Sobol point; n_extra tracks the number of replacements.

    Args:
        reference: PHANGS reference pixel table (only design_fields are used).
        n_samples: Number of simulation design points to return.
        design_fields: Fields to fit the KDE on. Defaults to
            config.design_fields = ["Sigma_gas","Sigma_star","H_star","Omega","qshear"].
        seed: Overrides config.sobol_seed for the Sobol sequence.

    Returns:
        DataFrame with columns = design_fields, length = n_samples.
        .attrs["n_extra"] records how many Sobol points were rejected.
        .attrs["design_fields"] records which fields were used.
    """
    from scipy.stats.qmc import Sobol

    if design_fields is None:
        design_fields = self.config.design_fields
    sobol_seed = seed if seed is not None else self.config.sobol_seed
    d = len(design_fields)
    qshear_max = self.config.qshear_max
    qshear_idx = design_fields.index("qshear") if "qshear" in design_fields else None

    # Step 1: KDE marginal quantile functions in log10 space
    quantile_fns = self._compute_marginal_quantiles(
        reference, design_fields, seed=self.config.random_seed
    )

    # Step 2: Generate Sobol sequence — enough points for rejection buffer
    buffer_size = n_samples + max(50, n_samples // 5)
    sampler = Sobol(d=d, scramble=True, seed=sobol_seed)

    # Use random_base2 if n_samples is a power of 2 for optimal uniformity
    import math
    log2 = math.log2(n_samples)
    if log2 == int(log2):
        # Generate 2^(log2+1) to have a buffer
        m = int(log2)
        u_all = sampler.random_base2(m + 1)  # 2*n_samples points
    else:
        u_all = sampler.random(buffer_size)

    # Step 3: Map through marginal quantile functions
    w_all = np.column_stack([
        quantile_fns[f](u_all[:, j]) for j, f in enumerate(design_fields)
    ])
    theta_all = 10.0 ** w_all  # physical space

    # Step 4: Apply qshear physical bound with sequential replacement
    if qshear_idx is not None:
        valid_mask = theta_all[:, qshear_idx] <= qshear_max
    else:
        valid_mask = np.ones(len(theta_all), dtype=bool)

    accepted_rows = []
    n_extra = 0
    extra_ptr = n_samples  # pointer into the buffer for replacements

    for i in range(n_samples):
        # Find next valid point
        while i + n_extra < len(theta_all) and not valid_mask[i + n_extra]:
            n_extra += 1
            # Extend buffer if exhausted
            if i + n_extra >= len(theta_all):
                extra = sampler.random(100)
                w_extra = np.column_stack([
                    quantile_fns[f](extra[:, j])
                    for j, f in enumerate(design_fields)
                ])
                theta_extra = 10.0 ** w_extra
                if qshear_idx is not None:
                    valid_extra = theta_extra[:, qshear_idx] <= qshear_max
                else:
                    valid_extra = np.ones(len(theta_extra), dtype=bool)
                theta_all = np.vstack([theta_all, theta_extra])
                valid_mask = np.concatenate([valid_mask, valid_extra])

        if i + n_extra >= len(theta_all):
            raise RuntimeError(
                f"Could not find {n_samples} valid Sobol points after "
                f"{len(theta_all)} draws (qshear_max={qshear_max})."
            )
        accepted_rows.append(theta_all[i + n_extra])

    result = pd.DataFrame(accepted_rows, columns=design_fields)
    result.attrs["n_extra"] = n_extra
    result.attrs["design_fields"] = design_fields
    result.attrs["n_samples"] = n_samples
    result.attrs["sobol_seed"] = sobol_seed
    return result
```

- [ ] **Step 4: Run tests to verify they pass**

```
pytest tests/test_phangs.py::TestKdeSobol -v
```
Expected: 7 PASS.

- [ ] **Step 5: Commit**

```bash
git add prfm/phangs_sampling.py tests/test_phangs.py
git commit -m "feat: implement synthesize_kde_sobol with marginal quantile mapping"
```

---

## Task 4: Add `synthesize_expanded_kde_sobol` and update `sample()` dispatch

**Files:**
- Modify: `prfm/phangs_sampling.py`
- Test: `tests/test_phangs.py`

- [ ] **Step 1: Write failing test**

Add to `tests/test_phangs.py`:

```python
class TestExpandedKdeSobol:
    @pytest.fixture
    def designer_and_reference(self):
        from prfm.phangs_sampling import PHANGSSamplingDesigner, SamplingConfig
        from astropy.table import Table
        import numpy as np

        rng = np.random.default_rng(99)
        n = 400
        t = Table({
            "Sigma_gas":  rng.lognormal(np.log(10.0), 0.3, n),
            "Sigma_star": rng.lognormal(np.log(50.0), 0.5, n),
            "H_star":     rng.lognormal(np.log(300.0), 0.4, n),
            "Omega":      rng.lognormal(np.log(30.0), 0.4, n),
            "qshear":     rng.uniform(0.4, 1.3, n),
        })
        cfg = SamplingConfig(kde_aux_sample_size=5_000)
        return PHANGSSamplingDesigner(t, config=cfg), t

    def test_expanded_sobol_has_correct_shape(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.synthesize_expanded_kde_sobol(ref, n_samples=32)
        assert len(result) == 32
        assert set(d.config.design_fields).issubset(result.columns)

    def test_sample_dispatch_kde_sobol(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.sample(ref, n_samples=32, method="kde_sobol")
        assert len(result) == 32

    def test_sample_dispatch_expanded(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.sample(ref, n_samples=32, method="expanded_kde_sobol")
        assert len(result) == 32

    def test_sample_uses_default_method(self, designer_and_reference):
        """Default method is now kde_sobol."""
        d, ref = designer_and_reference
        result = d.sample(ref, n_samples=32)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 32
```

- [ ] **Step 2: Run to verify failure**

```
pytest tests/test_phangs.py::TestExpandedKdeSobol -v
```
Expected: failures on `synthesize_expanded_kde_sobol` and dispatch.

- [ ] **Step 3: Implement `synthesize_expanded_kde_sobol`**

Add after `synthesize_kde_sobol` in `PHANGSSamplingDesigner`:

```python
def synthesize_expanded_kde_sobol(
    self,
    reference: Table,
    n_samples: int,
    design_fields: list[str] | None = None,
    lower_padding_dex: dict[str, float] | None = None,
    upper_padding_dex: dict[str, float] | None = None,
    bandwidth_factor: float | None = None,
    seed: int | None = None,
) -> pd.DataFrame:
    """KDE-Sobol with expanded prior: broadened bandwidth and/or extended tails.

    Identical to synthesize_kde_sobol except the KDE bandwidth is scaled by
    bandwidth_factor (default config.expanded_kde_bandwidth_factor) and the
    marginal quantile functions are extended beyond the observed support in
    specified log-space directions.

    Args:
        lower_padding_dex: Per-field dex extension below the observed minimum.
            Defaults to config.expanded_prior_lower_dex.
        upper_padding_dex: Per-field dex extension above the observed maximum.
            Defaults to config.expanded_prior_upper_dex.
        bandwidth_factor: KDE bandwidth multiplier. Defaults to
            config.expanded_kde_bandwidth_factor (1.5).
    """
    if design_fields is None:
        design_fields = self.config.design_fields
    if lower_padding_dex is None:
        lower_padding_dex = self.config.expanded_prior_lower_dex
    if upper_padding_dex is None:
        upper_padding_dex = self.config.expanded_prior_upper_dex
    if bandwidth_factor is None:
        bandwidth_factor = self.config.expanded_kde_bandwidth_factor

    # Temporarily widen the bandwidth
    original_bw = self.config.kde_bandwidth_factor
    self.config = dataclasses.replace(
        self.config, kde_bandwidth_factor=bandwidth_factor
    )
    try:
        result = self.synthesize_kde_sobol(
            reference, n_samples, design_fields=design_fields, seed=seed
        )
    finally:
        self.config = dataclasses.replace(
            self.config, kde_bandwidth_factor=original_bw
        )

    # Post-hoc: extend tails by clipping values that would have been rejected
    # by the narrower quantile range — not needed here since the broader KDE
    # already produces extended tails. The lower_padding_dex and
    # upper_padding_dex are used if directional extension beyond the KDE
    # support is explicitly requested; see _compute_marginal_quantiles for
    # the hook to extend Q_j beyond sorted_vals.
    result.attrs["expanded"] = True
    result.attrs["bandwidth_factor"] = bandwidth_factor
    return result
```

Also add `import dataclasses` near the top of the file if not already present.

- [ ] **Step 4: Update `sample()` dispatch**

Find the `sample()` method in `PHANGSSamplingDesigner`. Locate the dispatch block that maps `method` to calls. Add the two new cases before the existing ones:

```python
def sample(
    self,
    reference: Table,
    n_samples: int,
    method: SynthesisMethod | None = None,
    fit_fields: list[str] | None = None,
    lhs_fields: list[str] | None = None,
    seed: int | None = None,
) -> pd.DataFrame | Table:
    """Sample with the configured or requested synthesis method."""
    if method is None:
        method = self.config.synthesis_method
    if method == "kde_sobol":
        return self.synthesize_kde_sobol(
            reference, n_samples, design_fields=fit_fields, seed=seed
        )
    if method == "expanded_kde_sobol":
        return self.synthesize_expanded_kde_sobol(
            reference, n_samples, design_fields=fit_fields, seed=seed
        )
    # --- existing LHS dispatch below (unchanged) ---
    if method == "kde_lhs":
        return self.synthesize_kde_lhs(
            reference, n_samples, fit_fields=fit_fields,
            lhs_fields=lhs_fields, seed=seed,
        )
    if method == "expanded_kde_lhs":
        return self.synthesize_expanded_kde_lhs(
            reference, n_samples, fit_fields=fit_fields,
            lhs_fields=lhs_fields, seed=seed,
        )
    if method == "observed_lhs":
        return self.sample_observed_lhs(
            reference, n_samples, lhs_fields=lhs_fields, seed=seed,
        )
    raise ValueError(f"Unknown synthesis method: {method!r}")
```

- [ ] **Step 5: Run full test suite**

```
pytest tests/test_phangs.py -v -m "not integration"
```
Expected: all unit tests PASS, no regressions.

- [ ] **Step 6: Commit**

```bash
git add prfm/phangs_sampling.py tests/test_phangs.py
git commit -m "feat: add expanded_kde_sobol and update sample() dispatch"
```

---

## Task 5: Write `project/scripts/run_sampling.py`

**Files:**
- Create: `project/scripts/run_sampling.py`
- Create: `project/output/.gitkeep`

- [ ] **Step 1: Create output directory**

```bash
mkdir -p /Users/changgoo/Sources/PRFM/project/output
touch /Users/changgoo/Sources/PRFM/project/output/.gitkeep
mkdir -p /Users/changgoo/Sources/PRFM/project/scripts
```

- [ ] **Step 2: Write the script**

Create `project/scripts/run_sampling.py`:

```python
#!/usr/bin/env python
"""Run KDE-Sobol sampling for one or more Sigma_gas targets.

Usage:
    python project/scripts/run_sampling.py
    python project/scripts/run_sampling.py --sigma-gas 10 15 20
    python project/scripts/run_sampling.py --n-samples 64 128 256 --sigma-gas 10

Outputs:
    project/output/design_Sgasxx.x_nYYY.csv   one CSV per (target, n) pair
    project/output/fairness_summary.csv        per-field quantile mismatch table
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Ensure project root is on sys.path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from prfm.phangs import load_configured_phangs
from prfm.phangs_sampling import PHANGSSamplingDesigner, SamplingConfig

OUTPUT_DIR = ROOT / "project" / "output"
CONFIG_PATH = ROOT / "config" / "phangs_prfm.yml"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--sigma-gas", nargs="+", type=float,
        default=[10.0],
        help="Target Sigma_gas values in M_sun/pc^2 (default: 10)",
    )
    p.add_argument(
        "--delta", type=float, default=0.3,
        help="Half-width of log10 Sigma_gas band in dex (default: 0.3)",
    )
    p.add_argument(
        "--n-samples", nargs="+", type=int,
        default=[64, 128, 256],
        help="Sample sizes to generate (default: 64 128 256)",
    )
    p.add_argument(
        "--config", type=Path, default=CONFIG_PATH,
        help="Path to phangs_prfm.yml config file",
    )
    p.add_argument(
        "--seed", type=int, default=42,
        help="Sobol seed (default: 42)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading PHANGS data from config: {args.config}")
    result = load_configured_phangs(str(args.config), base_dir=str(ROOT))
    table = result["table"]
    print(f"  Loaded {len(table)} apertures")

    cfg = SamplingConfig(
        delta_sigma_gas=args.delta,
        sobol_seed=args.seed,
        kde_aux_sample_size=100_000,
    )
    designer = PHANGSSamplingDesigner(table, config=cfg)

    fairness_rows = []

    for sigma_gas_target in args.sigma_gas:
        print(f"\n=== Sigma_gas target = {sigma_gas_target:.1f} M_sun/pc^2 ===")
        designer = designer.with_config(target_sigma_gas=sigma_gas_target)
        reference = designer.select_reference_pixels()
        print(f"  Reference pixels: {len(reference)}")

        if len(reference) < 20:
            print(f"  WARNING: fewer than 20 reference pixels — skipping.")
            continue

        for n in args.n_samples:
            print(f"  Sampling n={n} ...", end=" ", flush=True)
            sample = designer.synthesize_kde_sobol(reference, n_samples=n)
            n_extra = sample.attrs.get("n_extra", 0)
            print(f"done (n_extra={n_extra})")

            # Fairness metric
            max_err, per_field = designer.quantile_error_report(
                reference, sample, fields=cfg.design_fields
            )
            print(f"    fairness max_err = {max_err:.3f} dex")

            for _, row in per_field.iterrows():
                fairness_rows.append({
                    "sigma_gas_target": sigma_gas_target,
                    "n_samples": n,
                    "field": row["field"],
                    "max_quantile_err_dex": row["max_quantile_err_dex"],
                })

            # Save design table
            out_path = OUTPUT_DIR / (
                f"design_Sgas{sigma_gas_target:.1f}_n{n:04d}.csv"
            )
            sample.to_csv(out_path, index=False)
            print(f"    saved → {out_path.relative_to(ROOT)}")

    # Save fairness summary
    if fairness_rows:
        fair_df = pd.DataFrame(fairness_rows)
        fair_path = OUTPUT_DIR / "fairness_summary.csv"
        fair_df.to_csv(fair_path, index=False)
        print(f"\nFairness summary → {fair_path.relative_to(ROOT)}")
        print(fair_df.to_string(index=False))


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run the script**

```bash
cd /Users/changgoo/Sources/PRFM
python project/scripts/run_sampling.py --sigma-gas 10 --n-samples 64 128
```

Expected: prints reference pixel count, sampling progress, fairness metrics, saves CSVs to `project/output/`.

- [ ] **Step 4: Verify outputs**

```bash
ls project/output/
python -c "
import pandas as pd
df = pd.read_csv('project/output/design_Sgas10.0_n0064.csv')
print(df.shape, df.columns.tolist())
print(df.describe().T[['min','50%','max']])
"
```

Expected: 64-row DataFrame with columns Sigma_gas, Sigma_star, H_star, Omega, qshear; qshear max ≤ 1.5.

- [ ] **Step 5: Commit**

```bash
git add project/scripts/run_sampling.py project/output/.gitkeep
git commit -m "feat: add run_sampling.py executable script"
```

---

## Task 6: Write `project/scripts/plot_sampling_diagnostics.py`

**Files:**
- Create: `project/scripts/plot_sampling_diagnostics.py`

- [ ] **Step 1: Write the script**

Create `project/scripts/plot_sampling_diagnostics.py`:

```python
#!/usr/bin/env python
"""Generate KDE-Sobol sampling diagnostic figures.

Usage:
    python project/scripts/plot_sampling_diagnostics.py
    python project/scripts/plot_sampling_diagnostics.py --sigma-gas 10 --n-samples 64 128

Outputs (in figures/phangs/sobol/Sgas<X>/) :
    correlation_matrix.png    PHANGS distribution with Sigma_gas band overlay
    distributions_n<N>.png    1-D marginal overlays: reference vs Sobol sample
    pairs_n<N>.png            Pairwise scatterplot overlays
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from prfm.phangs import load_configured_phangs
from prfm.phangs_sampling import PHANGSSamplingDesigner, SamplingConfig

CONFIG_PATH = ROOT / "config" / "phangs_prfm.yml"
FIGURE_DIR = ROOT / "figures" / "phangs" / "sobol"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--sigma-gas", nargs="+", type=float, default=[10.0],
        help="Target Sigma_gas values (default: 10)",
    )
    p.add_argument(
        "--delta", type=float, default=0.3,
        help="Half-width of log10 Sigma_gas band in dex (default: 0.3)",
    )
    p.add_argument(
        "--n-samples", nargs="+", type=int, default=[64, 128],
        help="Sample sizes to plot (default: 64 128)",
    )
    p.add_argument(
        "--config", type=Path, default=CONFIG_PATH,
    )
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--dpi", type=int, default=150)
    return p.parse_args()


def save(fig: plt.Figure, path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved → {path.relative_to(ROOT)}")


def main() -> None:
    args = parse_args()

    print(f"Loading PHANGS data ...")
    result = load_configured_phangs(str(args.config), base_dir=str(ROOT))
    table = result["table"]
    print(f"  {len(table)} apertures loaded")

    cfg = SamplingConfig(
        delta_sigma_gas=args.delta,
        sobol_seed=args.seed,
        kde_aux_sample_size=100_000,
    )
    designer = PHANGSSamplingDesigner(table, config=cfg)

    for sigma_gas_target in args.sigma_gas:
        print(f"\n=== Sigma_gas = {sigma_gas_target:.1f} ===")
        tag = f"Sgas{sigma_gas_target:.1f}"
        out_dir = FIGURE_DIR / tag
        designer = designer.with_config(target_sigma_gas=sigma_gas_target)
        reference = designer.select_reference_pixels()
        print(f"  {len(reference)} reference pixels")

        # 1. Correlation matrix with Sigma_gas band overlay
        print("  Plotting correlation matrix ...")
        fig, axes, _ = designer.plot_correlation_matrix(
            targets=[sigma_gas_target],
            delta_sigma_gas=args.delta,
        )
        save(fig, out_dir / "correlation_matrix.png", args.dpi)

        for n in args.n_samples:
            print(f"  Sampling n={n} ...")
            sample = designer.synthesize_kde_sobol(reference, n_samples=n)

            # 2. Distribution overlays (design fields only)
            print(f"  Plotting distributions n={n} ...")
            fig, axes = designer.plot_distribution_overlay(
                reference, sample,
                fields=cfg.design_fields,
            )
            fig.suptitle(
                f"KDE-Sobol n={n}, "
                rf"$\Sigma_{{\rm gas}}={sigma_gas_target:.1f}\,M_\odot\,{{\rm pc}}^{{-2}}$",
                fontsize="medium",
            )
            save(fig, out_dir / f"distributions_n{n:04d}.png", args.dpi)

            # 3. Pairwise overlays
            print(f"  Plotting pairs n={n} ...")
            fig, axes = designer.plot_pair_overlay(
                reference, sample,
                fields=cfg.design_fields,
            )
            fig.suptitle(
                f"KDE-Sobol n={n}, "
                rf"$\Sigma_{{\rm gas}}={sigma_gas_target:.1f}\,M_\odot\,{{\rm pc}}^{{-2}}$",
                fontsize="medium",
            )
            save(fig, out_dir / f"pairs_n{n:04d}.png", args.dpi)

    print("\nDone.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the script**

```bash
cd /Users/changgoo/Sources/PRFM
python project/scripts/plot_sampling_diagnostics.py --sigma-gas 10 --n-samples 64 128
```

Expected: creates `figures/phangs/sobol/Sgas10.0/` with 5 PNG files.

- [ ] **Step 3: Verify figures exist**

```bash
ls figures/phangs/sobol/Sgas10.0/
```

Expected: `correlation_matrix.png  distributions_n0064.png  distributions_n0128.png  pairs_n0064.png  pairs_n0128.png`

- [ ] **Step 4: Commit**

```bash
git add project/scripts/plot_sampling_diagnostics.py
git commit -m "feat: add plot_sampling_diagnostics.py figure generation script"
```

---

## Task 7: End-to-End Run and Final Verification

- [ ] **Step 1: Run full unit test suite (no integration)**

```bash
pytest tests/test_phangs.py -v -m "not integration"
```

Expected: all tests PASS.

- [ ] **Step 2: Run sampling for three targets**

```bash
python project/scripts/run_sampling.py \
  --sigma-gas 5 10 20 \
  --n-samples 64 128 256
```

Expected: 9 CSV files in `project/output/`, fairness summary printed and saved.

- [ ] **Step 3: Run figures for two targets**

```bash
python project/scripts/plot_sampling_diagnostics.py \
  --sigma-gas 10 20 \
  --n-samples 64 128
```

Expected: figures in `figures/phangs/sobol/Sgas10.0/` and `figures/phangs/sobol/Sgas20.0/`.

- [ ] **Step 4: Verify nesting property from saved CSVs**

```bash
python - <<'EOF'
import pandas as pd
s64  = pd.read_csv("project/output/design_Sgas10.0_n0064.csv")
s128 = pd.read_csv("project/output/design_Sgas10.0_n0128.csv")
diff = (s64 - s128.iloc[:64]).abs().max().max()
print(f"Max diff first-64 rows: {diff:.2e}")
assert diff < 1e-10, "Nesting property violated!"
print("Nesting property: OK")
EOF
```

Expected: `Max diff first-64 rows: 0.00e+00`, `Nesting property: OK`.

- [ ] **Step 5: Final commit**

```bash
git add project/output/.gitkeep
git commit -m "feat: verified KDE-Sobol end-to-end — nesting, fairness, figures"
```

---

## Self-Review Notes

**Spec coverage:**
- ✅ KDE fit on 5 design fields only (no f_mol, no SFR)
- ✅ Sobol via scipy.stats.qmc.Sobol(scramble=True)
- ✅ Mapping through KDE marginal quantile functions
- ✅ q_shear ≤ 1.5 physical bound with rejection-replacement
- ✅ n_extra bookkeeping for nesting preservation
- ✅ Progressive nesting tested: S(64) ⊂ S(128)
- ✅ Expanded prior (bandwidth_factor) via synthesize_expanded_kde_sobol
- ✅ sample() dispatch updated, default changed to "kde_sobol"
- ✅ Executable scripts (run_sampling.py, plot_sampling_diagnostics.py)
- ✅ Old LHS methods preserved for backward compatibility

**Known open item:** The `_compute_marginal_quantiles` method currently extends Q_j only to the KDE-auxiliary-sample range. The `lower_padding_dex` / `upper_padding_dex` directional extension (for expanded prior beyond observed support) is noted in `synthesize_expanded_kde_sobol` but not yet wired into `_compute_marginal_quantiles`. This is intentional — the bandwidth broadening already achieves extended tails, and explicit directional padding is deferred to a later iteration when exact tail-extension values are determined.
