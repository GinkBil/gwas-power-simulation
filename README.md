# GWAS Power Simulation Pipeline

**Simulation-based evaluation of GWAS detection under population structure using Python, SLiM, and GEMMA.**

This research pipeline generates simulated genomic and quantitative-phenotype data, performs genome-wide association analyses using a linear mixed model (LMM), and records association statistics for a known causal locus alongside counts of significant associations outside a specified genomic window. It was developed in the context of research on GWAS statistical power and sample-size requirements.

The files in this repository represent **one population-structure scenario** within that broader research: an ancestral population splits into two populations, followed by a one-time admixture event that establishes a third population. The current notebook samples individuals from this admixed population at multiple sample sizes.

## Research questions

- How does the number of sampled individuals affect the detectability of a causal genetic variant?
- How can population structure and relatedness be represented in simulated GWAS data?
- How many significant associations occur at the causal locus versus elsewhere in the genome under a specified significance threshold?

The scripts produce **per-replicate results** that can be aggregated across independent simulation runs to estimate empirical detection power.

## Workflow

1. **Forward simulation (SLiM):** Simulate a chromosome and introduce a selected mutation at a known genomic position. Begin with a single population (`p1`), split off `p2`, then establish an admixed population (`p3`) through a single 50:50 contribution from `p1` and `p2`.
2. **Tree-sequence processing (pyslim, tskit, msprime):** Recapitate ancestral genealogies, overlay neutral mutations, and extract population-level genetic summaries.
3. **Sampling and phenotype simulation (Python):** Sample individuals across a range of sizes; calculate genotypes at the causal site; generate quantitative phenotypes from a genetic effect plus normally distributed noise.
4. **GWAS data preparation (PLINK):** Export sampled genotypes to VCF and convert them to PLINK binary files with associated phenotype values.
5. **Association analysis (GEMMA):** Compute a kinship matrix and run a linear mixed-model association test.
6. **Result summaries (pandas):** Record association results for the causal locus and summarize significant associations inside and outside a user-defined window around that locus.

### Phenotype model

The simulation uses the additive model

\[
y_i = \beta g_i + \varepsilon_i, \qquad \varepsilon_i \sim \mathcal{N}(0, \sigma^2),
\]

where \(g_i\) is an individual's dosage at the causal locus. The code calculates an effect size from a specified heritability target, genotype variance, and environmental-noise variance. The current implementation uses genotype variance from population `p1` to calibrate that effect, while the default GWAS samples individuals from `p3`.

## Repository contents

| File | Description |
| --- | --- |
| `hybrid_constant_freqsel.txt` | SLiM recipe for a three-phase demographic scenario with frequency-dependent selection on a known mutation. |
| `gwas_simulation_functions.py` | Helper functions for simulation execution, tree-sequence processing, sampling, phenotype generation, data conversion, GWAS execution, and result summaries. |
| `gwas_run_simulations_vMay1.ipynb` | Jupyter notebook containing parameters, orchestration of the pipeline, and a one-replicate example call. |

> The notebook also contains a copy of the helper-function definitions. The simulation driver calls the functions imported from `gwas_simulation_functions.py` as `gsim`.

## Default simulation configuration

The supplied notebook currently uses the following settings:

| Parameter | Value |
| --- | --- |
| Population size (`N`) | 500 |
| Chromosome length | 90,000,000 bases |
| Causal mutation position | 45,000,000 |
| Duration before split | 500 generations |
| Duration between split and contact | 500 generations |
| Duration after contact | 100 generations |
| Frequency-dependent selection target (`ftar`) | 0.5 |
| Frequency-response coefficient (`sel_speed`) | 1.0 |
| Admixture at contact | 50% `p1`, 50% `p2` |
| Sampled population | `p3` only |
| Sample sizes | 20, 50, 80, 100, 150, 200, 250, 300, 350, 400, 450, 500 |
| Target heritability input (`h2`) | 0.05 |
| Environmental noise standard deviation | 1 |
| Significance threshold | 5 × 10⁻⁸ |
| Causal-region window | ±100,000 bases |
| Recapitation ancestral effective size | 500 |
| Recapitation recombination rate | 1 × 10⁻⁸ |
| Neutral mutation rate | 1 × 10⁻⁸ |

**Selection note:** The active SLiM `mutationEffect` implementation uses a frequency-dependent coefficient, `s = (ftar - f) * sel_speed`, during all three phases. The Python parameter dictionary also contains `s_p1`, `s_p2`, and `mig_to_p3`, but the supplied SLiM recipe does **not** currently use these values; the contact event is explicitly set to a 50:50 admixture pulse.

## Requirements

**Python packages**

- `numpy`
- `pandas`
- `tskit`
- `pyslim`
- `msprime`
- `jupyter` (to open the notebook)

**External command-line tools**

- [SLiM](https://messerlab.org/slim/) — forward-in-time genomic simulation
- [PLINK](https://www.cog-genomics.org/plink/) — VCF to binary genotype conversion
- [GEMMA](https://github.com/genetics-statistics/GEMMA) — kinship estimation and LMM association analysis

Install compatible releases of the external tools and ensure that their executables (`slim`, `plink`, and `gemma`) are accessible on your `PATH`, or modify the executable paths in the notebook. The source files do not include pinned dependency versions or a verified environment specification.

## Running the notebook

1. Place the three source files together in a working directory.
2. Install the Python dependencies in an appropriate Python environment, for example:

   ```bash
   python -m pip install numpy pandas tskit pyslim msprime jupyter
   ```

3. Install SLiM, PLINK, and GEMMA separately.
4. Open the notebook:

   ```bash
   jupyter notebook gwas_run_simulations_vMay1.ipynb
   ```

5. Review the parameters in the configuration cell, especially `BASE_DIR`, `sample_sizes`, `params_dict`, and the executable paths.
6. Review the code notes below before running. The notebook's example calls `run_one_simulation(1)`.

The notebook has **not been independently verified in a clean environment** for this public documentation; successful execution depends on compatible packages, external executables, and resolution of the noted issue below.

## Outputs

With `BASE_DIR = "test_april9v1"`, the simulation is intended to create a directory such as:

```text
test_april9v1/
└── sim_001/
    ├── sim_summary.tsv
    ├── gemma_summary.tsv
    └── false_positive_summary.tsv
```

- **`sim_summary.tsv`**: A per-simulation record of causal allele frequencies, mutation/site counts, nucleotide diversity, divergence, and pairwise FST estimates for the simulated populations.
- **`gemma_summary.tsv`**: Association statistics for the causal position, including allele frequency, effect estimate, standard error, Wald-test p-value, and sample-size information.
- **`false_positive_summary.tsv`**: Counts of genome-wide significant sites, exact causal-position hits (`tp_abs`), other significant positions (`fp_abs`), hits within the causal window (`tp_rel`), and significant sites outside that window (`fp_rel`).

**Interpretation:** In this code, `fp_abs` means *significant sites other than the exact causal position*. It may include variants linked to the causal locus. `fp_rel` excludes significant positions inside the specified causal window and is therefore a different, region-based measure. Neither number alone establishes the statistical type-I error rate without further analysis.

Temporary per-sample working directories and the forward tree-sequence output are removed by the notebook after processing. The supplied notebook does not currently generate plots.

## Research context

Developed by **Diego Salazar** as part of computational research into GWAS statistical power and sample-size requirements in evolutionary genomics at the University of Oslo.
