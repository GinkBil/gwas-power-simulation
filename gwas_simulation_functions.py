#!/usr/bin/env python
# coding: utf-8

# In[1]:

import os
import csv
import shutil
import subprocess
from pathlib import Path

import numpy as np
import tskit
import pyslim
import msprime
import pandas as pd

def build_params_flat_list(params_dict):
    params = []
    for p, v in params_dict.items():
        if isinstance(v, (list, tuple)):
            v = "c(" + ",".join(map(str, v)) + ")"
        params.extend(["-d", f"{p}={v}"])
    return params
    
def run_slim_simulation(SLIM_BIN, slim_script, params, sim_dir):
    subprocess.run(
        [SLIM_BIN, *params, str(Path(slim_script).resolve())],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
        cwd=sim_dir
    )
    tree_file = sim_dir/"forward.trees"
    return tree_file

def recapitate_ts(ts, RECAP_NE, RECAP_RECOMB_RATE, RNG_SEED, sim_iteration):
    ts = pyslim.recapitate(
        ts,
        ancestral_Ne=RECAP_NE,
        recombination_rate=RECAP_RECOMB_RATE,
        random_seed=RNG_SEED + 100000 + sim_iteration
    )
    return ts

def overlay_neutral_mutations(ts, RECAP_MUT_RATE, RNG_SEED, sim_iteration):
    next_id = pyslim.next_slim_mutation_id(ts)

    ts = msprime.sim_mutations(
        ts,
        rate=RECAP_MUT_RATE,
        random_seed=RNG_SEED + 100000 + sim_iteration,
        keep=True,
        model=msprime.SLiMMutationModel(type=0, next_id=next_id)
    )
    return ts

def set_fixed_ancestral_state(ts, value: str = "Anc"):
    """
    Returns a new TreeSequence where every site's ancestral_state is set to 'value'.
    Does NOT modify mutations/derived states.
    """
    tables = ts.dump_tables()
    n = tables.sites.num_rows
    tables.sites.packset_ancestral_state([value] * n)
    return tables.tree_sequence()

def get_m2_positions(ts, mutation_type):
    mutations_type2 = []
    for m in ts.mutations():
        if m.metadata['mutation_list'][0]['mutation_type'] == mutation_type:
            site_id = m.site
            pos = int(ts.site(site_id).position)
            mutations_type2.append(pos)
    return mutations_type2

# Sort ind. IDs by population
def sort_inds_by_pop(ts, inds_sample):
    pop_inds = {}
    for pop in range(0,ts.num_populations-1):
        pop_inds[pop] = [
            ind_id for ind_id in inds_sample
            if ts.node(ts.individual(ind_id).nodes[0]).population == pop
        ]
    return pop_inds

# Get nodes IDs using ind. IDs
def get_sample_nodes(ts, sampled_inds):
    sample_nodes = {}
    for pop in range(0, ts.num_populations-1):
        sample_nodes[pop] = [node 
        for ind_id in sampled_inds[pop]
        for node in ts.individual(ind_id).nodes]
    return sample_nodes

def get_site_frequencies(ts, m2_positions, sample_nodes):
    site_id = np.searchsorted(ts.sites_position, m2_positions[0])
    freqs = []
    # Assumes 1 extra population after recapitation, and pop0 is pop1 after simplification.
    for pop in range(0,ts.num_populations-1):
        samples = sample_nodes[pop]
        if len(samples) == 0:
            freqs.append(np.float64(0))
            continue
        v = tskit.Variant(ts, samples=samples)
        v.decode(site_id)
        freqs.append(v.frequencies()[v.alleles[1]])
    return {
        "freq_pop1": freqs[0],
        "freq_pop2": freqs[1],
        "freq_pop3": freqs[2],
    }

def compute_ts_population_descriptors(ts):

    s1 = ts.samples(population=0)
    s2 = ts.samples(population=1)
    s3 = ts.samples(population=2)

    sets=[s1, s2, s3]
    idxs=[(0, 1), (0, 2), (1, 2)]

    pi = ts.diversity(sample_sets=sets)
    dxy = ts.divergence(sample_sets=sets,indexes=idxs)
    fst = ts.Fst(sample_sets=sets,indexes=idxs)
    dxy_b = ts.divergence(sample_sets=sets,indexes=idxs,mode="branch")
    return {
        "p1_diversity": pi[0],
        "p2_diversity": pi[1],
        "p3_diversity": pi[2],
        "p1p2_divergence": dxy[0],
        "p1p3_divergence": dxy[1],
        "p2p3_divergence": dxy[2],
        "p1p2_fst": fst[0],
        "p1p3_fst": fst[1],
        "p2p3_fst": fst[2],
        "p1p2_dxy_branch": dxy_b[0],
        "p1p3_dxy_branch": dxy_b[1],
        "p2p3_dxy_branch": dxy_b[2],
    }

def write_sim_summary_row(ts, m2_positions, sample_inds, sim_iteration, sim_dir):

    sample_nodes = get_sample_nodes(ts, sample_inds)
    freqs = get_site_frequencies(ts, m2_positions, sample_nodes)
    stats = compute_ts_population_descriptors(ts)
    row = {
        "sim": sim_iteration,
        **freqs,
        "sim_mutations": ts.num_mutations,
        "sim_sites": ts.num_sites,
        **stats,
    }
    headers = row.keys()
    with open(sim_dir / "sim_summary.tsv", 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

def sample_individuals(ts, pop_inds, sample_sizes):
    sampled_inds = {}
    # Sampling individuals in each population
    for pop in range(0, ts.num_populations-1):
        sampled_inds[pop] = np.random.choice(pop_inds[pop], 
                                             size=sample_sizes[pop],
                                             replace=False)
    return sampled_inds
    
def get_genotypes(ts, sample_nodes, m2_positions):
    site_id = np.searchsorted(ts.sites_position, m2_positions[0])
    genotypes = {}

    for pop in range(0, ts.num_populations-1):
        v = tskit.Variant(ts, samples=sample_nodes[pop])
        v.decode(site_id)
        genotypes[pop] = [sum(v.genotypes[i:i+2]) 
            for i in range(0, len(v.genotypes), 2)]
    return genotypes

def get_phenotypes(ts, genotypes, beta, noise_sd):
    phenotypes = {}
    for pop in range(0, ts.num_populations-1):
        phenotypes[pop] = [genotype*beta for genotype in genotypes[pop]]
        noise = np.random.normal(0, noise_sd, len(phenotypes[pop]))
        phenotypes[pop] = phenotypes[pop] + noise
    return phenotypes

def write_sample_vcf(ts_vcf, sampled_inds, path):
    inds_list = [i for pop in sampled_inds for i in sampled_inds[pop]]
    names = [f'ind_{i}' for i in inds_list]
    with open(path, "w") as f:
        ts_vcf.write_vcf(f, individuals=inds_list, individual_names=names)

def convert_vcf_to_bed(PLINK_BIN, sim_dir: Path):
    subprocess.run(
        [
            PLINK_BIN,
            "--vcf", "geno.vcf",
            "--double-id",
            "--allow-no-sex",
            "--set-missing-var-ids", "@:#",
            "--make-bed",
            "--out", "data"
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
        cwd=sim_dir
    )

def write_bedformat_phenotypes(sim_dir, phenotypes):
    df = pd.read_csv(sim_dir / "data.fam", sep=" ", header=None)
    pheno_list = [i for pop in phenotypes for i in phenotypes[pop]]
    df[5] = pheno_list
    df.to_csv(sim_dir / "data.fam", sep=" ", header=False, index=False)

def run_gemma_kinship(GEMMA_BIN, sim_dir: Path):
    output_dir = sim_dir / "output"
    output_dir.mkdir(exist_ok=True)

    subprocess.run(
        [
            GEMMA_BIN,
            "-bfile", "data",
            "-gk",
            "-o", "kinship"
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
        cwd=sim_dir
    )

def run_gemma_lmm(GEMMA_BIN, sim_dir: Path):
    subprocess.run(
        [
            GEMMA_BIN,
            "-bfile", "data",
            "-k", "output/kinship.cXX.txt",
            "-lmm",
            "-o", "gwas"
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
        cwd=sim_dir
    )

def get_sample_m2_assoc(gwas_results, m2_positions, sample_size, sim_iteration):
    df = pd.read_csv(
        gwas_results, sep="\t",
        usecols=["chr", "ps", "allele1", "allele0", "af", "beta", "se", "p_wald"])
    df_m2 = df.loc[df['ps'].isin(m2_positions)].copy()
    for pop, n_sample in enumerate(sample_size, start=1):
        df_m2[f"sample_p{pop}"] = n_sample
    df_m2['sim'] = sim_iteration
    return df_m2

def false_positive_test(gwas_results, causal_ps, window_size, sample_size, sim_iteration):
    df = pd.read_csv(
        gwas_results, sep="\t",
        usecols=["ps", "p_wald"])

    significative = df["p_wald"] < 5e-8
    causal = df["ps"].isin(causal_ps)
    sig_snps = (significative).sum()
    tp_abs = (significative & causal).sum()
    fp_abs = sig_snps - tp_abs
    tp_rel = (
        df.loc[significative, "ps"]
          .apply(lambda pos: any(abs(pos - c) <= window_size for c in causal_ps))
          .sum()
    )
    fp_rel = sig_snps - tp_rel
    return {"sim": sim_iteration,
            "sample_size": sample_size, 
            "sig_snps": sig_snps,
            "tp_abs": tp_abs,
            "fp_abs": fp_abs, 
            "tp_rel": tp_rel,
            "fp_rel": fp_rel}

def calculate_var_genotypes(genotypes, pop=0):
    X = np.asarray(genotypes[pop], dtype=float)
    varX = np.var(X, axis=0, ddof=0)
    return varX

def calculate_beta(h2, varX):
    noise_sd = 1.0
    h2 = np.array(h2)

    beta = np.sqrt((h2 * noise_sd**2) / ((1 - np.sum(h2)) * varX))

    return beta, noise_sd
    

# In[ ]:




