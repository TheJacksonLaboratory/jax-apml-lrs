# Troubleshooting

Common issues and solutions for running `jax-apml-lrs` on Sumner2.

---

## 1. Pipeline must be submitted via `sbatch`

**Symptom:** Pipeline fails immediately or `outputDir` cannot be written to.

**Cause:** `/flashscratch` is only mounted on compute nodes, not all login nodes. Running the pipeline interactively from a login node will fail.

**Solution:** Always submit the pipeline as a SLURM job using `sbatch`. See the README for the full submission command.

---

## 2. `apptainer: command not found`

**Symptom:** A process fails with:
```
bash: line 1: apptainer: command not found
```

**Cause:** Apptainer is not in the default PATH on compute nodes.

**Solution:** Always include the following in your `--wrap` command before calling `run.sh`:
```bash
export PATH=/cm/local/apps/apptainer/current/bin:$PATH
```

This is already included in the submission commands documented in the README.

---

## 3. Transient container pull failures

**Symptom:** A process fails with:
```
Failed to pull singularity image
FATAL: While making image from oci registry: ...
```
with a status of 255.

**Cause:** Transient network or registry errors when pulling container images at runtime.

**Solution:** Resubmit with `--resume` — Nextflow will retry the failed pull and skip any already-completed steps:
```bash
./run.sh -w <workflow> -p hpc \
    --csv_path /path/to/samplesheet.csv \
    --outputDir /flashscratch/${USER}/jax-apml-lrs/ \
    --refs_path /path/to/refs/ \
    --resume
```

If the failure persists across multiple retries, check network connectivity from the compute nodes or verify that the container tag exists on the registry.

---

## 4. pbmm2 out-of-memory (OOM)

**Symptom:** The `PBMM2_ALIGN` process is killed with exit code 137, or you receive an email from the cluster stating the job was terminated for exceeding memory limits.

**Cause:** pbmm2 calculates sort memory based on total node RAM divided by sort threads, which can exceed the SLURM allocation.

**Solution:** Reduce `sort_threads` and `sort_memory` in `workflows/lrs_read/nextflow.config`:
```
pbmm2 {
    sort_threads = 4
    sort_memory  = "8G"
}
```

This limits the total sort memory footprint to ~32GB, well within a typical SLURM allocation.

---

## 5. hifiasm out-of-memory (OOM)

**Symptom:** The `ASSEMBLY_HIFIASM` or `ASSEMBLY_HIFIASM_TRIO` process is killed with exit code 137.

**Cause:** hifiasm is memory-intensive, especially for whole-genome assemblies. The default process memory may be insufficient.

**Solution:** The `withName` block in `workflows/lrs_asm_single/nextflow.config` and `workflows/lrs_asm_trio/nextflow.config` sets memory to 200GB by default. If this is still insufficient for your sample, increase it:
```
withName: 'ASSEMBLY_HIFIASM' {
    memory = '300 GB'
    cpus   = 96
}
```

---

## 6. `java.io.IOException: No locks available`

**Symptom:** A workflow fails immediately with:
```
java.io.IOException: No locks available
```

**Cause:** Two instances of the same workflow are running simultaneously. Each workflow (e.g. `lrs_read`) shares a Nextflow session cache, and concurrent runs conflict on the cache lock file.

**Solution:** Do not submit the same workflow more than once at a time. Wait for the first run to complete before resubmitting. If the error occurs after an abrupt job cancellation, clear the stale lock and resume:
```bash
find /path/to/jax-apml-lrs/.nextflow/cache -name "LOCK" -delete
./run.sh -w <workflow> -p hpc \
    --csv_path /path/to/samplesheet.csv \
    --outputDir /flashscratch/${USER}/jax-apml-lrs/ \
    --refs_path /path/to/refs/ \
    --resume
```

---

## 7. `download_refs.sh` — container runtime not found

**Symptom:** Running `download_refs.sh` fails with:
```
ERROR: 'apptainer' not found.
```

**Cause:** Apptainer is not in the default PATH on login nodes.

**Solution:** The script automatically resolves the full apptainer path at `/cm/local/apps/apptainer/current/bin/apptainer`. If this path has changed on your system, specify the runtime explicitly:
```bash
./download_refs.sh \
    --container-runtime /path/to/apptainer \
    --outdir /path/to/refs/
```

---

## 8. SvAnna out-of-memory (OOM)

**Symptom:** `PRIORITIZE_BY_SVANNA_UNFILTERED` fails with exit code 140.

**Cause:** SvAnna memory requirements scale with the number of variants in the input VCF. The unfiltered PAV SV VCF contains more variants than the RARE+UNIQUE filtered VCF, requiring more memory.

**Solution:** Increase the memory allocation for `PRIORITIZE_BY_SVANNA_UNFILTERED` in `workflows/lrs_asm_single/subworkflow/05_PRIORITIZE_SV_SVANNA.nf` and/or `workflows/lrs_asm_trio/subworkflow/05_PRIORITIZE_SV_SVANNA.nf`:
```
process PRIORITIZE_BY_SVANNA_UNFILTERED {
    memory '400 GB'
    ...
}
```

---

## 9. Running multiple samples

To process multiple samples, add all samples to a single samplesheet and submit one workflow run. Nextflow will parallelize across samples automatically.

**Do not** submit separate workflow runs for each sample simultaneously — each workflow (e.g. `lrs_read`) shares a Nextflow session cache, and concurrent runs of the same workflow will conflict with a `No locks available` error.

If you must run samples separately, wait for one run to fully complete before submitting the next, or use `--resume` to continue a failed run.
