# Python 3.12 upgrade plan

Source: assembler audit, 2026-10-12 (`audit/assembler_audit.html`).
Completed items are struck through with a note on where the work landed. Nothing is deleted.

## Recommended next steps

1. ~~Target Python 3.12, not 3.10. Python 3.10 reaches end of life this month (October 2026); 3.12 is supported until October 2028.~~
   Done 2026-10-12: decided; implemented by step 3.
2. Fix compleasm with a one-line pull request to bioconda-recipes that removes `dendropy <4.6.0` from recipes/compleasm/meta.yaml. Until it merges, install compleasm with `--no-deps` after the environment is built, the same way the Dockerfile installs tiara-entheome. Its real dependencies (hmmer, miniprot, sepp, pandas) are already in the environment.
   Progress 2026-10-12: the `--no-deps` workaround is in the Dockerfile, `entheome.sif.def` and CI on `chore/py312`. The same fix already exists as the stale bioconda-recipes PR [#60482](https://github.com/bioconda/bioconda-recipes/pull/60482). We commented with the sepp/dendropy evidence ([comment](https://github.com/bioconda/bioconda-recipes/pull/60482#issuecomment-6134894597)) and are waiting for it to be rebased and merged. Strike this step when that PR merges. **This blocks the bioconda `egap` release on 3.12** (meta.yaml still lists compleasm).
3. ~~Create the branch `chore/py312` and change the Python pin everywhere it appears: `meta.yaml`, the `Dockerfile`, `entheome.sif.def`, `EGAP_setup.sh`, and the CI test matrix (3.12, plus 3.13 to catch problems early).~~
   * ~~Add `pandas>=2.2,<3` and `numpy<3`. Left unpinned at 3.12, conda picks pandas 3.0, a major release whose changes in behaviour (copy-on-write, a new string type) could quietly alter EGAP's tables.~~
   * Move pandas 3 in a later, separate round. *(still open)*

   Done 2026-10-12 on `chore/py312` (branched from origin/main; uncommitted). Changes:
   * Python `>=3.12,<3.13` in all five places. CI runs 3.12 and 3.13, plus a floor job (pandas 2.2.0, numpy 1.26.4). `tests/requirements.txt` matches.
   * `numpy>=1.26,<3`, because 1.26 is the first numpy with 3.12 builds. `busco=5.8.3` pinned to the pre-3.12 version so Round 1 changes no tools.
   * compleasm is installed with `--no-deps` (step 2 workaround), with `hmmer` and `miniprot` added explicitly. The CI tools job drops upstream `tiara`, which needs numpy<1.20.
   * **Found and fixed: channel order.** The Dockerfile, `entheome.sif.def` and `EGAP_setup.sh` listed `-c bioconda` first. At 3.12, bioconda's stale biopython/tqdm/openpyxl copies win. The solve either fails (strict priority) or doesn't finish in 20 minutes (flexible priority, which is the Docker base image's default). All three now use `--strict-channel-priority -c conda-forge -c bioconda`, Bioconda's documented setup.
   * **Found and fixed: `runner` was never in the Docker image.** It was installed with `setup.py install --user` into `/root/.local`, outside the env that conda-pack ships. It now goes in with `pip install ./runner`.

   Verified:
   * The Dockerfile's exact EGAP_env list solves at 3.12 with every tool version unchanged (542 packages; masurca 4.1.4, spades 4.2.0, flye 2.9.5, hifiasm 0.25.0, quast 5.3.0, busco 5.8.3, kraken2 2.1.6, sra-tools 3.2.0, trimmomatic 0.40, nanoplot 1.46.2, bbmap 39.52, racon 1.5.0, pilon 1.24, ragtag 2.1.0; sepp 4.5.6 with dendropy 5.1.0).
   * Unit tests on Python 3.12.15, pandas 2.3.3, numpy 2.5.4: 132 passed, 14 skipped (fastb/pigz binaries not present).

   Not yet verified: an actual image build. That happens in step 4.
4. Round 1 A/B test: change only Python and the libraries, keeping every tool version the same. Assemblies should come out identical; if they don't, the change of Python is the cause. Run the existing E. coli tests plus one same-strain set, for example K. pneumoniae.
5. Round 2: raise the tool versions to SPAdes 4.3.0, Flye 2.9.6, compleasm 0.2.10 and BUSCO 6.1. Score both images with the same BUSCO version so the measuring stick stays fixed.

## MaSuRCA: patch the package, then fork if needed

The real MaSuRCA bugs are the ones EGAP currently works around:

* The Flye collision. MaSuRCA's package ships its own `bin/flye`, which starts Python 2.7 and brings minimap2 2.17. Those files overwrite the real Flye; the fix is the Dockerfile repair at lines 118–144.
* A hard-coded Flye path. MaSuRCA looks for Flye at `$MASURCA/../Flye/bin/flye`. EGAP rewrites its generated scripts with regular expressions to get past this (`assemble_masurca.py:222-370`).
* Calls to `python2.7`. EGAP patches these to point at Python 3.
* A leftover overlap store after an interrupted run (`*.ovlStore.BUILDING`) blocks every rerun until it's deleted. `clear_stale_overlap_store` does that cleanup.

The cheapest fix is a pull request to MaSuRCA's Bioconda recipe so it stops installing `bin/flye` and `flye-minimap2` into the shared `bin/` folder. That fixes the problem for every MaSuRCA user, not just EGAP. Send the Flye-path and Python 3 fixes upstream as well.
The upstream repo has had no commits for about 16 months. If it doesn't respond, fork it as `masurca-entheome`, the same way you forked tiara-entheome. With the package fixed, roughly 150 lines of workarounds in `assemble_masurca.py` and the Dockerfile repair could be deleted.

## Speed and memory gains to A/B test after the move

* SPAdes 4.3: uses less memory, only use `--frugal` when CPU threads is less than 8.
* MaSuRCA: EGAP currently forces its slower internal assembler (`flye assembly: 0`). Setting `FLYE_ASSEMBLY=1` makes MaSuRCA use Flye for the final assembly step, which is typically much faster and uses less memory. That's worth an A/B test.
* Pilon is a Java tool that needs a lot of memory. Polypolish or Pypolca do the same job with far less, and the fungal and bacterial benchmarks in the audit favour them.
* Racon: no release in nearly five years. `dorado polish` replaces it for nanopore data, and Medaka 2.2 is the other current option. Both need Python 3.12, so after this move they could run in the main environment.

## Side findings (not steps yet)

* The Dockerfile/`entheome.sif.def` comment claiming "Flye >=3.0 changes assembly graph format" was false (Flye 3 does not exist). Corrected on `chore/py312`.
* `EGEP_env` (Dockerfile and `entheome.sif.def`) already asks for `python>=3.10` with `compleasm>=0.2.8`, so it hits the same dendropy conflict as step 2.
