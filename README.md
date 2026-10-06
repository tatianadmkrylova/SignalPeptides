# SignalPeptides

This repository contains the datasets, scripts, SignalP 6.0 predictions, and analyses used to study signal peptides in viral, bacterial, and mammalian proteins.

The project combines data from:

- UniProt;
- signalpeptide.de;
- the SignalP training set;
- SignalP 6.0 predictions generated from the prepared protein datasets.

The repository is organized to preserve the datasets used in the analyses and to make the main processing and analysis steps reproducible.

For more details:

- repository structure: `docs/project_structure.md`
- script execution and arguments: `docs/run_scripts.md`

## Workflow

The project consists of four main stages:

1. Dataset acquisition
2. UniProt metadata retrieval
3. FASTA preparation and SignalP 6.0 prediction
4. Statistical analysis and visualization

## 1. Dataset acquisition

The initial datasets are obtained from UniProt and signalpeptide.de.

### UniProt

UniProt datasets are retrieved using:

```text
scripts/01_data_acquisition/uniprot_dataset.sh
```

The viral negative dataset containing proteins without annotated signal peptides is retrieved using:

```text
scripts/01_data_acquisition/nosignalpeptide_data.sh
```

Raw UniProt datasets are stored in:

```text
data/raw/uniprot/
```

### signalpeptide.de

Datasets from signalpeptide.de are retrieved using:

```text
scripts/01_data_acquisition/signalpeptide.py
```

Raw signalpeptide.de datasets are stored in:

```text
data/raw/signalpeptide_de/
```

## 2. UniProt metadata retrieval

Additional protein metadata are retrieved from the UniProt REST API.

For UniProt-derived datasets:

```text
scripts/02_metadata/uniprot_request_metadata.py
```

For proteins obtained from signalpeptide.de:

```text
scripts/02_metadata/uniprot_request_metadata_signalpeptide.py
```

Processed metadata are stored in:

```text
data/processed/
```

The processed data are separated according to their source:

```text
data/processed/uniprot/
data/processed/signalpeptide_de/
data/processed/signalp_training_set/
```

## 3. FASTA preparation and SignalP 6.0 prediction

Full protein sequences are retrieved from UniProt using:

```text
scripts/03_prepare_signalp/seq_full_prot_fasta_sec_acc.py
```

FASTA files can then be processed to remove duplicate protein sequences using:

```text
scripts/03_prepare_signalp/fasta_for_signalp.sh
```

Intermediate FASTA files are stored in:

```text
data/interim/fasta/
```

These generated FASTA files are not tracked by Git.

### SignalP 6.0

SignalP 6.0 is run separately using the prepared FASTA files.

Example:

```bash
signalp6 --fastafile /path/to/input.fasta --organism other --output_dir /path/to/output --format txt --mode slow-sequential
```

For viral proteins, the SignalP organism category used in this project is:

```text
other
```

SignalP output files used in the analyses are stored in:

```text
results/signalp/
```

They are organized by organism and data source, including:

```text
results/signalp/viruses_uniprot/
results/signalp/viruses_signalpeptide_de/
results/signalp/bacteria_uniprot/
results/signalp/bacteria_signalpeptide_de/
results/signalp/mammalia_uniprot/
results/signalp/mammalia_signalpeptide_de/
results/signalp/negative_dataset/
```

## 4. Statistical analysis and visualization

The R scripts used for the analyses are located in:

```text
scripts/04_analysis/
```

### Dataset comparison

```text
scripts/04_analysis/Datasets.R
```

This script compares UniProt, signalpeptide.de, and the SignalP training set.

The analyses include:

- intersections between UniProt and signalpeptide.de datasets;
- Venn diagrams;
- comparison of signal peptide lengths;
- comparison of signal peptide length distributions;
- analysis of UniProt ECO evidence codes;
- comparison with the SignalP training set;
- comparison of UniProt signal peptide lengths between viruses, bacteria, and mammals.

### SignalP comparison across datasets

```text
scripts/04_analysis/SignalP_UniProt_signalpeptide_analysis.R
```

This script compares SignalP 6.0 predictions with signal peptide annotations from UniProt and signalpeptide.de for viruses, bacteria, and mammals.

The analyses include:

- comparison of SignalP and UniProt signal peptide lengths;
- comparison of SignalP and signalpeptide.de signal peptide lengths;
- comparison between UniProt and signalpeptide.de;
- evaluation of SignalP performance on the viral dataset using a negative dataset;
- calculation of classification metrics and a confusion matrix.

### Detailed viral SignalP analysis

```text
scripts/04_analysis/SignalP_virus_UniProt_prediction_analysis.R
```

This script performs a more detailed analysis of SignalP 6.0 predictions for viral proteins with signal peptides annotated in UniProt.

The analyses include:

- proteins classified as `OTHER` by SignalP 6.0;
- SignalP class probabilities for non-predicted proteins;
- comparison of UniProt signal peptide lengths between predicted and non-predicted proteins;
- identification of high-confidence SignalP predictions.

Generated figures are stored in:

```text
results/figures/
```

## Data organization

The main data directories are:

```text
data/raw/        Original downloaded datasets
data/external/   Data obtained from external resources
data/interim/    Generated intermediate files
data/processed/  Metadata and processed datasets
```

SignalP results and generated figures are stored separately under:

```text
results/
```

## Data version note

The UniProt datasets used for the analyses in this project correspond to versions downloaded **before the UniProt update of 18 July 2026**.

The repository preserves these datasets so that the analyses and figures correspond to the data used during the project.

Raw source datasets are stored in `data/raw/`, while metadata and derived tables generated by the project scripts are stored in `data/processed/`.

## Running the scripts

Detailed commands, required arguments, and expected output locations are documented in:

```text
docs/run_scripts.md
```

Scripts using relative paths should be executed according to the working-directory instructions provided in that document.