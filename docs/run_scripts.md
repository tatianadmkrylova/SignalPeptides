# Running the scripts

This document describes how to run the scripts included in the project, their required arguments, working directories, and expected output locations.

The commands below assume that the repository has been cloned locally.

## Important note about script execution

The Python and Bash scripts contain a shebang line specifying the appropriate interpreter. Therefore, executable scripts can be run directly using `./` without explicitly typing `python3` or `bash`.

For example:

```bash
./signalpeptide.py viruses viruses
```

instead of:

```bash
python3 signalpeptide.py viruses viruses
```

and:

```bash
./uniprot_dataset.sh
```

instead of:

```bash
bash uniprot_dataset.sh
```

The commands below therefore use the `./script_name` syntax.

## Important note about working directories

Some scripts use relative paths to read or save files. These scripts must therefore be executed from the directory in which they are located.

For example, scripts in:

```text
scripts/01_data_acquisition/
```

should be run after moving to this directory:

```bash
cd scripts/01_data_acquisition
```

Similarly, metadata scripts should be run from:

```text
scripts/02_metadata/
```

and FASTA preparation scripts from:

```text
scripts/03_prepare_signalp/
```

The R analysis scripts use paths relative to the repository root and should therefore be executed from the root directory of the repository.

---

## 1. Data acquisition

The scripts in `scripts/01_data_acquisition/` retrieve the initial datasets from UniProt or signalpeptide.de.

Move to the corresponding directory:

```bash
cd scripts/01_data_acquisition
```

### 1.1. UniProt datasets

Script:

```text
uniprot_dataset.sh
```

This script retrieves UniProt datasets containing proteins with signal peptide annotations for viruses, mammals, and bacteria according to the queries defined in the script.

Run:

```bash
./uniprot_dataset.sh
```

No arguments are required.

The downloaded tables are automatically saved in:

```text
data/raw/uniprot/
```

The UniProt API returns tab-separated data. The script converts tab separators to commas and saves the resulting tables as CSV files.

---

### 1.2. Viral proteins without signal peptide annotations

Script:

```text
nosignalpeptide_data.sh
```

This script retrieves viral proteins without a UniProt signal peptide annotation according to the query defined in the script.

This dataset is used to construct the viral negative dataset for SignalP evaluation.

Run:

```bash
./nosignalpeptide_data.sh
```

No arguments are required.

The resulting CSV file is automatically saved in:

```text
data/raw/uniprot/
```

---

### 1.3. signalpeptide.de datasets

Script:

```text
signalpeptide.py
```

This script retrieves signal peptide records from signalpeptide.de for a selected organism group.

Run:

```bash
./signalpeptide.py <output_name> <organism_class>
```

Arguments:

- `<output_name>` — name used to construct the output file name;
- `<organism_class>` — organism class to retrieve, for example `viruses`, `mammalia`, or `bacteria`.

Example:

```bash
./signalpeptide.py viruses viruses
```

The output file is automatically saved in:

```text
data/raw/signalpeptide_de/
```

with a name of the form:

```text
signalpeptide_<output_name>.csv
```

For example:

```text
signalpeptide_viruses.csv
```

---

## 2. UniProt metadata retrieval

The scripts in `scripts/02_metadata/` retrieve additional protein information from UniProt using accession numbers contained in previously generated datasets.

Move to the corresponding directory:

```bash
cd scripts/02_metadata
```

### 2.1. Metadata for UniProt datasets

Script:

```text
uniprot_request_metadata.py
```

Run:

```bash
./uniprot_request_metadata.py <input.csv> <output.csv>
```

Arguments:

- `<input.csv>` — input CSV file containing UniProt accession numbers in the `Entry` column;
- `<output.csv>` — name of the metadata output file.

Example:

```bash
./uniprot_request_metadata.py ../../data/raw/uniprot/uniprotkb_taxonomy_viruses_AND_ft_sign_entry.csv uniprotkb_taxonomy_viruses_AND_ft_sign_metadata.csv
```

The output file is automatically saved in:

```text
data/processed/uniprot/
```

The script retrieves information including:

- UniProt primary accession;
- protein length;
- organism information;
- host information;
- signal peptide annotations;
- transmembrane annotations;
- annotation evidence when available.

---

### 2.2. Metadata for signalpeptide.de datasets

Script:

```text
uniprot_request_metadata_signalpeptide.py
```

Run:

```bash
./uniprot_request_metadata_signalpeptide.py <input.csv> <output.csv>
```

Arguments:

- `<input.csv>` — signalpeptide.de CSV file containing protein accession numbers;
- `<output.csv>` — name of the metadata output file.

Example:

```bash
./uniprot_request_metadata_signalpeptide.py ../../data/raw/signalpeptide_de/signalpeptide_viruses.csv signalpeptide_viruses_all_metadata.csv
```

The output file is automatically saved in:

```text
data/processed/signalpeptide_de/
```

---

## 3. FASTA preparation for SignalP

The scripts in `scripts/03_prepare_signalp/` prepare protein sequences for SignalP 6.0 analysis.

Move to the corresponding directory:

```bash
cd scripts/03_prepare_signalp
```

### 3.1. Retrieve full protein sequences from UniProt

Script:

```text
seq_full_prot_fasta_sec_acc.py
```

This script reads UniProt accession numbers from a metadata CSV file and retrieves the corresponding full protein sequences in FASTA format.

Run:

```bash
./seq_full_prot_fasta_sec_acc.py <input_metadata.csv> <output.fasta>
```

Arguments:

- `<input_metadata.csv>` — metadata CSV file containing a `PrimaryAccession` column;
- `<output.fasta>` — FASTA file in which the retrieved protein sequences will be saved.

Example:

```bash
./seq_full_prot_fasta_sec_acc.py ../../data/processed/uniprot/uniprotkb_taxonomy_viruses_AND_ft_sign_before_update_metadata.csv ../../data/interim/fasta/viruses_uniprot.fasta
```

Generated FASTA files should be stored in:

```text
data/interim/fasta/
```

These files are intermediate files and are not tracked by Git.

---

### 3.2. Remove duplicate protein sequences

Script:

```text
fasta_for_signalp.sh
```

This script converts a FASTA file to a table, removes duplicate protein sequences, and creates a FASTA file containing unique sequences.

Run:

```bash
./fasta_for_signalp.sh <input.fasta> <output_all.csv> <output_unique.csv> <output_unique.fasta>
```

Arguments:

- `<input.fasta>` — original FASTA file;
- `<output_all.csv>` — CSV table containing all protein identifiers and sequences;
- `<output_unique.csv>` — CSV table containing only unique sequences;
- `<output_unique.fasta>` — FASTA file containing only unique sequences.

Example:

```bash
./fasta_for_signalp.sh ../../data/interim/fasta/viruses_uniprot.fasta ../../data/interim/fasta/viruses_uniprot_all.csv ../../data/interim/fasta/viruses_uniprot_unique.csv ../../data/interim/fasta/viruses_uniprot_unique.fasta
```

The unique FASTA file can then be used as input for SignalP 6.0.

---

## 4. Running SignalP 6.0

SignalP 6.0 itself is not included in this repository and must be installed separately.

A prepared FASTA file can be analysed using a command such as:

```bash
signalp6 --fastafile /path/to/input.fasta --organism other --output_dir /path/to/output --format txt --mode slow-sequential
```

For the viral datasets analysed in this project, the SignalP organism category used is:

```text
other
```

SignalP outputs used in the project are stored in:

```text
results/signalp/
```

They are organized according to organism group and source dataset:

```text
results/signalp/
├── bacteria_signalpeptide_de/
├── bacteria_uniprot/
├── mammalia_signalpeptide_de/
├── mammalia_uniprot/
├── negative_dataset/
├── viruses_signalpeptide_de/
└── viruses_uniprot/
```

Each directory contains SignalP prediction results used by the downstream R analysis scripts.

---

## 5. Statistical analysis and visualization

The R scripts are located in:

```text
scripts/04_analysis/
```

Unlike the Python and Bash scripts described above, these analysis scripts use paths relative to the repository root.

Therefore, return to the root directory of the repository before running them.

For example, if the current directory is `scripts/03_prepare_signalp/`:

```bash
cd ../..
```

The scripts can then be executed with `Rscript`.

### 5.1. Dataset comparison

Script:

```text
scripts/04_analysis/Datasets.R
```

Run from the repository root:

```bash
Rscript scripts/04_analysis/Datasets.R
```

This script analyses and compares datasets from UniProt, signalpeptide.de, and the SignalP training set.

The analyses include:

- intersections between UniProt and signalpeptide.de datasets;
- Venn diagrams;
- signal peptide length distributions;
- comparison of signal peptide lengths between UniProt and signalpeptide.de;
- UniProt ECO evidence distributions;
- comparison with the SignalP training set;
- comparison of signal peptide length distributions between bacterial, mammalian, and viral UniProt proteins.

Generated figures are saved in:

```text
results/figures/
```

---

### 5.2. SignalP comparison with UniProt and signalpeptide.de

Script:

```text
scripts/04_analysis/SignalP_UniProt_signalpeptide_analysis.R
```

Run from the repository root:

```bash
Rscript scripts/04_analysis/SignalP_UniProt_signalpeptide_analysis.R
```

This script compares SignalP 6.0 predictions with signal peptide annotations from UniProt and signalpeptide.de for viruses, bacteria, and mammals.

The analyses include:

- comparison of SignalP and UniProt signal peptide lengths;
- comparison of SignalP and signalpeptide.de signal peptide lengths;
- comparison of UniProt and signalpeptide.de annotations;
- evaluation of SignalP performance for the viral dataset;
- use of the viral negative dataset;
- calculation of classification metrics;
- generation of a confusion matrix.

Generated figures are saved in:

```text
results/figures/
```

The script reads the required SignalP prediction files from:

```text
results/signalp/
```

---

### 5.3. Detailed viral SignalP prediction analysis

Script:

```text
scripts/04_analysis/SignalP_virus_UniProt_prediction_analysis.R
```

Run from the repository root:

```bash
Rscript scripts/04_analysis/SignalP_virus_UniProt_prediction_analysis.R
```

This script performs a detailed analysis of SignalP 6.0 predictions for viral proteins with signal peptides annotated in UniProt.

The analyses include:

- analysis of proteins classified as `OTHER` by SignalP;
- analysis of SignalP prediction probabilities;
- comparison of signal peptide lengths between predicted and non-predicted proteins;
- identification of high-confidence SignalP predictions.

Generated figures are saved in:

```text
results/figures/
```

---

## 6. Generated files

The main generated files are organized as follows:

```text
data/interim/fasta/
```

Intermediate FASTA and sequence-processing files. These files are not tracked by Git.

```text
data/processed/
```

Metadata and processed datasets generated by the project scripts.

```text
results/signalp/
```

SignalP 6.0 prediction outputs used in the analyses.

```text
results/figures/
```

Figures generated by the R analysis scripts.

---

## Data version note

The analyses in this repository use UniProt datasets downloaded before the UniProt update of 18 July 2026.

These versions are intentionally retained so that the source datasets, SignalP outputs, statistical analyses, and figures remain consistent with the data used during the project.
