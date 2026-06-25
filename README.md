# SignalPeptides

# Creation of datasets for SignalP 6.0 analysis

1. ## To downloaded annotated proteins with signal peptide from UniProt (Primary accession number), run the script **Datasets/UniProt/uniprot_dataset.sh**
   ## To have metadata from UniProt, run the script **Datasets/Uniprot/UniProt_request_metadata.py**
   ## To downloaded proteins (Accession numbers) with signal peptide from the site signalpeptide.de, run the script **Datasets/signalpeptide/signalpeptide.py**
   ## To have metadata from Uniprot for the site signalpeptide dataset, run the script **Datasets/signalpeptide/uniprot_request_metadata_signalpeptide.py**

2. ## To create images and make a statistical analysis, run the script **Datasets/Analysis/DatasetsR.r**

# Creation of fasta files to run SignalP 6.0

3. ## To download a full protein sequence in fasta format for each PrimaryAccession number from files .csv with UniProt metadata, run the script **SignalP/seq_full_prot_fasta_sec_acc.py**

4. ## For removing the duplicates from downloaded fasta files, run the script **SignalP/fasta_for_signalp.sh**

5. ## Run "SignalP 6.0" using fasta files like an input (recommended to use an environement like conda)

signalp6 --fastafile /path/to/input.fasta --organism other --output_dir path/to/be/saved --format txt --mode slow-sequencial

# Results from SignalP analysis

6. ## To create images and make a statistical analysis from SignalP 6.0 results, run the script **SignalP/Analysis_results_SignalP/SignalPR.r** 

__________________________________________________________________________________________________________________________________________________



Run script taxon_search/**organism_taxonID.py**, an input file needs to contain the protein accession number

# To find a class name of organism (or host organism) via id of taxon

Rin script taxon_search/***taxon_class.sh**  with a file containing taxID like an input

# To download the sequences of confirmed signal peptides in tsv format, run the folowing command (taxonomy_id:40674 - for mammalia) :

# To download the taxonID of hosts

Rin script taxon_search/**host_taxonID.py**  with a file containing ProteinID like an input 

# To download host name via proteinId
 Run script  taxon_search/**host_vir_sigprot.py**
 
 #To download the verified protein sequences from UniProt
 
 Run **all_virSP_confirmed_uniprot.py**

 #To download metadata from UniProt
 
  Run **uniprot_request_metadata.py**
