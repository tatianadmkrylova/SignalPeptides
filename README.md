## SignalPeptides

1. Make a list of interested proteins ("primaryAccession" of protein or "id" in UniProt)

2. Downloaded a full protein sequence in fasta format for each id -> script **seq_full_prot_fasta.py**

# if you need only confirmed: make a fasta file which contains only confirmed sequences  - full proteins (check if proteinID in a column "confirmed")

Rin script fasta_sequences/**seq_sp_confirmed.py**

3. Run "SignalP 6.0" using created file like an input (recommended to use an environement like conda)

signalp6 --fastafile /path/to/input.fasta --organism other --output_dir path/to/be/saved --format txt --mode fast

## Results analysis

# To convert a file gff3 to dataset in R, run this command in RStudio

**read.table("output.gff3", sep = "\t", header = "FALSE", comment.char = "#", stringsAsFactors = FALSE)**

# To download a taxId for each protein via ProteinID

Run script taxon_search/**organism_taxonID.py**, an input file needs to contain the protein accession number

# To find a class name of organism (or host organism) via id of taxon

Rin script taxon_search/***taxon_class.sh**  with a file containing taxID like an input

# To download the sequences of confirmed signal peptides in tsv format, run the folowing command (taxonomy_id:40674 - for mammalia) :

** curl -L -o uniprot_mammalia_signal_confirmed_exp.tsv \
"https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:40674%20AND%20ft_signal_exp:*&format=tsv&fields=accession,id,protein_name,gene_names,organism_name,length,ft_signal"**

# To download the taxonID of hosts

Rin script taxon_search/**host_taxonID.py**  with a file containing ProteinID like an input 

# To download host name via proteinId
 Run script  taxon_search/**host_vir_sigprot.py**
 
 #To download the verified protein sequences from UniProt
 
 Run **all_virSP_confirmed_uniprot.py**
