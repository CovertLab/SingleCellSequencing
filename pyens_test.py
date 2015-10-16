import pyensembl

mouse_genome = pyensembl.Genome(reference_name = 'GRCm38', gtf_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz', transcript_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz', protein_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz')
genes = mouse_genome.gene_names()

print genes