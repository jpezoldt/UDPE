module use /software/module/;
module add UHTS/Analysis/samtools/1.3;
module add Development/java_jdk/1.8.0_102;

# Human

wget -O Homo_sapiens.GRCh38.87.gtf.gz ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz
gzip -d Homo_sapiens.GRCh38.87.gtf.gz
wget -O Homo_sapiens.GRCh38.87.dna.primary_assembly.fa.gz ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh38.87.dna.primary_assembly.fa.gz
wget -O Homo_sapiens.GRCh38.87.vcf.gz ftp://ftp.ensembl.org/pub/release-87/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz

samtools faidx Homo_sapiens.GRCh38.87.dna.primary_assembly.fa
java -jar picard.jar CreateSequenceDictionary R=Homo_sapiens.GRCh38.87.dna.primary_assembly.fa O=Homo_sapiens.GRCh38.87.dna.primary_assembly.dict

#Fly

wget -O Drosophila_melanogaster.BDGP6.87.gtf.gz ftp://ftp.ensembl.org/pub/release-87/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.87.gtf.gz
gzip -d Drosophila_melanogaster.BDGP6.87.gtf.gz
wget -O Drosophila_melanogaster.BDGP6.87.dna.top_level.fa.gz ftp://ftp.ensembl.org/pub/release-87/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz
gzip -d Drosophila_melanogaster.BDGP6.87.dna.top_level.fa.gz

samtools faidx Drosophila_melanogaster.BDGP6.87.dna.top_level.fa
java -jar picard.jar CreateSequenceDictionary R=Drosophila_melanogaster.BDGP6.87.dna.top_level.fa O=Drosophila_melanogaster.BDGP6.87.dna.top_level.dict

# Mouse

wget -O Mus_musculus.GRCm38.87.gtf.gz ftp://ftp.ensembl.org/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf.gz
gzip -d Mus_musculus.GRCm38.87.gtf.gz
wget -O Mus_musculus.GRCm38.87.dna.primary_assembly.fa.gz ftp://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gzip -d Mus_musculus.GRCm38.87.dna.primary_assembly.fa.gz

samtools faidx Mus_musculus.GRCm38.87.dna.primary_assembly.fa
java -jar picard.jar CreateSequenceDictionary R=Mus_musculus.GRCm38.87.dna.primary_assembly.fa O=Mus_musculus.GRCm38.87.dna.primary_assembly.dict
