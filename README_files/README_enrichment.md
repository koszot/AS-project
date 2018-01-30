# Armillaria ostoyae
### Input:
- __p3_i2_t47428_Arm_ostoy_v2.prot__
- __isoforms.fpkm_tracking__ : A CuffDiff analízis output fájlja ami az izoformák expressziós értékeit tartalmazza
- __genes.fpkm_tracking__ : A CuffDiff analízis output fájlja ami a gének expressziós értékeit tartalmazza
### Output:
- __p3_i2_t47428_Arm_ostoy_v2.prot.tsv__ : Az InterProScan domaineket tartalmazó output fájlja
- __aostoyae_GO_dict.tsv__ : A géneket és a GO-kat tartalmazó táblázat

Lefuttatjuk az InterProScan-t a GeneCatalog génjeire.
```
interproscan.sh -i p3_i2_t47428_Arm_ostoy_v2.prot -f tsv --iprlookup --goterms
```
Készítünk egy táblázatot ami tartalmazza a génekhez tartozó GO azonosítókat az enrichment számára.
```
ENRICHMENT_preparations_aostoyae.R
```
FONTOS: A következő két scriptet újra kell futtatni mert megváltozott a definíciója a DEVREG-nek és az AS-nek is.
```
ENRICHMENT_groups_aostoyae.R
ENRICHMENT_fisher_classic_aostoyae.R
```
# Coprinopsis cinerea (AmutBmut)
### Input:
- __Auramp1_GeneCatalog_proteins_20160719.aa.fasta__
- __Auramp1_GeneCatalog_genes_20160719.gff__
- __isoforms.fpkm_tracking__ : A CuffDiff analízis output fájlja ami az izoformák expressziós értékeit tartalmazza
- __genes.fpkm_tracking__ : A CuffDiff analízis output fájlja ami a gének expressziós értékeit tartalmazza
### Output:
- __Auramp1_GeneCatalog_proteins_20160719.aa.fasta.tsv__ : Az InterProScan domaineket tartalmazó output fájlja
- __aampla_GO_dict.tsv__ : A géneket és a GO-kat tartalmazó táblázat
- __aampla_dictionary.tsv__ : A protein ID/Transcript ID megfeleltetéshez szükséges táblázat

Kitöröljük a *-okat, mivel az InterProScan nem tudja értelmezni
```
perl -pi -e 's/\*//g' Auramp1_GeneCatalog_proteins_20160719.aa.fasta
```
Mivel nem fognak stimmelnek az ID-k a fasta fájl proteinjei és az alap GTF fájl ID-jai között ezért meg kell egymással feleltetni őket. Ehhez először kivágjuk a __Auramp1_GeneCatalog_proteins_20160719.aa.fasta__ átalakítatlan headerjeit -> __Auramp1_GeneCatalog_proteins_20160719.aa.fasta.headers__

Majd átalakítjuk a headereket tartalmazó fájlt táblázattá
```
>jgi\|Auramp1\| --> ""
\| --> \t
```
Lefuttajuk a szótárt elkészítő R scriptet
```
ENRICHMENT_dictionary_aampla.R
```
Átalakítjuk a proteineket tartalmazó fasta fájlt az InterProScan számára a könyebb áttekinthetőségért
```
# jgi\|Auramp1\| --> ""
# \|.*$ --> ""
```
Lefuttatjuk az InterProScan-t a GeneCatalog génjeire.
```
interproscan.sh -i Auramp1_GeneCatalog_proteins_20160719.aa.fasta -f tsv --iprlookup --goterms -dp
```
Készítünk egy táblázatot ami tartalmazza a génekhez tartozó GO azonosítókat az enrichment számára.
```
ENRICHMENT_preparations_aampla.R
```
FONTOS: A következő két scriptet újra kell futtatni mert megváltozott a definíciója a DEVREG-nek és az AS-nek is.
```
ENRICHMENT_groups_aampla.R
ENRICHMENT_fisher_classic_aampla.R
```

### Coprinopsis cinerea

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta

# lekell futtatni az IPR scant a genecatalog génekre

interproscan.sh -i Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta -f tsv --iprlookup --goterms

# mivel nem stimmelnek az ID-k a fasta fájl proteinjei és az alap GTF fájl ID-jai között ezért meg kell egymással feleltetni őket
# ehhez először kivágjuk a Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta átalakítatlan headerjeit -> Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta.headers

# >jgi\|Copci_AmutBmut1\| --> ""
# \| --> \t

ccinerea_AmutBmut_enrichment_dictionary.R

ccinerea_AmutBmut_enrichment_preparations.R

ccinerea_AmutBmut_genelvl_enrichment_groups.R
```

### Cryptococcus neoformans

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' Cryptococcus_neoformans_H99.proteins.fasta

# jgi\|Cryne_H99_1\| --> ""
# \|.*$ --> ""

# lekell futtatni az IPR scant a genecatalog génekre

interproscan.sh -i Cryptococcus_neoformans_H99.proteins.fasta -f tsv --iprlookup --goterms -dp

# mivel nem stimmelnek az ID-k a fasta fájl proteinjei és az alap GTF fájl ID-jai között ezért meg kell egymással feleltetni őket
# ehhez először kivágjuk a Cryptococcus_neoformans_H99.proteins.fasta átalakítatlan headerjeit -> Cryptococcus_neoformans_H99.proteins.fasta.headers

# >jgi\|Cryne_H99_1\| --> ""
# \| --> \t

cneoformans_enrichment_dictionary.R

cneoformans_enrichment_preparations.R

cneoformans_genelvl_enrichment_groups.R
```

### Lentinus tigrinus

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' Sisbr1_GeneCatalog_proteins_20130805.aa.fasta

# lekell futtatni az IPR scant a genecatalog génekre

interproscan.sh -i Sisbr1_GeneCatalog_proteins_20130805.aa.fasta -f tsv --iprlookup --goterms

# mivel nem stimmelnek az ID-k a fasta fájl proteinjei és az alap GTF fájl ID-jai között ezért meg kell egymással feleltetni őket
# ehhez először kivágjuk a Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta átalakítatlan headerjeit -> Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta.headers

# >jgi\|Sisbr1\|--> ""
# \| --> \t

ltigrinus_enrichment_dictionary.R

ltigrinus_enrichment_preparations.R

ltigrinus_genelvl_enrichment_groups.R
```

### Phanerochaete chrysosporium

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' Phchr2_GeneCatalog_proteins_20131210.aa.fasta

# >jgi\|Phchr2\| --> ""
# \|.*$ --> ""

# lekell futtatni az IPR scant a genecatalog génekre

interproscan.sh -i Phchr2_GeneCatalog_proteins_20131210.aa.fasta -f tsv --iprlookup --goterms -dp

# mivel nem stimmelnek az ID-k a fasta fájl proteinjei és az alap GTF fájl ID-jai között ezért meg kell egymással feleltetni őket
# ehhez először kivágjuk a Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta átalakítatlan headerjeit -> Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta.headers

# >jgi\|Phchr2\| --> ""
# \| --> \t

pchrysosporium_enrichment_dictionary.R

pchrysosporium_enrichment_preparations.R

pchrysosporium_genelvl_enrichment_groups.R
```

### Rickenella mellea

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' Ricmel1_GeneCatalog_proteins_20151108.aa.fasta

# jgi\|Ricmel1\| --> ""
# \|.*$ --> ""

# lekell futtatni az IPR scant a genecatalog génekre

interproscan.sh -i Ricmel1_GeneCatalog_proteins_20151108.aa.fasta -f tsv --iprlookup --goterms -dp

# mivel nem stimmelnek az ID-k a fasta fájl proteinjei és az alap GTF fájl ID-jai között ezért meg kell egymással feleltetni őket
# ehhez először kivágjuk a Ricmel1_GeneCatalog_proteins_20151108.aa.fasta átalakítatlan headerjeit -> Ricmel1_GeneCatalog_proteins_20151108.aa.fasta.headers

# >jgi\|Ricmel1\| --> ""
# \| --> \t

rmellea_enrichment_dictionary.R

rmellea_enrichment_preparations.R

rmellea_genelvl_enrichment_groups.R
```

### Schizophyllum commune

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' Schco3_GeneCatalog_proteins_20130812.aa.fasta

# lekell futtatni az IPR scant a genecatalog génekre

interproscan.sh -i Schco3_GeneCatalog_proteins_20130812.aa.fasta -f tsv --iprlookup --goterms

# mivel nem stimmelnek az ID-k a fasta fájl proteinjei és az alap GTF fájl ID-jai között ezért meg kell egymással feleltetni őket
# ehhez először kivágjuk a Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta átalakítatlan headerjeit -> Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta.headers

# >jgi\|Schco3\| --> ""
# \| --> \t

scommune_enrichment_dictionary.R

scommune_enrichment_preparations.R

scommune_genelvl_enrichment_groups.R
```

### Ustilago maydis

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' p3_t237631_Ust_maydi_v2GB.prot

interproscan.sh -i p3_t237631_Ust_maydi_v2GB.prot -f tsv --iprlookup --goterms -dp

umaydis_enrichment_preparations.R

umaydis_genelvl_enrichment_groups.R
```

# MCL Clustering a géncsaládok megkeresésére

```
# megkeressük a géncsaládokat amikből tovább tudunk haladni a primary transcript megtalálása felé

aampla 			-> 	Auramp1_GeneCatalog_proteins_20160719.aa.fasta
aostoyae 		->	p3_i2_t47428_Arm_ostoy_v2.prot
ccinerea 		->	Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta
ltigrinus 		->	Sisbr1_GeneCatalog_proteins_20130805.aa.fasta
pchrysosporium 	->	Phchr2_GeneCatalog_proteins_20131210.aa.fasta
rmellea			->	Ricmel1_GeneCatalog_proteins_20151108.aa.fasta
scommune 		->	Schco3_GeneCatalog_proteins_20130812.aa.fasta

# orthomclAdjustFasta

cd compliantFasta
orthomclAdjustFasta aampla ../aampla.fasta 1
orthomclAdjustFasta aostoyae ../aostoyae.fasta 1
orthomclAdjustFasta ccinerea ../ccinerea.fasta 1
orthomclAdjustFasta ltigrinus ../ltigrinus.fasta 1
orthomclAdjustFasta pchrysosporium ../pchrysosporium.fasta 1
orthomclAdjustFasta rmellea ../rmellea.fasta 1
orthomclAdjustFasta scommune ../scommune.fasta 1

# orthomclFilterFasta

orthomclFilterFasta ./compliantFasta 10 20

# all VS all BLAST

mpiformatdb -N 32 -i goodProteins.fasta -o T
mv goodProteins.fasta.* ~/share/
nohup mpirun -n 32 mpiblast -d goodProteins.fasta -i goodProteins.fasta -p blastp -m 8 -o goodProteins_blasted --copy-via=none &

# orthomclBlastParser

orthomclBlastParser goodProteins_blasted ./compliantFasta >> similarSequences.txt


# mivel nem stimmelnek a transcript_id-k a protein_id-s baszakodás miatt ezért az egyes fajok dictionaryje alapján át kell alakítani a családokat tartalmazó fájlt

cluster_transformation.R
```

# Silix/Hifix és enrichment (OUTDATED)

```
# a clusterezés megkezdéséhez egy fileba kell olvasztani az összes különálló faj fasta fájlját
# hozzá kell adni az adott faj headerjéhez a faj nevét hogy meg lehessen különböztetni

aampla 			-> 	Auramp1_GeneCatalog_proteins_20160719.aa.fasta
aostoyae 		->	p3_i2_t47428_Arm_ostoy_v2.prot
ccinerea 		->	Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta
cneoformans		->	Cryptococcus_neoformans_H99.proteins.fasta
ltigrinus 		->	Sisbr1_GeneCatalog_proteins_20130805.aa.fasta
pchrysosporium 	->	Phchr2_GeneCatalog_proteins_20131210.aa.fasta
rmellea			->	Ricmel1_GeneCatalog_proteins_20151108.aa.fasta
scommune 		->	Schco3_GeneCatalog_proteins_20130812.aa.fasta
umaydis			->	p3_t237631_Ust_maydi_v2GB.prot

# mpiblast database

mpiformatdb -N 48 -i all_species_protein.fasta -o T

# áthelyezzük a fájlokat ahogy kell

mpirun -n 48 mpiblast -d all_species_protein.fasta -i all_species_protein.fasta -p blastp -m 8 -o all_species_protein_blasted --copy-via=none

# Silix a Hifix-hez

silix all_species_protein.fasta all_species_protein_blasted --net > all_species_protein_SLX.fnodes

# Hifix KITERJESZTÉSNEK .fasta-nak kell lennie!

hifix -t 48 all_species_protein.fasta all_species_protein_blasted.net all_species_protein_SLX.fnodes > all_species_protein_HFX.fnodes

# p3_i2_t47428_Arm_ostoy_v2_ --> semmire csere

# mivel nem stimmelnek a transcript_id-k a protein_id-s baszakodás miatt ezért az egyes fajok dictionaryje alapján át kell alakítani a családokat tartalmazó fájlt

cluster_transformation.R

# fajonkénti enrichment

enrichment_scripts/aampla_genelvl_enrichment_fisher_classic.R
enrichment_scripts/aostoyae_genelvl_enrichment_fisher_classic.R
enrichment_scripts/ccinerea_AmutBmut_genelvl_enrichment_fisher_classic.R
enrichment_scripts/cneoformans_genelvl_enrichment_fisher_classic.R
enrichment_scripts/ltigrinus_genelvl_enrichment_fisher_classic.R
enrichment_scripts/pchrysosporium_genelvl_enrichment_fisher_classic.R
enrichment_scripts/rmellea_genelvl_enrichment_fisher_classic.R
enrichment_scripts/scommune_genelvl_enrichment_fisher_classic.R
enrichment_scripts/umaydis_genelvl_enrichment_fisher_classic.R

# összesítés

enrichment_scripts/families_count.R
enrichment_scripts/GO_count.R
```