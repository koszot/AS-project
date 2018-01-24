

####################################################
#     MCL CLUSTERING - GÉNCSALÁDOK MEGKERESÉSE     #
####################################################

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








========= Step 9: orthomclLoadBlast  ===========
Input:
  - similarSequences.txt
Output:
  - SimilarSequences table in the database

This step loads the BLAST results into the orthomcl database. 

Use the orthomclLoadBlast program for this. 

NOTE: You might get the following error when you run this command:

  "The used command is not allowed with this MySQL version."

The SQL that causes this is LOAD DATA LOCAL INFILE.  MySql needs specific configuration to enable this command.  See these two pages:

http://dev.mysql.com/doc/refman/5.1/en/load-data-local.html
http://dev.mysql.com/doc/refman/5.0/en/loading-tables.html

In sum:
  1) you need to start the MySql server with the option:  --local_infile=1
  2) during installation, MySql needs to be have been compiled with:  --enable-local-infile.  

It is possible that your MySql was compiled with that flag.  It is included by default in some distributions of MySql.   The only way we know of to find out what compile flags were used is to try the mysqlbug command.  This command opens an email so you can report a bug.  Apparently at the bottom of the mail is a list of compile flags.  Once you see them you can abort the mail.  If your mysql was not compiled with that flag it will need to be reinstalled.


Benchmark time: 4 hours


========= Step 10: orthomclPairs =========
Input:
  - SimilarSequences table in the database
Output:
  - PotentialOrthologs table
  - PotentialInParalogs table
  - PotentialCoOrthologs table

This is a computationally major step that finds protein pairs.  It executes the algorithm described in the OrthoMCL Algorithm Document (docs.google.com/Doc?id=dd996jxg_1gsqsp6), using a relational database.  The program proceeds through a series of internal steps, each creating an intermediate database table or index.  There are about 20 such tables created. Finally, it populates the output tables.

The cleanup= option allows you to control the cleaning up of the intermediary tables.  The 'yes' option drops the intermediary tables once they are no longer needed.  The 'no' option keeps the intermediary tables in the database.  In total, they are expected to be about 50 percent of the SimilarSequences table. They are useful mostly for power users or developers who would like to query them. They can be removed afterwards with the 'only' or 'all' options.  The latter also removes the final tables, and should only be done after Step 11 below has dumped them to files.

The startAfter= option allows you to pick up where the program left off, if it stops for any reason.  Look in the log to find the last completed step, and use its tag as the value for startAfter=

Because this program will run for many hours, we recommend you run it using the UNIX 'screen' program, so that it does not abort in the middle.  (If it does, use startAfter=).

Benchmark time: 16 hours


========== Step 11: orthomclDumpPairsFiles ========
Input:
  - database with populated pairs tables
Output
  - pairs/ directory.  
  - mclInput file

Run the orthomclDumpPairsFiles.

The pairs/ directory contains three files: ortholog.txt, coortholog.txt, inparalog.txt.  Each of these has three columns:
   - protein A
   - protein B
   - their normalized score (See the Orthomcl Algorithm Document).

Benchmark time: 5 min


========== Step 12: mcl ========
Input:
  - mclInput file
Output:
  - mclOutput file

mcl my_orthomcl_dir/mclInput --abc -I 1.5 -o my_orthomcl_dir/mclOutput

Benchmark time: 3 hours


========== Step 13: orthomclMclToGroups ==========
Input:
  - mclOutput file
Output:
  - groups.txt

Change to my_orthomcl_dir and run:
  orthomclMclToGroups my_prefix 1000 < mclOutput > groups.txt

my_prefix is an arbitrary string to use as a prefix for your group IDs.

1000 is an arbitrary starting point for your group IDs.

Benchmark time: 1 min

























# mivel nem stimmelnek a transcript_id-k a protein_id-s baszakodás miatt ezért az egyes fajok dictionaryje alapján át kell alakítani a családokat tartalmazó fájlt

cluster_transformation.R







###############################################
#     SILIX/HIFIX CLUSTERING ÉS ENRICHMENT    #
###############################################

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

##################
#     SUMMARY    #
##################

# egyenként megcsináljuk a különböző fajok summaryjeit

summary_scripts/aampla_summary.R
summary_scripts/aostoyae_summary.R
summary_scripts/ccinerea_AmutBmut_summary.R
summary_scripts/cneoformans_summary.R
summary_scripts/ltigrinus_summary.R
summary_scripts/pchrysosporium_summary.R
summary_scripts/rmellea_summary.R
summary_scripts/scommune_summary.R
summary_scripts/umaydis_summary.R

# összesítjük a summaryket

summary_scripts/summary_all.R





