# Dermentzaki_CellReports_2024
This repository contains the code required to reproduce the analyses presented in Figure 4 and Figure S2 of _Dermentzaki et al - Cell Reports - 2024_.
Noticeably, these scripts are also part of the supplemental files of the paper; Data S1.

Content of the repository:

- dataAnalysis.R: R script to reproduce all the analyses done on the nanopore dRNA-seq data except for those at single molecule resolution which require also the code reported in the singleMolecule scripts.
- singleMolecule.R: R script to identify a specific read matching dwell-time and current intensity from the Nanocompore database with the counterparts from nanopolish eventalign.
- singleMolecule.py: python script to extract the modification probability of TARDBP m6A sites.
- counts.rds: Matrix reporting transcripts counts.
- TARDBP_sitesData.rds: clustering data from Nanocompore for TARDBP m6A sites at single molecule resolution (result of singleMolecule.R).
- TARDBP_sitesProbabilities.rds: modification probability for TARDBP m6A sites at single molecule resolution (result of singleMolecule.R).
- memeResults.zip: archive with the results of the meme online suite.
- m6Anet.R: R script to reproduce the analyses on the results provided by m6Anet.
- m6Anet_KO_data.site_proba.csv: results from m6Anet on one KO sample.
- m6Anet_WT_data.site_proba.csv: results from m6Anet on one WT sample.