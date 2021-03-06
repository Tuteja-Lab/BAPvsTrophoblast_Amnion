# BAPvsTrophoblast_Amnion

## [Tuteja Lab](https://www.tutejalab.org), Iowa State University, Ames IA-50011

**The HTML version of this documentation is available [here](https://tuteja-lab.github.io/BAPvsTrophoblast_Amnion)**:

### Documentation

1. [Section 1: Preparing datasets for anlayses](1_input-data-prep_v2.Rmd)
2. [Section 2: DESeq2 analyses](2_de-analyses_v2.Rmd)
3. [Section 3: PCA plots for the datasets](3_pca-plots_v2.Rmd)
4. [Section 4: Finding novel amnion markers](4_heatmaps_v2.Rmd)
5. [Section 5: sn/scRNAseq data analyses](5_snRNAseq_v2.Rmd)
6. [Section 6: Custom encrichment dataset preparation](6_custom-PCE_v2.Rmd)

### Data availability

The data used in this study are publicly available:
| BioProject   | Assay Type | source_name                                                | replicates | treatment                                       | Publication     | Abbreviation          |
|--------------|------------|------------------------------------------------------------|------------|-------------------------------------------------|-----------------|-----------------------|
| PRJNA605646  | RNA-Seq    | Naïve H9 ESC                                               | 2          | Naive human PSC   condition                     | Io et al.       | nH9_Io                |
| PRJNA605646  | RNA-Seq    | Naïve H9 ESC-derived   nTE_D1                              | 2          | ABPJ treatment for 1   day                      | Io et al.       | H9_nTE_D1_Io          |
| PRJNA605646  | RNA-Seq    | Naïve H9 ESC-derived   nTE_D2                              | 2          | ABPJ treatment for 2   days                     | Io et al.       | H9_nTE_D2_Io          |
| PRJNA605646  | RNA-Seq    | term Amnion   (Cont2,3,4,4b,6,7,8,9,11,12,14,15 RNA-seq)   | 2          | ABPJ treatment for 3   days                     | Io et al.       | H9_nTE_D3_Io          |
| PRJNA605646  | RNA-Seq    | Naïve H9 ESC-derived   nCT_P3                              | 2          | ACE treatment                                   | Io et al.       | H9_nCT_P3_Io          |
| PRJNA605646  | RNA-Seq    | Naïve H9 ESC-derived   nCT_P10                             | 2          | ACE treatment                                   | Io et al.       | H9_nCT_P10_Io         |
| PRJNA605646  | RNA-Seq    | Naïve H9 ESC-derived   nCT_P15                             | 2          | ACE treatment                                   | Io et al.       | H9_nCT_P15_Io         |
| PRJNA605646  | RNA-Seq    | Naïve PSC (Guo et al.   Development 2017)-derived nCT_P15  | 2          | ACE treatment                                   | Io et al.       | cR_nCT_P15_Io         |
| PRJNA605646  | RNA-Seq    | 409B2 hiPSC-derived   nCT_P15                              | 2          | ACE treatment                                   | Io et al.       | iPS_nCT_P15_Io        |
| PRJNA605646  | RNA-Seq    | 1st Trimester_TSC   (Okae et al) (CT30)                    | 2          | Okae et al. Cell Stem   Cell 2018               | Io et al.       | CT_TSC_Okae           |
| PRJNA605646  | RNA-Seq    | nST from Naïve H9   ESC-derived nCT                        | 2          | modified from Okae et   al. Cell Stem Cell 2018 | Io et al.       | H9_nSCT_Io            |
| PRJNA605646  | RNA-Seq    | nEVT from Naïve H9   ESC-derived nCT                       | 2          | modified from Okae et   al. Cell Stem Cell 2018 | Io et al.       | H9_nEVT_Io            |
| PRJNA605646  | RNA-Seq    | Primed H9 ESC                                              | 2          | primed human PSC   condition                    | Io et al.       | pH9_Io                |
| PRJNA605646  | RNA-Seq    | pBAP_D1 from Primed   H9 ESC                               | 2          | MEF-CM with BAP   treatment                     | Io et al.       | H9_pBAP_D1_Io         |
| PRJNA605646  | RNA-Seq    | pBAP_D2 from Primed   H9 ESC                               | 2          | MEF-CM with BAP   treatment                     | Io et al.       | H9_pBAP_D2_Io         |
| PRJNA605646  | RNA-Seq    | pBAP_D3 from Primed   H9 ESC                               | 2          | MEF-CM with BAP   treatment                     | Io et al.       | H9_pBAP_D3_Io         |
| PRJNA605646  | RNA-Seq    | Cytotrophoblast at 7   gestational weeks                   | 2          | 7 gestational weeks                             | Io et al.       | CT_7wk                |
| PRJNA605646  | RNA-Seq    | Cytotrophoblast at 9   gestational weeks                   | 1          | 9 gestational weeks                             | Io et al.       | CT_9wk                |
| PRJNA605646  | RNA-Seq    | Cytotrophoblast at 11   gestational weeks                  | 1          | 11 gestational weeks                            | Io et al.       | CT_11wk               |
| PRJNA294733  | RNA-Seq    | H1 ESC-derived TB_BAP   D8 (large cells, >70 μm)           | 3          | BAP treatment                                   | Yabe et al.     | H1_BAP_D8_>70_Yabe    |
| PRJNA294733  | RNA-Seq    | H1 ESC-derived TB_BAP D8 (between   40-70 μm)              | 3          | BAP treatment                                   | Yabe et al.     | H1_BAP_D8_40-70_Yabe  |
| PRJNA294733  | RNA-Seq    | H1 ESC-derived TB_BAP   D8 (small cells, <40 μm)           | 3          | BAP treatment                                   | Yabe et al.     | H1_BAP_D8_<40_Yabe    |
| PRJNA294733  | RNA-Seq    | undifferentiated H1   ESC                                  | 3          | MEF-CM with FGF2   treatment                    | Yabe et al.     | H1_Yabe               |
| E-MTAB-10429 | RNA-Seq    | Blastocyst_TSC (Okae   et al) (BTS5, BTS11)                | 2          | Okae et al. Cell Stem   Cell 2018               | Sheridan et al. | BTS_TSC_Okae          |
| E-MTAB-10429 | RNA-Seq    | 1st Trimester_TSC   (Okae et al) (CT27, CT29, CT30)        | 3          | Okae et al. Cell Stem   Cell 2018               | Sheridan et al. | CT_TSC_Okae           |
| E-MTAB-10429 | RNA-Seq    | Blastocyst_EVT (Okae   et al) (BTS5, BTS11)                | 2          | Okae et al. Cell Stem   Cell 2018               | Sheridan et al. | BTS_EVT_Okae          |
| E-MTAB-10429 | RNA-Seq    | 1st Trimester_EVT (Okae et al) (CT27, CT29, CT30)          | 3          | Okae et al. Cell Stem   Cell 2018               | Sheridan et al. | CT_EVT_Okae           |
| E-MTAB-10429 | RNA-Seq    | Blastocyst_SCT (Okae   et al) (BTS5, BTS11)                | 2          | Okae et al. Cell Stem   Cell 2018               | Sheridan et al. | BTS_SCT_Okae          |
| E-MTAB-10429 | RNA-Seq    | 1st Trimester_TSC (Okae et al) (CT27, CT29, CT30)          | 3          | Okae et al. Cell Stem   Cell 2018               | Sheridan et al. | CT_SCT_Okae           |
| E-MTAB-10429 | RNA-Seq    | Organoid_villous   (Turco et al) (R053,   X011, X021,X035) | 4          | Turco et al. Nature   2018                      | Sheridan et al. | org_villous           |
| E-MTAB-10429 | RNA-Seq    | Organoid_EVT   (Turco et al) (R053,   X011, X021, x035)    | 4          | Turco et al. Nature   2018                      | Sheridan et al. | org_EVT               |
| PRJNA414247  | RNA-Seq    | undifferentiated   hESC line H9 cells                      | 2          | untreated                                       | Krendl et al.   | H9_Krendl             |
| PRJNA414247  | RNA-Seq    | hESC line H9 cells   after 8 hours of BMP4 treatment       | 2          | BMP4 exposed H9 cells                           | Krendl et al.   | H9_BMP4_8h_Krendl     |
| PRJNA414247  | RNA-Seq    | hESC line H9 cells   after 16 hours of BMP4 treatment      | 2          | BMP4 exposed H9 cells                           | Krendl et al.   | H9_BMP4_16h_Krendl    |
| PRJNA414247  | RNA-Seq    | hESC line H9 cells   after 24 hours of BMP4 treatment      | 2          | BMP4 exposed H9 cells                           | Krendl et al.   | H9_BMP4_D1_Krendl     |
| PRJNA414247  | RNA-Seq    | hESC line H9 cells   after 48 hours of BMP4 treatment      | 2          | BMP4 exposed H9 cells                           | Krendl et al.   | H9_BMP4_D2_Krendl     |
| PRJNA414247  | RNA-Seq    | hESC line H9 cells   after 72 hours of BMP4 treatment      | 2          | BMP4 exposed H9 cells                           | Krendl et al.   | H9_BMP4_D3_Krendl     |
| PRJNA352339  | RNA-Seq    | undifferentiated H9   ESC, feeder free                     | 3          | mTeSR                                           | Shao et al.     | H9_ESC_Shao           |
| PRJNA352339  | RNA-Seq    | H9 ESC-derived   amnion-like cells (3D amniogenesis)       | 3          | Shao et al. Nature   Mater 2017                 | Shao et al.     | H9_amnion_Shao        |
| PRJNA276463  | RNA-Seq    | CC2 Amnion_18_female                                       | 1          | 18 gestational weeks                            | Roost et al.    | amnion_18wk           |
| PRJNA276463  | RNA-Seq    | CJ1 Amnion_9_male                                          | 2          | 9 gestational weeks                             | Roost et al.    | amnion_9wk            |
| PRJNA276463  | RNA-Seq    | CI1 Amnion_16_male                                         | 2          | 16 gestational weeks                            | Roost et al.    | amnion_16wk           |
| PRJNA276463  | RNA-Seq    | CR1 Amnion_22_male                                         | 2          | 22 gestational weeks                            | Roost et al.    | amnion_22wk           |
| PRJNA316992  | RNA-Seq    | Term Amnion   (Cont2,3,4,4b,6,7,8,9,11,12,14,15 RNA-seq)   | 12         | Term fetal amnion                               | Suzuki et al.   | amnion_term           |
| PRJNA720858  | snRNA-Seq  | H1 ESC-derived TB_BAP   D8 (20% O2)_snRNA-seq              | 2          | BAP treatment                                   | Khan et al.     | H1_BAP_D8_20pcO2_Khan |
| PRJNA720858  | snRNA-Seq  | H1 ESC-derived TB_BAP   D8 (5% O2)_snRNA-seq               | 2          | BAP treatment                                   | Khan et al.     | H1_BAP_D8_5pcO2_Khan  |
| PRJNA705594  | scRNA-Seq  | Naïve H9 ESC-derived   nTE_D2                              | 1          | ABPJ treatment                                  | Io et al.       | H9_sc.nTE_D2_Io       |
| PRJNA705594  | scRNA-Seq  | Naïve H9 ESC-derived   nTE_D3                              | 1          | ABPJ treatment                                  | Io et al.       | H9_sc.nTE_D3_Io       |
| PRJNA705594  | scRNA-Seq  | Naïve H9 ESC-derived   nCT_D5                              | 1          | ACE treatment                                   | Io et al.       | H9_sc.nCT_D5_Io       |
| PRJNA705594  | scRNA-Seq  | Naïve H9 ESC-derived   nCT_D10                             | 1          | ACE treatment                                   | Io et al.       | H9_sc.nCT_D10_Io      |

### Other datasets

Syncytiotrophoblast (STB) genes expression were tested against [Zhou et. al.](https://pubmed.ncbi.nlm.nih.gov/31435013), dataset ([GSE109555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109555)).

Custom [PlacentaCellEnrich](https://placentacellenrich.gdcb.iastate.edu) were generated using: [Xiang _et. al._](https://pubmed.ncbi.nlm.nih.gov/31830756/)(NCBI BioProject: [PRJNA562548](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA562548") and [Castel _et al._](https://pubmed.ncbi.nlm.nih.gov/33238118) publication, using data from [Petropoulos _et al.,_](https://pubmed.ncbi.nlm.nih.gov/27062923) and [Zhou _et al.,_](https://pubmed.ncbi.nlm.nih.gov/31435013) downloaded from [BitBucket](https://gitlab.univ-nantes.fr/E137833T/Castel_et_al_2020)

### Contacts

1. Project related questions: Geetu Tuteja
2. Scripts and workflow related questions: please open an [issue](https://github.com/Tuteja-Lab/BAPvsTrophoblast_Amnion/issues/new) here on GitHub.

### Publication

If you use the scripts from this repository in your research, please cite this publication:

Seetharam, Arun S., Ha TH Vu, Sehee Choi, Teka Khan, Megan A. Sheridan, Toshihiko Ezashi, R. Michael Roberts, and Geetu Tuteja. "The product of BMP-directed differentiation protocols for human primed pluripotent stem cells is placental trophoblast and not amnion." _Stem Cell Reports_ (2022). Available online 19 May 2022.
[<img src="https://img.shields.io/badge/-Open_Access-blue?style=flat&logo=Open-Access"/>](https://doi.org/10.1016/j.stemcr.2022.04.014)
