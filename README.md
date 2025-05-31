When Age Tips the Balance: a Dual Mechanism Affecting Hemispheric
Specialization for Language
================

[![DOI](https://zenodo.org/badge/726223472.svg)](https://zenodo.org/doi/10.5281/zenodo.10253278)

------------------------------------------------------------------------

## Contents

- [Background](#background)
- [Reference](#reference)
- [Code Release](#code-release)
- [Atlas Used](#atlas-used)
- [Other related papers that might interest
  you](#other-related-papers-that-might-interest-you)
- [Questions](#questions)

------------------------------------------------------------------------

## Background

**Aging** is accompanied by changes in brain architecture that *alter
the lateralization of functional networks*. In this study, we examined
how hemispheric specialization changes across the adult lifespan by
analyzing resting-state fMRI and structural MRI data from 728 typical
adults aged 18 to 88. Using the **Language-and-Memory Network** atlas,
we quantified regional asymmetries in functional connectivity along the
cortex’s principal gradient, and normalized regional volumes across 37
bilateral regions. We identified two distinct age-related asymmetry
trajectories: one pattern revealed a bilateralization of
language-dominant regions, while the other showed increasing leftward
specialization in multimodal regions associated with memory and
language. **These opposing patterns emerged around midlife** and were
*linked to performance in language production tasks*. By integrating
connectivity gradients, structural asymmetries, and behavioral data, our
findings provide new evidence for a dual mechanism reshaping functional
brain lateralization with age and demonstrate the utility of
resting-state metrics in tracking these shifts.

<p align="center">
<img src="readme_files/summary_LanguAging.png" width="70%" height="70%" />
</p>

------------------------------------------------------------------------

## Reference

For usage of the ***manuscript***, please cite:

- Roger, E.$^{†}$, **Labache, L.$^{†}$**, Hamlin, N., Kruse, J., Baciu,
  M., & Doucet, G. E. (2023). When age tips the balance: a dual
  mechanisms affecting hemispheric specialization for language.
  *BioRxiv* (2023). DOI:
  [10.1101/2023.12.04.569978](https://doi.org/10.1101/2023.12.04.569978).
  $^{†}$ these authors contributed equally.

For usage of the associated ***code***, please also cite:

- **Labache, L.**, Roger, E., Hamlin, N., Kruse, J., Baciu, M., &
  Doucet, G. E. (2023). When age tips the balance: a dual mechanisms
  affecting hemispheric specialization for language.
  loiclabache/RogerLabache_2023_LanguAging. DOI:
  [10.5281/zenodo.10253278](https://zenodo.org/doi/10.5281/zenodo.10253278)
- The implementations of the *Generalized Additive Mixed Models* and
  *Partition Around Medoids algorithm* were adapted from Roe, J., et
  al. 2021. DOI:
  [10.1038/s41467-021-21057-y](https://doi.org/10.1038/s41467-021-21057-y),
  and related GitHub repository:
  [AgeSym](https://github.com/jamesmroe/AgeSym).

------------------------------------------------------------------------

## Code Release

The `Script` folder includes three `R` scripts. The three `R` scripts
are designed to facilitate the replication of results as detailed in the
`Method Section` of the **manuscript**.

- `1_GAMM_hROIs.R`: `R` script to model gradient asymmetry trajectories
  throughout life using factor-smooth Generalized Additive Mixed Models.
  The script allows to compute the asymmetry trajectories underlying the
  interaction *Hemisphere×Age* and their confidence intervals. This
  script also assesses the significance of the smooth *Hemisphere×Age*
  interaction by testing for a difference in the smooth term of *Age*
  between hemispheres. We applied a False Discovery Rate correction to
  control for the number of tests conducted.
- `2_PAM_Clustering.R`: `R` script to classify regions in the
  **Language-and-Memory** network that demonstrate a significant
  *Hemisphere×Age* interaction, based on their functional asymmetry
  skewness profiles. This script also allows to compute the intersection
  point between the two average clusters curves.
- `3_CCA_BrainCognitionAssociation.R`: `R` script to proceed with the
  Canonical Correlation Analysis to assess brain–behavior Associations.

The `Data` folder contains output files generated from the *Generalized
Additive Mixed Models* analysis.

Note that the `R` scripts also contain the code **to reproduce the
figures found in the manuscript**. The brain renderings in the paper
require a customized version of [Surf
Ice](https://www.nitrc.org/projects/surfice/) that we will be happy to
share on demand.

------------------------------------------------------------------------

## Atlas Used

The atlas used in the paper is available in the `Atlas` folder.

- The **Language-and-Memory atlas** provides an atlas in standardized
  MNI volume space of 74 sentence- and memory-related areas (37 by
  hemisphere). The Language-and-Memory atlas encompasses the core
  regions that compose the stable components for language and memory.
  The Language-and-Memory atlas is composed of multiple brain regions
  provided by task-fMRI: one cross-sectional study for language (see
  [Labache, L., et al. 2019](https://doi.org/10.1007/s00429-018-1810-2),
  Github repository:
  [SENSAAS](https://github.com/loiclabache/SENSAAS_brainAtlas)) and one
  meta-analysis for memory (see [Spaniol, J., et
  al. 2009](https://doi.org/10.1016/j.neuropsychologia.2009.02.028)).
  The compilation of the Language-and-Memory atlas was initially
  undertaken in the following paper: [Roger, E., et
  al. 2020](https://doi.org/10.1002/hbm.24839).
  - The file `Atlas/language_memory_atlas.nii.gz` contains the
    `Volumetric` Language-and-Memory atlas (in MNI ICBM 152 space).
  - `Atlas/language_memory_atlas.txt`: text file containing a full
    description of each Language-and-Memory areas. The first column
    *Abbreviation* is the abbreviation of a region. The second column
    *Region* is the full anatomical label of a region. *Hemisphere*
    refers to the cerebral hemisphere to which a region belongs (“L” for
    left, “R” for right). *Function* indicates if a regions process
    language (“L”), memory (“M”), or language and memory (“LM”). *Index*
    is the index of each region that is used in the `NIfTI` file.
    *Number of Voxels* is the number of voxels of each region for the
    2mm version of the atlas. The MNI coordinate (columns *Xmm*, *Ymm*,
    *Zmm*) of each regions centroid is also provided.

<p align="center">
<img src="readme_files/LuM.gif" width="50%" height="50%" />
</p>

------------------------------------------------------------------------

## Other related papers that might interest you

- Language-and-Memory Network seminal paper: Roger, E., et al. 2020.
  DOI: [10.1002/hbm.24839](https://doi.org/10.1002/hbm.24839)
- Influence of Language Lateralisation on Gradient Asymmetry: Labache,
  L., et al. 2023. DOI:
  [10.1038/s41467-023-39131-y](https://doi.org/10.1038/s41467-023-39131-y),
  and related GitHub repository:
  [Labache_2022_AO](https://github.com/loiclabache/Labache_2022_AO)
- Sentence Supramodal Areas Atlas; Labache, L., et al. 2019. DOI:
  [10.1007/s00429-018-1810-2](https://doi.org/10.1007/s00429-018-1810-2),
  and related GitHub repository:
  [SENSAAS](https://github.com/loiclabache/SENSAAS_brainAtlas)
- For additional reading on GAMMs, consult Gavin Simpson’s procedure for
  comparing smooth terms: [Comparing smooths in factor-smooth
  interactions
  (1/2)](https://fromthebottomoftheheap.net/2017/10/11/difference-splines-i/),
  and [Comparing smooths in factor-smooth interactions
  (2/2)](https://fromthebottomoftheheap.net/2017/12/14/difference-splines-ii/)

------------------------------------------------------------------------

## Questions

Please contact me (Loïc Labache) at: <loic.labache@yale.edu> and/or
<loic.labache@ensc.fr>, or Élise Roger at: <elise.roger@umontreal.ca>
