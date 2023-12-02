When Age Tips the Balance: a Dual Mechanism Affecting Hemispheric
Specialization for Language
================

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

**Aging** engenders neuroadaptations, generally *reducing specificity
and selectivity* in functional brain responses. Our investigation delves
into the **functional specialization** of brain hemispheres within
language-related networks across adulthood. In a cohort of 728 healthy
adults spanning ages 18 to 88, we modeled the trajectories of
inter-hemispheric asymmetry concerning the principal functional gradient
across 37 homotopic regions of interest of an extensive language network
known as the **Language-and-Memory Network**. Our findings reveal that
over two-thirds of Language-and-Memory Network homotopic regions of
interest undergo asymmetry changes with age, falling into two main
clusters. The first cluster evolves from left-sided specialization to
right-sided tendencies, while the second cluster transitions from
right-sided asymmetry to left-hemisphere dominance. These reversed
asymmetry shifts manifest around midlife, occurring after age 50, and
are associated with poorer language production performance. Our results
provide valuable insights into the **influence of functional brain
asymmetries on language proficiency** and present a **dynamic
perspective on brain plasticity during the typical aging process**.

<p align="center">
<img src="readme_files/summary_LanguAging.png" width="70%" height="70%" />
</p>

------------------------------------------------------------------------

## Reference

For usage of the ***manuscript***, please cite:

- Roger, E., **Labache, L.**, Baciu, M., & Doucet, G. E. (2023). When
  age tips the balance: a dual mechanisms affecting hemispheric
  specialization for language. X (2023). DOI:
  [X/X-X-X-X](https://doi.org/X.X/X-X-X-X)

For usage of the associated ***code***, please also cite:

- Roger, E., **Labache, L.**, Baciu, M., & Doucet, G. E. (2023). When
  age tips the balance: a dual mechanisms affecting hemispheric
  specialization for language. loiclabache/RogerLabache_2023_LanguAging.
  DOI: [10.5281/zenodo.C](https://doi.org/10.5281/zenodo.C)

------------------------------------------------------------------------

## Code Release

The `Script` folder includes 4 `R` scripts alongside X **simulated**
`data files`. These simulated files serve as templates, demonstrating
the required data structure for successfully executing the code.
**Please note that these files contain simulated, *not real*, data.**

The 4 `R` scripts are designed to facilitate the replication of results
as detailed in the `Method section` of the **manuscript**.

- `1_GAMM_hROIs.R`: `R` script to model gradient asymmetry trajectories
  throughout life using factor-smooth Generalized Additive Mixed Models.
  The script allows to compute the asymmetry trajectories underlying the
  interaction *Hemisphere×Age* and their confidence intervals.
- `2_GAMM_Significance.R`: `R` script to assess the significance of the
  smooth *Hemisphere×Age* interaction by testing for a difference in the
  smooth term of *Age* between hemispheres. The script applies a False
  Discovery Rate correction to control for the number of tests
  conducted.
- `3_PAM_Clustering.R`: `R` script to classify regions in the
  **Language-and-Memory** network that demonstrate a significant
  *Hemisphere×Age* interaction, based on their functional asymmetry
  skewness profiles. This script also allows to compute the intersection
  point between the two average clusters curves.
- `4_CCA_BrainCognitionAssociation.R`: `R` script to proceed with the
  Canonical Correlation Analysis to assess brain–behavior Associations.

Note that the `R` scripts also contains the code **to reproduce the
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
  - The file `Atlas/language_memory_atlas.txt` contains a full
    description of each Language-and-Memory areas. The first column
    *abbreviation* corresponds to the abbreviation of each region. The
    second column *region* is the anatomical labels of each regions. The
    third column *hemisphere* is the location of the region in the left
    or right hemisphere. The fourth column *index* corresponds to the
    index of each region that is used in the `NIfTI` file. The fifth
    column *vol_vox* corresponds to the number of voxels in each region.
    Finally, the MNI coordinate (columns *xmm*, *ymm*, and *zmm*) of
    each regions centroid is provided.

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

------------------------------------------------------------------------

## Questions

Please contact me (Loïc Labache) at: <loic.labache@yale.edu> and/or
<loic.labache@ensc.fr>, or Élise Roger at: <elise.roger@umontreal.ca>
