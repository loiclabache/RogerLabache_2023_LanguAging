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

The `Script` folder contains 2 sub-folders: `Analysis` and
`Visualization`.

The `Analysis` folder contains the scripts to reproduce the results
presented in the manuscript.

- `X.R`: `R` script to \[…\].

The `Visualization` folder contains `R` files (`FigX_script.R`) used to
generate each figures included in the manuscript. Each script
corresponds to a figure or a panel. The brain renderings in the paper
require a customized version of [Surf
Ice](https://www.nitrc.org/projects/surfice/) that we will be happy to
share on demand.

------------------------------------------------------------------------

## Atlas Used

The atlas used in the paper is available in the `Atlas` folder.

<!-- -   **SENSAAS** provide an atlas in standardized MNI volume space of 32 sentence-related areas based on a 3-step method combining the analysis of *activation and asymmetry during multiple language tasks* with hierarchical clustering of resting-state connectivity and graph analyses. The temporal correlations at rest between these 32 regions made it possible to detect their belonging to 3 networks. Among these networks, one, *including 18 regions*, contains the essential language areas (**SENT_CORE** network), *i.e.* those whose **lesion would cause an alteration in the understanding of speech**. Full description of the language atlas can be found there: [SENSAAS](https://github.com/loiclabache/SENSAAS_brainAtlas), and the related paper there: [Labache, L., et al. 2019](https://doi.org/10.1007/s00429-018-1810-2). -->
<!--     -  The *volumetric* (in the MNI ICBM 152 space) and *area* (32k_fs_LR space) atlas are available in the `Atlas/SENSAAS` folder. The sub-folder `Volumetric` contains the volumetric SENSAAS atlas: `SENSAAS_MNI_ICBM_152_2mm.nii`, and a CSV file containing a full description of each language areas: `SENSAAS_description.csv`. The sub-folder `Area` contains the area SENSAAS atlas in the left (`S1200_binarySentCore_L_surface.shape.gii`) and right hemisphere (`S1200_binarySentCore_R_surface.shape.gii`), as well as the hub atlas (*i.e.* regions STS3, STS4 and F3t only) in the left (`S1200_binaryHubsSentCore_L_surface.shape.gii`) and right hemisphere (`S1200_binaryHubsSentCore_R_surface.shape.gii`). -->
<!--         -  Briefly, the hub language network atlas corresponded to the inferior frontal gyrus (Broca’s area, F3t) and to the posterior aspect  of the superior temporal sulcus (corresponding to Wernicke’s area, STS3 and STS4). -->
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
