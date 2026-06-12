<img width="303" height="450" alt="image" src="https://github.com/user-attachments/assets/bd130cbb-ae88-4c15-9c08-a7eff09698c1" />


Hugo de Vries (1848-1935)  

---

# Paper Title.
This repository contains the data and analytical pipeline for Weir and Phillimore (2026). Title. Citation.

The paper explores the effect of asycnrhony on the performance of polyphagous spring caterpillars, across a food-web matrix of six consumer and eight host-plant taxa. Caterpillar hatch-time/diapause break was manipulated to induce various levels of asynchrony relative to bud-burst time on the host trees. Caterpillars were reared in captivity to assay variation in performance among trophic interactions at different levels of asynchrony. We consider these results in light of the resilience of these species under climate change and the extent to which such trophic interactions are buffered against asynchrony. 

## Data
Data ('all_synchrony_data.csv') consists of caterpillar (species) performance on different host-plant species (host) with simulated late-hatching asynchrony (day, with 0 as perfect synchrony, to a maximum of 65 days late-hatching). Performance was measured as: survival (survp, 0/1, survival from the start of the experiment to successful pupation); pupal mass (masspm, mg, surviving pupae measured during summer diapause); and development time (durationl, calendar days, time from the start of the experiment to pupation). Caterpillars were reared in uniquely numbered 'rearing cultures' (culu) of 10 individuals. of Female moths are typically larger than males, particularly species which do not feed as adults and must therefore carry a complete burden of eggs. Pupaae can be sexed (sex) by inspecting genital pores, shape, and size. The experiment was carried out over two years (2020 and 2021), with year one (2020) assaying a smaller subset of host-caterpillar interactions. Because different tree species have different spring phenologies (i.e., come into leaf at slightly different times) it was not possible to start all experimental treatments at exactly the same time. For example, Hawthorn comes into leaf consistently before Oak. Therefore, different host treatments were established on different dates (date_est), within a ~5 day time-frame. 

## Workflow pipeline
This pipeline contains the analytical milestones needed to reproduce the results of this study, including data visualisation (e.g., predictions and plotting).

- **Step 1** `1_modelling.R` : Estimating the caterpillar peak date and height (i.e. peak caterpillar biomass) in each site by year combination across available data. 
- **Step 2** `2_predictions.R` : Extract interpolated MetOffice temperature data for spring in each site by year combination relevant to our dataset.
- **Step 3** `3_plotting_code.R` : Model variation in spring temperature across all sites over time.
- **Step 4** `4_variance_comp_estimation.R` : Collate bird breeding performance data across sites and years, pairing data with relevant caterpillar abundance data and estimating trophic mismatch at the nest-level (i.e. difference in timing between peak chick demand within a nest and date of peak caterpillar abundance).
- **Step 5** `5_slope_extraction_and_plotting.R` : Model the relationship between temperature and trophic mismatch in three target bird species.

## Rights
The authors reserve all rights to the data, which cannot be utilised in whole or in part without their prior consent.

## Contact
JamieCWeir@outlook.com
