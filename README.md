# HSC_study

Could we reconstruct/interpolate sparsely sampled data (Subfig (a)) by adding in domain knowledge and mathematics?
See Subfig (b-d) for simulated sampling at higher frequencies.

![Alt text](figs/cartoons/sample_gaps.png?raw=true "Optional Title")


Situation:
    - Hematopoietic Stem Cells (~10^4-10^5) keep us alive by generating 10^9 
    blood cells daily (white/red blood cells, platelet). 
    Given that total renewal rate in peripheral blood is 10^11 cells per day, 
    each stem cell is responsible for 10^6~10^7 cell production daily. 
    - Important to understand HSC behaviors at single-cell level 
    (e.g. for transplant), but extremely hard. 
        - Hidden (bone marrow)
        - Complex (multi lineage)

Target: 
    - Have data matrix of about 1000 rows, 9 columns from peripheral blood samples. 
    Each row represents progenies of a unique stem cell. 
    Each col represents sampling at a specific time point (2, 3, 6, 11... months
    after transplantation).
    - Goal: are there different behaviors among different stem cells?

Action:
    - First look
        - Biologist: found signal of 7 pattern groups
        - Mathematician: Seem from noise
    - Null hypothesis -> Modeling:
        - Assume static. Does each row come from the same distribution? 
            - Eye balling: No
            - Test for mean: One-way ANOVA. 
            - Test for variance: Levene's
        - Does every stem cell starts the same?
            - Starts same -> Equilibrium
            - New Null hypothesis: A part of stem cells' behaviors can be explained by noise 
                - (1) What does noise do? 
                - (2) What fraction?
    - What does noise do?
        - Forward model (multistage): 
            - HSC birth death noise -> stem cell distribution in number (static) 
                in the bone marrow
            - HSC fire noise + progenies birth death noise -> cells in number 
            over time from bone marrow to peripheral blood
            - TODO
            - Result: a distribution of cells in the simulated sample
        - Comparing simulated sample with real data
            - Gaussian distribution. Also most biologically relevant: 
                Mean (number/activity of each stem cells) and 
                variance (correlates with mean; also relates to all noises)
            - However, cannot fit either individually, since no idea about the
                static distribution of stem cells
            - Fitting the functional of Coefficient of Variation
            - If null hypothesis is true, and all parameters are right, should
            be fitting the functional of Coefficient of Variation at any mean
    - Target is not testing if all clones have the same CV, but if it is possible
        to consistently model CVs of different stem cells by noise. 
        - By tuning noise from different stage, most sensitive to L. 
        - By tuning L, CV varies.  
        - Difference between simulated CV and experiment CV:
            - Binning, averaging, and differing
        - Best fixing best parameters, simulating the model, found that >90% 
            points in the 95% interval (TODO...)
            
        
            