Authors' response to: “PREreview of Predicting Relative Populations of Protein Conformations without a Physics Engine Using AlphaFold2” (10.5281/zenodo.8235542, 2023), by Ashraya Ravikumar and Sonya Lee with input from other Fraser Lab members at UCSF.

# **Preface:**

We thank the reviewers for their interest in our work and for the thoughtful, clearly-organized, and fair feedback. Besides increasing the reproducibility and rigor of our work, the points raised by the reviewers led us to explore previously-overlooked aspects of our study, revealing new insights with significant implications for future studies. In a field that moves at a break-neck pace, this form of feedback is invaluable for optimizing the rate of discovery and the accuracy of the findings. 

In this response, we address each of the major points raised in the review.

We reference a set of figures and a table that were prepared to address points raised by the reviewers. Unfortunately, these figures cannot be shared in text-only forms, but we have hosted them in a public repository, accessible at this link: <https://github.com/GMdSilva/10.5281-zenodo-8235542_response>, which also includes all the raw data referenced in this response.

# **Response:**

## **Major points:**

**1.0. The MSA subsampling approach that the authors have adapted in this work has been used by others previously (as cited by the authors themselves), albeit with some modifications. So it is important to see if the existing methodologies, for instance the DBSCAN based clustering and MSA subsampling by Wayment-Steele et al., are able to predict these relative state populations of variants.**

**R:** We agree that these comparisons would be important, especially in the context of benchmarking and understanding the strengths and weaknesses of each approach. Although we thought of running a few comparisons with the DBSCAN method, we ultimately decided to focus our study on showcasing the results of our approach in the context of the two test cases instead of as a comparison with previous methods. We are very curious to see how our results line up with results from alternative methods and are happy to assist other groups with further studies by sharing our data or running predictions.

**1.1. Also, the optimization of max\_seq:extra\_seq requires quite a bit of pre-existing experimental information. How is this method to be applied for a relatively new system? The authors could also provide some guidelines on how the max\_seq:extra\_seq numbers to be sampled are c0hosen and in general comment about the hyper-parameter space in their approach and how it compares to other schemes/approaches.**

**R:** This is an excellent point, which could have been more thoroughly addressed in our current manuscript and will be in subsequent versions. While our approach is **best suited for the high-throughput exploration** of the effects of changes in previously-explored systems, we argue that our method can also stud mostly unknown (in the context of structural or dynamics data) systems. 

To accomplish this, a user should first identify structural elements with significant conformational variation within the predicted ensemble (such as by calculating the root mean square fluctuation of alpha carbons within the ensemble when aligned with the prediction of the highest-ranking pLDDT, or more easily by evaluating per-residue pLDDT values, as lower pLDDT values are usually correlated with increased dynamics, see Figure 1’A for an example of this correlation.)

Following that, the user could plot the RMSD of that specific region for each prediction with respect to the highest-ranking (by average pLDDT) reference, and either analyze the distribution visually or via statistical methods.

This approach works exceptionally well for Abl1. The A-loop is the structural element with the largest amplitude of conformational changes in the ensemble, as evidenced in Figure 1’A. By plotting the A-Loop backbone RMSD with respect to the top-ranked prediction for each prediction within the Abl1 ensemble, we clearly identify patterns across the parameter optimization process. 

Specifically, we see that the distributions of predictions as a function of RMSD are similar and cluster around the ground state for the most extreme parameter values (32:64 and 2048:4096), but become increasingly wider and more diverse for intermediate parameter values(256:512 and 512:1024), where we observe an increased number of samples around 14 A. For Abl1, we already know the implications of these observations in regards to its structure and function, but in the absence of previous information, it would be reasonable to infer that the cluster at \~14 A represents a metastable conformation of lower population than the ground state and that the predictions that lie in-between them (more diverse in 256:512 and 512:1024) could represent low-occupancy intermediate states. Consequently, in systems for which there is limited information, choosing parameters that lead to the widest distribution of configurations and the largest number of intermediate predictions between two or more states of interest should lead to optimal results, without the need for experimental data. 

![](https://lh3.googleusercontent.com/BzVkXCX4GUHgOZSB0yldGSTWPPqyWrveaF-JWNDatsDZjJkc9uz45M17r3_9-MComqjYjOqEglaQvVaSUVCwA51qgwQMrrxDiyZoXGWqmApq5xdt457rRn8p7nhxs2B4nNG41UIDDFVAkYbZ9_ESzyk)

Figure 1’. Results of subsampling AF2 for the Abl1 kinase core prediction. A: Average per-residue pLDDT and per-residue root mean square fluctuations (when aligned with the top-ranked prediction) across a predicted ensemble are correlated: residues with lower-than-average pLDDT tend to fluctuate significantly. Data points are colored according to residue positions. As residue positions are 0-indexed in this plot, the activation loop roughly corresponds to residue positions 153-173. B: Distribution of activation loop RMSDs with respect to the top-ranked prediction for a range _max\_msa:extra\_msa_ parameter sets. Each dot indicates one prediction within the ensemble. The grey box indicates the cluster corresponding to ground-like predictions, the purple box indicates transitory/intermediate conformations and the black box indicates I2-like predictions. 

**2.0. Apart from the large change in A-loop from active to inactive state in Abl kinase, the other important structure change involves the 𝛼C helix moving out (as shown in Reference 22 cited in the preprint). The authors have not discussed this aspect. Does the AF2 subsampled ensemble reflect the change in the helix position?**

**R:** We agree with the reviewers that looking at aC helix conformations is important, and regret that we neglected it when writing the original manuscript. In response to this suggestion, we calculated the RMSD of the aC Helix backbone (residues 275-295) with respect to both the Ground and I2 state references (6XR6 and 6XRG, respectively) for each ensemble of wild-type Abl1 predictions at different subsampling labels, in the same way we did for the activation loop and plotted the results in figure 2’A. 

Importantly, we see a distribution of results that is strikingly similar to the activation loop results (Figure 1’A). When plotting these variables against each other (Figure 2’B), we observe clear clustering corresponding to Ground-like and I2-like states, suggesting that AF2 is not predicting just the activation loop conformational changes, but the aC-helix movement as well. This pattern holds for mutant forms of Abl1 (mutations that decrease the ground state population lead to predictions with increased aC-Helix motions and an increased frequency of predictions different than the ground state;an opposite effect is observed for those that increase the ground state population, as summarized in Figure 2’B). In summary, our approach accurately identifies that aC-Helix and A-loop motions are correlated in our ensemble of predictions, and can be used individually or in tandem to distinguish between conformational states.

![](https://lh3.googleusercontent.com/MhKcJ-CgzEXrg1dlJ5fjLAZK6SKfp5WTK9jZ9H-3EF7NUxF7X2QJ3OUB0xiPW-xMlwlK2C_edr9rGCXDBu__ORApnRUe3fUlUvk6-SQq_W97DW_Xz-_9_7a0ONRDz7rTWTgLzP5d5KjbDcveRIFR7sQ)

Figure 2’. Distribution of aC-Helix conformations within the predicted AF2 ensembles. A: Distribution of aC-Helix RMSDs with respect to the top-ranked prediction for a range _max\_msa:extra\_msa_ parameter sets. Each dot indicates one prediction within the ensemble. B: Correlation between aC-Helix and A-Loop RMSDs compared to the Ground state reference (PDB 6XR6) or the I2 state reference (PDB 6XRG). Data points are colored according to their average predicted local distance difference test score (pLDDT) as calculated by AF2, which is a metric of the confidence/accuracy of AF2 predictions

**2.1. The snapshots shown from enhanced MD does not seem to show this change either (upon visual examination of the snapshots shown in the figures). Hence, the biological relevance of the MD simulation becomes questionable.** 

**R:** We thank the reviewers for bringing this to our attention. We do observe the aC-Helix shifting from its initial position in our simulations, but due to a labeling error while making Figure 4 of the original manuscript, we miscommunicated the scope of our results, leading to the impression that there is no change in the aC-Helix position (which is also not helped by the camera angle of the structural renderings we included in the figure, which makes it hard to discern changes in the helix position).

As described in the supporting information, we ran two separate enhanced-sampling simulations using the Weighted Ensemble method to determine the timescales of the Ground to  I2 transition. This led to two representative trajectories: Simulation I spanning 13 ns corresponding to the Ground to I1 (DFG flip) transition, and Simulation II spanning 24 ns corresponding to the I1 to I2 (A-Loop collapse) transition. 

In Figure 4, we took snapshots from the I1 to I2 state simulation **(in contrast with the Ground to I2 simulation, as erroneously labeled in the paper)**, since the activation loop transition happens within that interval. At 0 ns for simulation II, the aC-Helix is already slightly shifted with respect to its ground state position, as some of the aC-Helix conformational change happens in the Ground to I1 transition. Since the shown snapshots are all from the I1 to I2 transition, Figure 4 in the manuscript doesn’t accurately communicate aC-Helix conformational changes.

We will remedy this oversight by adjusting the language in the figure and related discussions in the text. Further, for the sake of clarity, we demonstrate in Figure 3 that our collated trajectory (II appended to the end of I) captures the aC-Helix moving away from the kinase core. We have also updated the original plots to include the entire 37 ns simulation (Figure 3’B). 

Lastly, we also noticed an error regarding the time points at which snapshots were collected and adjusted the vertical bars to represent the actual time points from which they were taken. We will rework Figure 4 with the correct timestamp labels for each snapshot in a future version of the manuscript.

For the sake of transparency, our collated trajectory (solvent and ions removed) is available at <https://github.com/GMdSilva/10.5281-zenodo-8235542_response/blob/main/abl1_wt_gr_to_i1_to_i2_traj.dcd>, and the topology at: https\://github.com/GMdSilva/10.5281-zenodo-8235542\_response/blob/main/abl1\_wt\_gr\_to\_i1\_to\_i2\_topology.pdb.

![](https://lh3.googleusercontent.com/LTMO6GPfYSi6G8dUN3zG4On---3g3g9bRgzk70lbaZNLD1iBNhbmMzwfjVTKUsZ6UxQUjNDYbIV3aueUJ9rC468WtqKtt1AJJvgoGvA7KBZzfhBucbL8QOSfmBDwP9ilkdxGRfaQqoSIDiP6JRpcvAw)

Figure 3’. Results of enhanced-sampling molecular dynamics (MD) simulations of the Abl1 Ground to I2 transition. A: Relevant Abl1 kinase core elements referenced in this study. B: Overlay of the I2 reference (PDB 6XRG) and three representative snapshots of the simulation spanning the transition showcasing differences in the aC-Helix and A-Loop positions. C: Evolution of simulation observables as a function of time.

**3. The authors haven’t performed statistical analyses on the RMSD comparisons or the CSP comparisons of GMCSF to claim the differences to be significant or not. For example, the authors say their approach has worked “as the range of the distribution of RMSDs of residues 80-90 and 110-125 is significantly larger for most of the mutations tested at both of these sites”. What is this distribution of RMSD compared against? Are these differences statistically significant?**

**R:** We thank the reviewers for bringing this to our attention and realize in hindsight that the text could be clearer regarding this comparison. We are comparing to the wild-type prediction that looks the most like the ground state, defined as the conformation occupied by PDB ID: 1CSG.

We used the Kruskal-Wallis H-test to measure the significance of the differences we observe in the GMCSF distributions and report them in Table ‘1. This test was chosen because it is unlikely that RMSD values will follow a normal distribution considering the nature of energy landscapes (which are more likely to exhibit multimodal distributions). These results line up with our previous observations: the H83Y predictions led to statistically insignificant differences, while the H87 mutations lead to differences that are orders of magnitude less significant than H83 or H15 mutations. The H83N mutation also led to changes that are not as statistically significant as those induced by H83R, matching the observed NMR results. 

Table '1: Results of the Kruskal-Wallis H-test for the GMCSF ensemble predictions. The distributions for each observable described in "test" for each ensemble prediction described in "trial" were compared to the distribution of the same observable for the wild-type (GMCSF) for the hypothesis test.
| p_value 	| h_stat 	| trial 	| test 	| sample_size 	|
|---	|---	|---	|---	|---	|
| 1 	| 0 	| GMCSF 	| Distance 15/CA to 83/CA 	| 480 	|
| 0.070398 	| 3.273716 	| GMCSF H15Y 	| Distance 15/CA to 83/CA 	| 480 	|
| 0.000289 	| 13.14039 	| GMCSF H15R 	| Distance 15/CA to 83/CA 	| 480 	|
| 0.013501 	| 6.102222 	| GMCSF H15N 	| Distance 15/CA to 83/CA 	| 480 	|
| 0.047731 	| 3.919421 	| GMCSF H83Y 	| Distance 15/CA to 83/CA 	| 480 	|
| 0.000392 	| 12.57116 	| GMCSF H83R 	| Distance 15/CA to 83/CA 	| 480 	|
| 0.027057 	| 4.887163 	| GMCSF H83N 	| Distance 15/CA to 83/CA 	| 480 	|
| 0.553995 	| 0.350208 	| GMCSF H87Y 	| Distance 15/CA to 83/CA 	| 480 	|
| 0.5115 	| 0.430996 	| GMCSF H87R 	| Distance 15/CA to 83/CA 	| 480 	|
| 0.546459 	| 0.363699 	| GMCSF H87N 	| Distance 15/CA to 83/CA 	| 480 	|
| 1 	| 0 	| GMCSF 	| Backbone RMSD vs. Ref. 	| 480 	|
| 0.02232 	| 5.220667 	| GMCSF H15Y 	| Backbone RMSD vs. Ref. 	| 480 	|
| 9.55E-06 	| 19.59876 	| GMCSF H15R 	| Backbone RMSD vs. Ref. 	| 480 	|
| 3.55E-05 	| 17.0962 	| GMCSF H15N 	| Backbone RMSD vs. Ref. 	| 480 	|
| 0.205222 	| 1.60482 	| GMCSF H83Y 	| Backbone RMSD vs. Ref. 	| 480 	|
| 0.006538 	| 7.39582 	| GMCSF H83R 	| Backbone RMSD vs. Ref. 	| 480 	|
| 0.021759 	| 5.264923 	| GMCSF H83N 	| Backbone RMSD vs. Ref. 	| 480 	|
| 0.779695 	| 0.078241 	| GMCSF H87Y 	| Backbone RMSD vs. Ref. 	| 480 	|
| 0.246255 	| 1.344427 	| GMCSF H87R 	| Backbone RMSD vs. Ref. 	| 480 	|
| 0.311422 	| 1.024635 	| GMCSF H87N 	| Backbone RMSD vs. Ref. 	| 480 	|
| 1 	| 0 	| GMCSF 	| Residues 80-90 Backbone RMSD vs. Ref. 	| 480 	|
| 2.99E-26 	| 112.3543 	| GMCSF H15Y 	| Residues 80-90 Backbone RMSD vs. Ref. 	| 480 	|
| 4.18E-47 	| 207.7826 	| GMCSF H15R 	| Residues 80-90 Backbone RMSD vs. Ref. 	| 480 	|
| 7.70E-27 	| 115.0427 	| GMCSF H15N 	| Residues 80-90 Backbone RMSD vs. Ref. 	| 480 	|
| 6.27E-13 	| 51.75958 	| GMCSF H83Y 	| Residues 80-90 Backbone RMSD vs. Ref. 	| 480 	|
| 2.59E-37 	| 162.9281 	| GMCSF H83R 	| Residues 80-90 Backbone RMSD vs. Ref. 	| 480 	|
| 7.10E-09 	| 33.50829 	| GMCSF H83N 	| Residues 80-90 Backbone RMSD vs. Ref. 	| 480 	|
| 0.029922 	| 4.713746 	| GMCSF H87Y 	| Residues 80-90 Backbone RMSD vs. Ref. 	| 480 	|
| 0.703822 	| 0.144527 	| GMCSF H87R 	| Residues 80-90 Backbone RMSD vs. Ref. 	| 480 	|
| 2.33E-08 	| 31.20188 	| GMCSF H87N 	| Residues 80-90 Backbone RMSD vs. Ref. 	| 480 	|
| 1 	| 0 	| GMCSF 	| Residues 110-125 Backbone RMSD vs. Ref. 	| 480 	|
| 0.243007 	| 1.363067 	| GMCSF H15Y 	| Residues 110-125 Backbone RMSD vs. Ref. 	| 480 	|
| 0.002709 	| 8.993774 	| GMCSF H15R 	| Residues 110-125 Backbone RMSD vs. Ref. 	| 480 	|
| 0.012352 	| 6.259648 	| GMCSF H15N 	| Residues 110-125 Backbone RMSD vs. Ref. 	| 480 	|
| 0.716998 	| 0.131387 	| GMCSF H83Y 	| Residues 110-125 Backbone RMSD vs. Ref. 	| 480 	|
| 5.86E-05 	| 16.14757 	| GMCSF H83R 	| Residues 110-125 Backbone RMSD vs. Ref. 	| 480 	|
| 0.574285 	| 0.315566 	| GMCSF H83N 	| Residues 110-125 Backbone RMSD vs. Ref. 	| 480 	|
| 0.547621 	| 0.361596 	| GMCSF H87Y 	| Residues 110-125 Backbone RMSD vs. Ref. 	| 480 	|
| 0.716303 	| 0.132063 	| GMCSF H87R 	| Residues 110-125 Backbone RMSD vs. Ref. 	| 480 	|
| 0.826416 	| 0.048092 	| GMCSF H87N 	| Residues 110-125 Backbone RMSD vs. Ref. 	| 480 	|
| 1 	| 0 	| GMCSF 	| Residues 45-55 Backbone RMSD vs. Ref. 	| 480 	|
| 3.74E-06 	| 21.39144 	| GMCSF H15Y 	| Residues 45-55 Backbone RMSD vs. Ref. 	| 480 	|
| 2.03E-11 	| 44.93757 	| GMCSF H15R 	| Residues 45-55 Backbone RMSD vs. Ref. 	| 480 	|
| 1.41E-14 	| 59.21387 	| GMCSF H15N 	| Residues 45-55 Backbone RMSD vs. Ref. 	| 480 	|
| 7.11E-07 	| 24.58403 	| GMCSF H83Y 	| Residues 45-55 Backbone RMSD vs. Ref. 	| 480 	|
| 1.69E-13 	| 54.33155 	| GMCSF H83R 	| Residues 45-55 Backbone RMSD vs. Ref. 	| 480 	|
| 1.03E-12 	| 50.78773 	| GMCSF H83N 	| Residues 45-55 Backbone RMSD vs. Ref. 	| 480 	|
| 0.012865 	| 6.187631 	| GMCSF H87Y 	| Residues 45-55 Backbone RMSD vs. Ref. 	| 480 	|
| 0.192417 	| 1.699012 	| GMCSF H87R 	| Residues 45-55 Backbone RMSD vs. Ref. 	| 480 	|
| 0.061018 	| 3.509485 	| GMCSF H87N 	| Residues 45-55 Backbone RMSD vs. Ref. 	| 480 	|
| 1 	| 0 	| GMCSF 	| Residues 10-20 Backbone RMSD vs. Ref. 	| 480 	|
| 0.148456 	| 2.08804 	| GMCSF H15Y 	| Residues 10-20 Backbone RMSD vs. Ref. 	| 480 	|
| 0.024358 	| 5.068952 	| GMCSF H15R 	| Residues 10-20 Backbone RMSD vs. Ref. 	| 480 	|
| 0.022179 	| 5.231672 	| GMCSF H15N 	| Residues 10-20 Backbone RMSD vs. Ref. 	| 480 	|
| 0.498708 	| 0.457685 	| GMCSF H83Y 	| Residues 10-20 Backbone RMSD vs. Ref. 	| 480 	|
| 0.015353 	| 5.875485 	| GMCSF H83R 	| Residues 10-20 Backbone RMSD vs. Ref. 	| 480 	|
| 0.162539 	| 1.950442 	| GMCSF H83N 	| Residues 10-20 Backbone RMSD vs. Ref. 	| 480 	|
| 0.611471 	| 0.258039 	| GMCSF H87Y 	| Residues 10-20 Backbone RMSD vs. Ref. 	| 480 	|
| 0.531769 	| 0.39101 	| GMCSF H87R 	| Residues 10-20 Backbone RMSD vs. Ref. 	| 480 	|
| 0.954516 	| 0.003253 	| GMCSF H87N 	| Residues 10-20 Backbone RMSD vs. Ref. 	| 480 	|

The raw data used in the Kruskal-Wallis H-test is available at: https\://github.com/GMdSilva/10.5281-zenodo-8235542\_response/blob/main/gmcsf\_histidine\_triad\_predictions\_maxseq4extraseq8.csv

**4. Given that GMCSF has very limited sequence data in MSA to start with, does MSA subsampling actually help? The authors could try doing predictions using the traditional AF2 pipeline and compare those distributions against their approach.**

**R:** Although this is not made clear in our results, we did run GMCSF without subsampling, and show the results in Figure S8. Specifically, since the length of the GMCSF MSA is only about 120 sequences, any of the subsampling runs with _max\_seq_ and _extra\_seq_ values greater than 120 are, in essence, using the traditional AF2 pipeline (with drop-offs enabled) as there is no subsampling employed. We still included these results for the sake of transparency, but in a future version of this manuscript, we will be clearer about this distinction and what it means with respect to the effects of subsampling.

## **Minor points:**

**1. Although the authors are right in looking for only the ground and I2 states in Abl kinase predictions, it will be interesting to explore if there were any predictions that matched the I1 state and if not, to speculate why more extensively**

**R.** Originally, when conceptualizing this study, we did look at I1 state similarities, but the results were inconclusive and hard to make sense of, so we dropped that from the scope of the work. However, when looking at the data now, it is clear that AF2 is predicting I1-like elements in each ensemble, although the resolution of those predictions is significantly lower than that for the I2 predictions. 

Specifically, the I1 state is mostly defined by a flip of the position of the side chains of D381 and F382, part of the DFG motif. In the ground state, D381 points towards the active site and F382 points towards the aC-Helix, allowing productive phosphorylation. In the I1 state, these side-chain positions are inverted, and F382 occludes the active site. Previous studies have shown that two dihedrals clearly distinguish between the four potential permutations of the DFG motif side-chain orientations: D in/F in (DFG in), D out/F out (DFG out), D out/F in (DFG up), and D in/F out (DFG down). The dihedrals are defined as follows: 

Dihedral 1: A380CB/A380CA/D381CA/D381/CG

Dihedral 2: A380CB/A380CA/F382CA/F382/CG

In the Ground to I1 transition, dihedrals 1 and 2 invert. The DFG up and DFG down configurations are short-lived in Abl1 and are usually only observed in MD simulations (where femtosecond resolution is feasible). 

When plotting dihedrals 1 and 2 for our ensemble of predictions (Figure 4’), we observe a curious behavior for dihedral 1: **the distribution is strongly bimodal,** with two clearly-defined clusters, one at the expected position for dihedral 1 in the Ground state, and one much less populated cluster situated around 20-40 degrees. Although the latter is nowhere near the values expected of the dihedral 1 position in I2 (about -120 to -140), the fact that there is clustering at this position and that the cluster populations mimic that of the Ground and I1/I2 populations suggests that AF2 is picking up some amount of dynamical signal that corresponds to the Ground to I1 transition, although the signal is not strong enough and/or AF2’s resolution is not high enough to correctly replicate the expected result. 

Even more interestingly, the population of this D “out” cluster is increased in mutations known to dislocate the conformational equilibrium away from the ground state and decreased in mutations known to have the opposite effect (although with very low statistical power), similar to what is observed in the A-Loop or aC-Helix conformational predictions (Figure 5’). This suggests that the dihedral 1 position is correlated with A-loop and aC-Helix conformations in the predictions, which is in agreement with the defining elements of the I1 and I2 states.

Finally, analyzing the distribution of dihedral 2 leads to significantly murkier results, as the swarm plot takes a comet form with a single largely-populated cluster at the expected values for the F-in conformation, and a few predictions that deviate from it up to -160 (which is also significantly far from the expected F-out values of 40 to 60). Interestingly, the population of the “outlier” predictions negatively correlates with the ground state population for the Abl1 allelic series. These results hint that there could be more coevolutional signal for the D381 flip than for the F382 flip, which corroborates MD simulation results that show that it is mostly D381 interactions that drive the entire DFG flip. 

![](https://lh4.googleusercontent.com/BGy5PaVa8aVuZhR6WI5U0zj5lTYfhB_BVnbyvXHiLeNltqeDhvG2W2Z936Ye8K5lgs3wMhDkmzu9cpvP0F3tEYE3G5ujame1Cgfe_1NrjGwWAv1JH1yYGItJW_-BhXYrdbWH80SGJrTm-jz6z1yMmJc)

Figure 4’. A: Correlation between dihedral 1 and dihedral 2 values within three representative predicted ensembles. Data points are colored according to their average predicted local distance difference test score (pLDDT) as calculated by AF2, which is a metric of the confidence/accuracy of AF2 predictions. 

![](https://lh4.googleusercontent.com/82x8ceoOINcU7SbklT3yykjYcoNKv-dNe76j2eGGVTNnho62OCXkHaWocw4QjhNsQu0cwm6SrOFr5eLaHzRaY2DtS9IOC-qNXs1wKk75xEsuJY-sLM94l1mMrpNPQkmBZK7IGviqJB37Am193cEoDMM)

Figure 5’. Distribution of dihedral 1 and dihedral 2 values within ensembles predicted with subsampled AF2. A: Results pertaining to the Abl1 evolutionary line test. B: Results pertaining to the Abl1 allelic series test. C: Percentage of Dihedral 1 prediction with values lower than 40 in the evolutionary line test. D: Percentage of Dihedral 1 predictions with values lower than 40 in the allelic series test.

**2. The data on some of the max\_seq:extra\_seq optimizations discussed for Abl kinase is missing. For example, 512:8 or 8:1024**

**R:** The text is supposed to say 32:512 and 512:32 (matching Figure 5), but we neglected to change placeholder values for the data that were included in a previous figure.

**3. There is no citation provided for the single and double mutants whose relative ground state populations were tested for Abl Kinase.**

**R:** The relevant citation is #22: Xie, T., Saleh, T., Rossi, P. & Kalodimos, C. G. Conformational states dynamically populated by a kinase determine its function. _Science_ 370 (Oct. 2020). We reference it in this specific context in the caption of Figure 6:

“(B) A series of point mutations in wild-type Abl1 are known to increase or decrease the relative population of the enzyme’s active (ground) state **\[22].”**

We will add a reference to the main text and to Figure S2 for improved clarity.

**4. The nature of these mutations on Abl Kinase is not discussed. Are some of these mutations pathogenic or drug-resistant? It will be interesting to correlate the nature of mutation with its structural effects. The authors could provide more introduction of how these mutations were identified and add more discussion on trends.**

**R:** Three of the mutations that increase the ground-state population are associated with type II inhibitor resistance in Abl1 (E255K, T315I, E255K+T315I). All mutations were chosen due to the extensive data regarding their effects on the populations of the Ground or I2 state (Kalodimos et al., _Science_ 2020; Kern et al, _PNAS_ 2020). We agree with the reviewers that correlating the nature of mutations with structural effects has important implications, especially in the case of the mutations that lead to drug resistance. Our running hypothesis is that reduced affinity to type II inhibitors is caused by the further stabilization of the Ground state caused by the mutations E255K, T315I, and E255K+T315I, as type II inhibitors have a strong preference for the I2 form of Abl1. We will discuss these implications in more detail in a revised version of the manuscript. 

**5. What is the rationale behind choosing the PCs mentioned by the authors for Abl kinase enhanced sampling?**

**R:** 

Ground to I1 simulation:

PC1:  Based on the observation by the D. E. Shaw Research group that protonated D381 forms a transient hydrogen bond with V299 in an unbiased DFG flip (Shaw et al., _PNAS 2009_).

PC2.  Based on a metadynamics collective variable used in a DFG flip study of the Aurora Kinase B (Li et al, _AAPS J_. 2019). This PC2 allows for differentiation between DFG-out and DFG-up (PC1 is not sufficient for that, although it is sufficient to sample the flip).

I1 to I2 simulation:

PC1: Chosen based on visual inspection of two atoms defining a distance that should increase during the I1 to I2 transition.

PC2: Chosen based on visual inspection of two atoms defining a distance that should decrease during the I1 to I2 transition.

**6. Why have the authors not shown the RMSD distribution of Distance 2 in Figure 4C?**

**R:** We made the decision to not show the distribution of the Distance 2 deltas in Figure 4C because the resulting plot was cluttered when using the same figure aesthetics as the other distributions (Figure 6’), and we thought it to be largely redundant considering the other elements in the figure. We included the Distance 2 delta data in its entirety in Figures S2-5. 

![](https://lh4.googleusercontent.com/CipJjWMgUc0QKFwKUHSGqG6epU7kGU4sBmO1wbbxM2KE8qC4V-HpZi9ICRRYHMPq5AyWN9Z4Hqex-hD2gRzdgm80tDbDHcuMQIpwJVF4FZ__OaY69bPaGSkP9zvsz_nUeLUJjfZIFpfQVVy3U6CzM_o)

Figure 6’: Distribution of the Distance 2 delta observable in the analysis presented in Figure 4 of the original manuscript.

**7. How were the mutations on the histidine triad of GMCSF chosen?\
Sebollela et al. 2005 (**[**https://pubmed.ncbi.nlm.nih...**](https://disq.us/url?url=https%3A%2F%2Fpubmed.ncbi.nlm.nih.gov%2F16027123%2F%29%3A7XkDYyJBOndfuaRuiP3tq_nsZZc\&cuid=2634513)**, which is not cited in this paper specifically but cited in one of the papers (Cui et al. 2020 -** [**https://doi.org/10.1021/acs...**](https://disq.us/url?url=https%3A%2F%2Fdoi.org%2F10.1021%2Facs.biochem.0c00538%29%3AjLjBlPOl2RaJcfMjJDKDVoqMX1E\&cuid=2634513) **that this paper cites, substitutes H15 with alanine to demonstrate a decrease in heparin affinity**

**R:** The mutations were chosen based on sequence conservation data of GMCSF. Specifically, positions 15, 83, and 87 are mostly conserved as histidine, but some species tolerate mutations to asparagine, arginine, or tyrosine. We did not test the alanine substitution as that leads to a protein with poor expression.

**8. For the GMCSF system, do the authors see a relationship between the plDDT scores and the extent of RMSD?**

**R:** Yes. We are conducting further measurements for a potential follow-up study, but predictions that significantly deviate from the ground state reference have significantly lower average pLDDTs than those that mostly replicate the ground state. The amplitude of this change is reduced for some mutations, which we hypothesize is due to the alternative-conformation-stabilizing effect of these mutations.

**9. Prior work that uses AF2 to sample conformational ensembles has seen that AF2 is able to predict more diverse conformations when the protein is not part of AF2’s training dataset. Was GMSCF part of the training dataset? If yes, how would the author’s approach vary for a protein that is not part of the training dataset?**

**R:** Yes. GMCSF is a part of the training dataset as the crystal structures that were solved for it were deposited before 2018 (1CSG, 2GMF, 5C7X, among others). We are measuring/honing our approach in proteins that were not part of the training dataset and will publish the results in a public repository once this analysis is finished.

**10. Some of the figures are not informative/important enough to be main figures. For example, Figure 2 is mainly the AF2 pipeline, Figure 5 is just a pictorial representation of Supplementary Table S1. Also, Figures 6 and 7 could be combined into a single figure.**

**R:** We agree with this assessment and will remove redundancies in a future version of this manuscript.

**11. The CSP data for H15N is not shown in Figure 9B whereas its RMSD is shown in Figure 9C**

**R:** There is no CSP data for H15N at the time of publication of the pre-print or of this review response, because of challenges with expressing this mutant. We made the decision to include the prediction data for the sake of transparency and completeness. Accordingly, we don’t make any assessments regarding the precision of the H15N predictions and don’t count them when calculating accuracy rates. We will upload an amended version of Figure 9 when GMSCF H15N data is collected.

Figure 9B compares the RMSD of each mutant ensemble to the wild-type, ground state structure. 

**12. The cut-off values used for jackhmmer not mentioned.**

**R:** We used the default jackhmmr Significance E-values as cut-off values (0.01 for Sequence, 0.03 for hit).

**13. Residues are being addressed as codons in some places in the text**

**R:** We will fix these inconsistencies in a future version of this manuscript.\

**14. The authors may also want to include a few sentences contrasting their approach with this recently posted work:** [**https://www.biorxiv.org/con...**](https://disq.us/url?url=https%3A%2F%2Fwww.biorxiv.org%2Fcontent%2F10.1101%2F2023.08.06.552168v1%3APysS4uIBO-AeLcbCleCcvEGDfec\&cuid=2634513) **in the introduction or discussion.**

**R:** Agreed, especially given the similarity of the systems tested. We will discuss this study and how our works relate to it in a future version of this manuscript.

# **References:**

1\. Xie, T., Saleh, T., Rossi, P. & Kalodimos, C. G. Conformational states dynamically populated by a kinase determine its function. _Science_ **370** (Oct. 2020).

2\. Hoemberger, M., Pitsawong, W. & Kern, D. Cumulative mechanism of several major imatinib-resistant mutations in Abl kinase. _Proceedings of the National Academy of Sciences_ **117,** 19221–19227 (July 2020).

3\. Shan, Y. et al. A conserved protonation-dependent switch controls drug binding in the Abl kinase. _Proceedings of the National Academy of Sciences_ **106** 139–144 (2009).

4\. Lakkaniga, N. R., Balasubramaniam, M., Zhang, S., Frett, B. & Li, H. Structural Characterization of the Aurora Kinase B “DFG-flip” Using Metadynamics. The AAPS Journal **22** (2019).
