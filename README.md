# CFM-Phyto_Temp
Cell Flux Model with the added stress of temperature. Different light and nutrient (3) regimes available, along with variation in growth rate. 

The three different nutrient regimes and the corresponding python file are listed below:

1. Nutrient co-limitation: N_P-colimiting_scenario_FINAL-Github

2. Phosphorus limitation: N_P-Plimiting-Github

3. Nitrogen Limitation: N_P-Nlimiting-Github

There are also two light regimes provided as seen in the supplementary material:

1. N_P-colimiting_scenario_low-light

2. N_P-colimiting_scenario_high-light

There are two different growth rates are provided as seen in the supplementary material:

1. N_P-colimiting_scenario_0.15

2. N_P-colimiting_scenario_0.35

If you prefer to use the temperature formulation of Q10 rather than the Arrhenius form, a simulation in the co-limiting scenario is available:

1. N_P-colimiting_scenario_FINAL-Q10


You can change the set values of light intensity or growth rate by:

To change the growth rate, change the parameter on line 54 (for all 3 nutrient regimes) labeled "Dd". This unit is measured in dilution per day.

To change the light intensity, change the parameter on line 32 labeled "I". This is irradiance measured in umol of photons per meter^2 per second.


The resulting code will produce 4 figures. 

Figure 1: Carbon allocation with increasing temperature to 4 macromolecular pools: essential molecules ("Other"), biosynthesis, photosynthesis, and carbon storage. 

Figure 2: N:C ratio with increasing temperature. Color-coded by macromolecular allocation to photosynthesis, biosynthesis, essential molecules, and, in the case of P-limiting, N-storage. 

Figure 3: P:C ratio with increasing temperature. Color-coded by macromolecular allocation to photosynthesis, biosynthesis, essential molecules, and, in the case of N-limiting, P-storage. 

Figure 4: Model comparison to culture data. N:P ratio with increasing temperature

Lastly, our model results are reflected in the global surface ocean using the following:

  Data input:
  
  1. X
  2. Y
  3. Z
  
  Codes:
  
  1. final_verison2

  2. t002_reading_plot

The resulting code produces 9 figures total. 

Figure 1: N:P values on a global plot for current ocean temperatures

Figure 2: P:C values on a global plot for current ocean temperatures

Figure 3: N:C values on a global plot for current ocean temperatures

Figure 4: N:P values on a global plot for future ocean temperatures (+4C)

Figure 5: P:C values on a global plot for future ocean temperatures (+4C)

Figure 6: N:C values on a global plot for future ocean temperatures (+4C)

Figure 7: Change in N:P values between these two scenarios on a global plot

Figure 8: Change in P:C values between these two scenarios on a global plot

Figure 9: Change in N:C values between these two scenarios on a global plot
