# CFM-Phyto_Temp
Cell Flux Model with the added stress of temperature. Different light and nutrient (3) regimes available, along with variation in growth rate. 

The three different nutrient regimes and the corresponding python file are listed below:

1. Nutrient co-limitation: N_P-colimiting_scenario_FINAL-Github

2. Phosphorus limitation: N_P-Plimiting-Github

3. Nitrogen Limitation: N_P-Nlimiting-Github

To change the growth rate, change the parameter on line 54 (for all 3 nutrient regimes) labeled "Dd". This unit is measured in dilution per day.

To change the light intensity, change the parameter on line 32 labeled "I". This is irradiance measured in umol of photons per meter^2 per second.


The resulting code will produce 4 figures. 

Figure 1: Carbon allocation with increasing temperature to 4 macromolecular pools: essential molecules ("Other"), biosynthesis, photosynthesis, and carbon storage. 

Figure 2: N:C ratio with increasing temperature. Color-coded by macromolecular allocation to photosynthesis, biosynthesis, essential molecules, and, in the case of P-limiting, N-storage. 

Figure 3: P:C ratio with increasing temperature. Color-coded by macromolecular allocation to photosynthesis, biosynthesis, essential molecules, and, in the case of N-limiting, P-storage. 

Figure 4: Model comparison to culture data. N:P ratio with increasing temperature

