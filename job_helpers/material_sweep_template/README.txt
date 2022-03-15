Process
1) generate input parameters in matlab
2) In each material directory (/xxxYield_xxxPinf/), edit the local input_baseline.dat to have the proper inputs for the desired shock:
	&PROBINPUT
	p_infty =
	gamma   =
	u_L     =
	u_R     =
	v_L     =
	v_R     =
	rho_L   =
	rho_R   =
	p_L     =
	p_R     =
	ge11_L  =
	ge11_R  =
	ge21_L  =
	ge21_R  =
3) Edit functions.sh reflect the desired sweeps in LAD coefficients for all materials.
4) Run make_sweeps.sh to create all the input files to sweep over LAD coefficients for all materials.
5) Run run.sh to submit the jobs for all materials.
