function [baseMVA, bus, gen, branch, areas, gencost] = WB2
%CASE4GS  Power flow data for 2 bus, 1 gen case from Tate & Overbye
%   Please see 'help caseformat' for details on the case file format.
%
%   This is the 2 bus example from pp. 1667-1674 of "IEEE TRANSACTIONS ON POWER SYSTEMS",
%   A Comparison of the Optimal Multiplier in Polar and Rectangular Coordinates, 2005.

%   MATPOWER
%   $Id: case2gs.m,v 1.1 2008/09/16 18:00:00 $

%%-----  Power Flow Data  -----%%
%% system MVA base
baseMVA = 100;

%% bus data
%	bus_i	type	Pd	    Qd   Gs	Bs	area	Vm	    Va	  baseKV	zone	Vmax	Vmin
bus = [
	1	    3	    0.00	    0.00    0	0	1	    1.0    	0	    230	    1	    1.05	.95;
	2	    1	    3.5         -3.5   0	0	1	    1.0     0       230	    1	    1.05	.95;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	    Vg	mBase	status	Pmax	Pmin
gen = [
	1   0	0	10   	-10	        1 	100  	1	    10	    -10;
];

%% branch data
%	fbus	tbus	r	 x	    b	    rateA	rateB	rateC	ratio	angle	status
branch = [
	1	    2	   .04e2	 .2e2	0.000       250	    250	    250	    0	    0	    1;
];

%% area data
areas = [
	1	1;
];

%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0

gencost = [
	2	0	0	3	0	2	0;
];

return;
