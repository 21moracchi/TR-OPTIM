function mpc = case5_critical_opf

%   MATPOWER

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
    1    3     81.00    37.00    0    0    1     1.0    0.0    345    1    1.05    0.95;
    2    2    143.00    75.00    0    0    1     1.0    0.0    345    1    1.05    0.95;
    3    2     94.00    51.00    0    0    1     1.0    0.0    345    1    1.05    0.95;
    4    2     87.00    37.00    0    0    1     1.0    0.0    345    1    1.05    0.95;
    5    2    106.00    50.00    0    0    1     1.0    0.0    345    1    1.05    0.95;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
    1    135     35.00      230.00    -170.00    1.0    100    1    270.00    0.0    0 0   0   0   0   0   0   0   0   0   0;
    2    0.0     0.0          0.00          0    1.0    100    1         0    0.0    0 0   0   0   0   0   0   0   0   0   0;
    3     90     48.65      150.00    -140.00    1.0    100    1    180.00    0.0    0 0   0   0   0   0   0   0   0   0   0;
    4     90     35.85      150.00    -140.00    1.0    100    1    180.00    0.0    0 0   0   0   0   0   0   0   0   0   0;
    5    0.0     0.0          0.00          0    1.0    100    1         0    0.0    0 0   0   0   0   0   0   0   0   0   0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
    1    2    0.0420    0.1680    0.0300    230.0     0.0     0.0     0.0     0.0    1    -360    360;
    2    3    0.0310    0.1260    0.0200    230.0     0.0     0.0     0.0     0.0    1    -360    360;
    3    5    0.0530    0.2100    0.0150    230.0     0.0     0.0     0.0     0.0    1    -360    360;
    3    4    0.0840    0.3360    0.0120    230.0     0.0     0.0     0.0     0.0    1    -360    360;
    5    4    0.0630    0.2520    0.0110    230.0     0.0     0.0     0.0     0.0    1    -360    360;
    5    1    0.0310    0.1260    0.0100    230.0     0.0     0.0     0.0     0.0    1    -360    360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
    2   1500    0   3   10.0    200.0  0.0 ;
    2   2000    0   3    0.0      0.0  0.0 ;
    2   3000    0   3   20.0    400.0  0.0 ;
    2   3000    0   3   30.0    100.0  0.0 ;
    2   3000    0   3    0.0      0.0  0.0 ;
];
