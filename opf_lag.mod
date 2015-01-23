reset;
param baseMVA;
set bus_i;
param type  {bus_i};
param Pd    {bus_i};
param Qd    {bus_i};
param Gs    {bus_i};
param Bs    {bus_i};
param area  {bus_i};
param Vm    {bus_i};
param Va    {bus_i};
param baseKV{bus_i};
param zone  {bus_i};
param Vmax  {bus_i};
param Vmin  {bus_i};

set gen_i within bus_i;
param Pg        {gen_i};
param Qg        {gen_i};
param Qmax      {gen_i};
param Qmin      {gen_i};
param Vg        {gen_i};
param mBase     {gen_i};
param gen_status{gen_i};
param Pmax      {gen_i};
param Pmin      {gen_i};



set branch dimen 2;
param r            {branch};
param x            {branch};
param b            {branch};
param rateA        {branch};
param rateB        {branch};
param rateC        {branch};
param ratio        {branch};
param angle        {branch};
param branch_status{branch};


param startup {gen_i};
param shutdown{gen_i};
param n       {gen_i};
param c0      {gen_i};
param c1      {gen_i};
param c2      {gen_i};

###
# other parameter
###
param line_G{(l,m) in branch}   :=  r[l,m]/(r[l,m]^2+x[l,m]^2);
param line_B{(l,m) in branch}   := -x[l,m]/(r[l,m]^2+x[l,m]^2);

param pst_ratio{(l,m) in branch} := 1;
param pst_angle{(l,m) in branch} := 0;
param pst_cos{(l,m) in branch} := cos(pst_angle[l,m]);
param pst_sin{(l,m) in branch} := sin(pst_angle[l,m]);

param y11_re{(l,m) in branch}:=  line_G[l,m]/(pst_ratio[l,m]^2);
param y11_im{(l,m) in branch}:=  (line_B[l,m]+b[l,m]/2)/(pst_ratio[l,m]^2);

param y22_re{(l,m) in branch}:=  line_G[l,m];
param y22_im{(l,m) in branch}:=  line_B[l,m]+b[l,m]/2;
# -ys/(rho exp +jTheta)
param y12_re{(l,m) in branch}:=  (-1/pst_ratio[l,m])*( line_G[l,m]*pst_cos[l,m]-line_B[l,m]*pst_sin[l,m]);
param y12_im{(l,m) in branch}:=  (-1/pst_ratio[l,m])*( line_B[l,m]*pst_cos[l,m]+line_G[l,m]*pst_sin[l,m]);
# -ys/(rho exp -jTheta)
param y21_re{(l,m) in branch}:=  (-1/pst_ratio[l,m])*( line_G[l,m]*pst_cos[l,m]+line_B[l,m]*pst_sin[l,m]);
param y21_im{(l,m) in branch}:=  (-1/pst_ratio[l,m])*( line_B[l,m]*pst_cos[l,m]-line_G[l,m]*pst_sin[l,m]);
###


###
# OPF formulation
###
  
param dual_p{bus in bus_i} default 0;
param dual_q{bus in bus_i} default 0;

#param dual_vmin{bus in bus_i} default 0;
#param dual_vmax{bus in bus_i} default 0;
param dual_convexity default 0;  
 
param nIte default 0;
 
param V_re{ite in 1..nIte, bus in bus_i};
param V_im{ite in 1..nIte, bus in bus_i};
 
param  P_or{ite in 1..nIte, (bus1,bus2) in branch} = 100*(
	+y11_re[bus1, bus2]*V_re[ite,bus1]*V_re[ite,bus1]
	+y11_re[bus1, bus2]*V_im[ite,bus1]*V_im[ite,bus1]
	+y12_re[bus1, bus2]*V_re[ite,bus1]*V_re[ite,bus2]
	+y12_re[bus1, bus2]*V_im[ite,bus1]*V_im[ite,bus2]
	+y12_im[bus1, bus2]*V_im[ite,bus1]*V_re[ite,bus2]
	-y12_im[bus1, bus2]*V_re[ite,bus1]*V_im[ite,bus2]);
	
param  P_ex{ite in 1..nIte, (bus1,bus2) in branch} = 100*(
	+y22_re[bus1,bus2]*V_re[ite,bus2]*V_re[ite,bus2]
	+y22_re[bus1,bus2]*V_im[ite,bus2]*V_im[ite,bus2]
	+y21_re[bus1,bus2]*V_re[ite,bus1]*V_re[ite,bus2]
	+y21_re[bus1,bus2]*V_im[ite,bus1]*V_im[ite,bus2]
	+y21_im[bus1,bus2]*V_re[ite,bus1]*V_im[ite,bus2]
	-y21_im[bus1,bus2]*V_im[ite,bus1]*V_re[ite,bus2]);

param  Q_or{ite in 1..nIte, (bus1,bus2) in branch} = 100*(
	  -y11_im[bus1,bus2]*V_re[ite,bus1]*V_re[ite,bus1]
	  -y11_im[bus1,bus2]*V_im[ite,bus1]*V_im[ite,bus1]
	  +y12_re[bus1,bus2]*V_im[ite,bus1]*V_re[ite,bus2]
	  -y12_re[bus1,bus2]*V_re[ite,bus1]*V_im[ite,bus2]
	  -y12_im[bus1,bus2]*V_re[ite,bus1]*V_re[ite,bus2]
	  -y12_im[bus1,bus2]*V_im[ite,bus1]*V_im[ite,bus2]);
	
param  Q_ex{ite in 1..nIte, (bus1,bus2) in branch} = 100*(
	-y22_im[bus1,bus2]*V_re[ite,bus2]*V_re[ite,bus2]
	  -y22_im[bus1,bus2]*V_im[ite,bus2]*V_im[ite,bus2]
	  +y21_re[bus1,bus2]*V_re[ite,bus1]*V_im[ite,bus2]
	  -y21_re[bus1,bus2]*V_im[ite,bus1]*V_re[ite,bus2]
	  -y21_im[bus1,bus2]*V_re[ite,bus1]*V_re[ite,bus2]
	  -y21_im[bus1,bus2]*V_im[ite,bus1]*V_im[ite,bus2]);

set master_basis default {};
###
# master
###
problem master;
param FEASIBILITY default 1;
var p_gen{gen in gen_i}
  >= Pmin[gen]*100, <= Pmax[gen]*100
  ;
var q_gen{gen in gen_i} 
  >= Qmin[gen]*100, <= Pmax[gen]*100
  ;
var lambda{ite in 1..nIte} >= 0;
var fake_lambda >= 0;

var p_pos{bus in bus_i} >= 0;
var p_neg{bus in bus_i} >= 0;

var q_pos{bus in bus_i} >= 0;
var q_neg{bus in bus_i} >= 0;

#var v_min_slack{bus in bus_i} >= 0;
#var v_max_slack{bus in bus_i} >= 0;
subject to convexity: 
	+(if FEASIBILITY == 1 then fake_lambda)+sum{ite in 1..nIte} lambda[ite]=1;

subject to p_balance{bus in bus_i}:
	+(if FEASIBILITY == 1 then +p_pos[bus]-p_neg[bus])
 	-(if bus in gen_i then p_gen[bus])/100
  	+sum{ite in 1..nIte, (bus1, bus2) in branch:bus==bus1} P_or[ite,bus1,bus2]*lambda[ite]
  	+sum{ite in 1..nIte, (bus1, bus2) in branch:bus==bus2} P_ex[ite,bus1,bus2]*lambda[ite]
  	=
  	-Pd[bus]
  ;
subject to q_balance{bus in bus_i}:
	+(if FEASIBILITY == 1 then +q_pos[bus]-q_neg[bus])
  	-(if bus in gen_i then q_gen[bus])/100
  	+sum{ite in 1..nIte, (bus1, bus2) in branch:bus==bus1} Q_or[ite,bus1,bus2]*lambda[ite]
  	+sum{ite in 1..nIte, (bus1, bus2) in branch:bus==bus2} Q_ex[ite,bus1,bus2]*lambda[ite]
  	=
  	-Qd[bus]
  ;

minimize master_obj:
	+if FEASIBILITY == 1 then(
		+10e5*fake_lambda
		+10e5*sum{bus in bus_i} (
			+p_pos[bus]+p_neg[bus]
			+q_pos[bus]+q_neg[bus]
			#+v_min_slack[bus]+v_max_slack[bus]
		)
	)else
		+sum{gen in gen_i}(c0[gen]+c1[gen]*p_gen[gen]+c2[gen]*p_gen[gen]^2)
	;
#subject to ctr_v_min{bus in bus_i}: sum{ite in 1..nIte}(V_re[ite,bus]^2+V_im[ite,bus]^2)*lambda[ite] + v_min_slack[bus] >= Vmin[bus]^2 ;
#subject to ctr_v_max{bus in bus_i}: sum{ite in 1..nIte}(V_re[ite,bus]^2+V_im[ite,bus]^2)*lambda[ite] - v_max_slack[bus] <= Vmax[bus]^2;
###
# sub
###
problem subproblem: p_gen, q_gen, lambda, p_pos, p_neg, q_pos, q_neg;

var v_re{bus in bus_i}, >= -Vmax[bus], <= Vmax[bus];
var v_im{bus in bus_i}, >= -Vmax[bus], <= Vmax[bus];

# convex
var ips_lambda >= 0;
subject to ref_bus_re{bus in bus_i:type[bus]==3}:v_re[bus]>= 0;
subject to ref_bus_im{bus in bus_i:type[bus]==3}:v_im[bus] = 0;

subject to ctr_v_min{bus in bus_i}: v_re[bus]*v_re[bus]+v_im[bus]*v_im[bus] >= Vmin[bus]^2;
subject to ctr_v_max{bus in bus_i}: v_re[bus]*v_re[bus]+v_im[bus]*v_im[bus] <= Vmax[bus]^2;

minimize sub_obj:
-sum{bus in bus_i, (bus1, bus2) in branch:bus==bus1} dual_p[bus]*100*(
		+y11_re[bus1, bus2]*v_re[bus1]*v_re[bus1]
		+y11_re[bus1, bus2]*v_im[bus1]*v_im[bus1]
		+y12_re[bus1, bus2]*v_re[bus1]*v_re[bus2]
		+y12_re[bus1, bus2]*v_im[bus1]*v_im[bus2]
		+y12_im[bus1, bus2]*v_im[bus1]*v_re[bus2]
		-y12_im[bus1, bus2]*v_re[bus1]*v_im[bus2])
	-sum{bus in bus_i, (bus1, bus2) in branch:bus==bus2} dual_p[bus]*100*(
		+y22_re[bus1,bus2]*v_re[bus2]*v_re[bus2]
		+y22_re[bus1,bus2]*v_im[bus2]*v_im[bus2]
		+y21_re[bus1,bus2]*v_re[bus1]*v_re[bus2]
		+y21_re[bus1,bus2]*v_im[bus1]*v_im[bus2]
		+y21_im[bus1,bus2]*v_re[bus1]*v_im[bus2]
		-y21_im[bus1,bus2]*v_im[bus1]*v_re[bus2])
	-sum{bus in bus_i, (bus1, bus2) in branch:bus==bus1} dual_q[bus]*100*(
 		-y11_im[bus1,bus2]*v_re[bus1]*v_re[bus1]
		 -y11_im[bus1,bus2]*v_im[bus1]*v_im[bus1]
		 +y12_re[bus1,bus2]*v_im[bus1]*v_re[bus2]
		 -y12_re[bus1,bus2]*v_re[bus1]*v_im[bus2]
		 -y12_im[bus1,bus2]*v_re[bus1]*v_re[bus2]
		 -y12_im[bus1,bus2]*v_im[bus1]*v_im[bus2])
	-sum{bus in bus_i, (bus1, bus2) in branch:bus==bus2} dual_q[bus]*100*(
		-y22_im[bus1,bus2]*v_re[bus2]*v_re[bus2]
  		-y22_im[bus1,bus2]*v_im[bus2]*v_im[bus2]
  		+y21_re[bus1,bus2]*v_re[bus1]*v_im[bus2]
  		-y21_re[bus1,bus2]*v_im[bus1]*v_re[bus2]
  		-y21_im[bus1,bus2]*v_re[bus1]*v_re[bus2]
  		-y21_im[bus1,bus2]*v_im[bus1]*v_im[bus2])
	-dual_convexity;

subject to ips_convexity{tmp in 1..1:FEASIBILITY==0 and card(master_basis)>0}: 
	ips_lambda+sum{ite in master_basis} lambda[ite]=1;
	
subject to ips_p_balance{bus in bus_i: FEASIBILITY==0 and card(master_basis)>0}:
 	-(if bus in gen_i then p_gen[bus])/100
  	+sum{ite in master_basis, (bus1, bus2) in branch:bus==bus1} P_or[ite,bus1,bus2]*lambda[ite]
  	+sum{ite in master_basis, (bus1, bus2) in branch:bus==bus2} P_ex[ite,bus1,bus2]*lambda[ite]
	+sum{(bus1, bus2) in branch:bus==bus1}ips_lambda*100*(
		+y11_re[bus1, bus2]*v_re[bus1]*v_re[bus1]
		+y11_re[bus1, bus2]*v_im[bus1]*v_im[bus1]
		+y12_re[bus1, bus2]*v_re[bus1]*v_re[bus2]
		+y12_re[bus1, bus2]*v_im[bus1]*v_im[bus2]
		+y12_im[bus1, bus2]*v_im[bus1]*v_re[bus2]
		-y12_im[bus1, bus2]*v_re[bus1]*v_im[bus2])
	+sum{(bus1, bus2) in branch:bus==bus2} ips_lambda*100*(
		+y22_re[bus1,bus2]*v_re[bus2]*v_re[bus2]
		+y22_re[bus1,bus2]*v_im[bus2]*v_im[bus2]
		+y21_re[bus1,bus2]*v_re[bus1]*v_re[bus2]
		+y21_re[bus1,bus2]*v_im[bus1]*v_im[bus2]
		+y21_im[bus1,bus2]*v_re[bus1]*v_im[bus2]
		-y21_im[bus1,bus2]*v_im[bus1]*v_re[bus2])
  	=
  	-Pd[bus]
  ;
subject to ips_q_balance{bus in bus_i: FEASIBILITY==0 and card(master_basis)>0}:
  	-(if bus in gen_i then q_gen[bus])/100
  	+sum{ite in master_basis, (bus1, bus2) in branch:bus==bus1} Q_or[ite,bus1,bus2]*lambda[ite]
  	+sum{ite in master_basis, (bus1, bus2) in branch:bus==bus2} Q_ex[ite,bus1,bus2]*lambda[ite]
	+sum{(bus1, bus2) in branch:bus==bus1} ips_lambda*100*(
 		-y11_im[bus1,bus2]*v_re[bus1]*v_re[bus1]
		 -y11_im[bus1,bus2]*v_im[bus1]*v_im[bus1]
		 +y12_re[bus1,bus2]*v_im[bus1]*v_re[bus2]
		 -y12_re[bus1,bus2]*v_re[bus1]*v_im[bus2]
		 -y12_im[bus1,bus2]*v_re[bus1]*v_re[bus2]
		 -y12_im[bus1,bus2]*v_im[bus1]*v_im[bus2])
	+sum{(bus1, bus2) in branch:bus==bus2} ips_lambda*100*(
		-y22_im[bus1,bus2]*v_re[bus2]*v_re[bus2]
  		-y22_im[bus1,bus2]*v_im[bus2]*v_im[bus2]
  		+y21_re[bus1,bus2]*v_re[bus1]*v_im[bus2]
  		-y21_re[bus1,bus2]*v_im[bus1]*v_re[bus2]
  		-y21_im[bus1,bus2]*v_re[bus1]*v_re[bus2]
  		-y21_im[bus1,bus2]*v_im[bus1]*v_im[bus2])
  	=
  	-Qd[bus]
  ;

