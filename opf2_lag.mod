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
# Shared
###
param IS_PHASE_ONE default 0;
param dual_convexity       default 0;
param dual_re{bus in bus_i}default 0;
param dual_im{bus in bus_i}default 0;

param nIte default 0;


param P_gen{ite in 1..nIte, gen in gen_i}default 0;
param Q_gen{ite in 1..nIte, gen in gen_i}default 0;

param V_re{ite in 1..nIte, bus in bus_i}default 0;
param V_im{ite in 1..nIte, bus in bus_i}default 0;

param I_re{ite in 1..nIte, bus in bus_i}default 0;
param I_im{ite in 1..nIte, bus in bus_i}default 0;

param I_or_re{ite in 1..nIte, (l,m) in branch} := 100*(+y11_re[l,m]*V_re[ite, l]-y11_im[l,m]*V_im[ite, l]+y12_re[l,m]*V_re[ite, m]-y12_im[l,m]*V_im[ite, m]);
param I_or_im{ite in 1..nIte, (l,m) in branch} := 100*(+y11_re[l,m]*V_im[ite, l]+y11_im[l,m]*V_re[ite, l]+y12_re[l,m]*V_im[ite, m]+y12_im[l,m]*V_re[ite, m]);

param I_ex_re{ite in 1..nIte, (l,m) in branch} := 100*(+y21_re[l,m]*V_re[ite, l]-y21_im[l,m]*V_im[ite, l]+y22_re[l,m]*V_re[ite, m]-y22_im[l,m]*V_im[ite, m]);
param I_ex_im{ite in 1..nIte, (l,m) in branch} := 100*(+y21_re[l,m]*V_im[ite, l]+y21_im[l,m]*V_re[ite, l]+y22_re[l,m]*V_im[ite, m]+y22_im[l,m]*V_re[ite, m]);

param Cost{ite in 1..nIte} := sum{gen in gen_i}(c0[gen]+c1[gen]*P_gen[ite, gen]+c2[gen]*P_gen[ite, gen]^2);
###
# Master
###
problem master;
var fake_lambda >= 0;
var lambda{ite in 1..nIte} >=0;

subject to i_balance_re{bus in bus_i}:
  +sum{ite in 1..nIte}lambda[ite]*I_re[ite, bus]
  +sum{ite in 1..nIte, (bus1, bus2) in branch:bus==bus1}lambda[ite]*I_or_re[ite, bus1, bus2]
  +sum{ite in 1..nIte, (bus1, bus2) in branch:bus==bus2}lambda[ite]*I_ex_re[ite, bus1, bus2]
  =
  0
  ;
subject to i_balance_im{bus in bus_i}:
  +sum{ite in 1..nIte}lambda[ite]*I_im[ite, bus]
  +sum{ite in 1..nIte, (bus1, bus2) in branch:bus==bus1}lambda[ite]*I_or_im[ite, bus1, bus2]
  +sum{ite in 1..nIte, (bus1, bus2) in branch:bus==bus2}lambda[ite]*I_ex_im[ite, bus1, bus2]
  =
  0
  ;
subject to convexity: (if IS_PHASE_ONE == 1 then fake_lambda)+sum{ite in 1..nIte}lambda[ite]=1;

minimize master_obj: 
	+sum{ite in 1..nIte:IS_PHASE_ONE==0}lambda[ite]*Cost[ite]
	+(if IS_PHASE_ONE==1 then fake_lambda)
	;


###
# Slave
###
problem slave;
var p_gen{gen in gen_i};
var q_gen{gen in gen_i};

var v_re{bus_i};
var v_im{bus_i};

var i_re{bus in bus_i};
var i_im{bus in bus_i};

subject to lb_p_gen{gen in gen_i}: p_gen[gen] >= Pmin[gen]*100;
subject to ub_p_gen{gen in gen_i}: p_gen[gen] <= Pmax[gen]*100;

subject to lb_q_gen{gen in gen_i}: q_gen[gen] >= Qmin[gen]*100;
subject to ub_q_gen{gen in gen_i}: q_gen[gen] <= Qmax[gen]*100;

subject to ctr_i_re{bus in bus_i}:
  +i_re[bus]*v_re[bus]
  +i_im[bus]*v_im[bus] 
  = 
  +(+Pd[bus]-(if bus in gen_i then p_gen[bus])/100)
  ;
subject to ctr_i_im{bus in bus_i}:  
  +i_re[bus]*v_im[bus]
  -i_im[bus]*v_re[bus]
  =
  +(+Qd[bus]-(if bus in gen_i then q_gen[bus])/100)
  ;
  
subject to ref_bus_re{bus in bus_i:type[bus]==3}:v_re[bus]>= 0;
subject to ref_bus_im{bus in bus_i:type[bus]==3}:v_im[bus] = 0;

subject to ctr_v_min{bus in bus_i}: Vmin[bus]^2 <= v_re[bus]*v_re[bus]+v_im[bus]*v_im[bus];
subject to ctr_v_max{bus in bus_i}: Vmax[bus]^2 >= v_re[bus]*v_re[bus]+v_im[bus]*v_im[bus];
#subject to ctr_v_min{bus in bus_i}: -Vmax[bus] <= v_re[bus] <= Vmax[bus];
#subject to ctr_v_max{bus in bus_i}: -Vmax[bus] <= v_im[bus] <= Vmax[bus];

minimize slave_obj: 
  +sum{gen in gen_i:IS_PHASE_ONE==0}(c0[gen]+c1[gen]*p_gen[gen]+c2[gen]*p_gen[gen]^2)
  -sum{bus in bus_i}(dual_re[bus]*(
    +i_re[bus]
    +sum{(l, m) in branch:bus==l}(100*(+y11_re[l,m]*v_re[l]-y11_im[l,m]*v_im[l]+y12_re[l,m]*v_re[m]-y12_im[l,m]*v_im[m]))
    +sum{(l, m) in branch:bus==m}(100*(+y21_re[l,m]*v_re[l]-y21_im[l,m]*v_im[l]+y22_re[l,m]*v_re[m]-y22_im[l,m]*v_im[m]))
  ))
  -sum{bus in bus_i}(dual_im[bus]*(
    +i_im[bus]
    +sum{(l, m) in branch:bus==l}(100*(+y11_re[l,m]*v_im[l]+y11_im[l,m]*v_re[l]+y12_re[l,m]*v_im[m]+y12_im[l,m]*v_re[m]))
    +sum{(l, m) in branch:bus==m}(100*(+y21_re[l,m]*v_im[l]+y21_im[l,m]*v_re[l]+y22_re[l,m]*v_im[m]+y22_im[l,m]*v_re[m]))  	
  ))
  -dual_convexity
;

# PQ
#subject to ref_bus_re{bus in bus_i:type[bus]==1}:v_re[bus]>=0;
#subject to ref_bus_im{bus in bus_i:type[bus]==1}:v_im[bus]= 0;
# PV
#subject to ref_bus_re{bus in bus_i:type[bus]==2}:v_re[bus]=Vm[bus]*cos(Va[bus]);
#subject to ref_bus_im{bus in bus_i:type[bus]==2}:v_im[bus]=Vm[bus]*sin(Va[bus]);

