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
var v_re{bus_i};
var v_im{bus_i};

var p_gen{gen in gen_i}
  >= Pmin[gen]*100, <= Pmax[gen]*100
  ;
var q_gen{gen in gen_i} 
  >= Qmin[gen]*100, <= Qmax[gen]*100
  ;

# linear equation I=YV
# power flow equation
# s_or = v_or i_or* = v_or(y_11* v_or* + y_12* v_ex*) = y_11* v_or v_or* + y_12* v_or v_ex*
# s_ex = v_ex i_ex* = v_ex(y_21* v_or* + y_22* v_ex*) = y_21* v_ex v_or* + y_22* v_ex v_ex*
# s_or = y_11*(v_or_re^2 +v_or_im^2) + y_12* v_or v_ex*)
# s_ex = y_21* v_ex v_or* + y_22* (v_ex_re^2 +v_ex_im^2)
var i_or_re{(l,m) in branch};
var i_or_im{(l,m) in branch};

var i_ex_re{(l,m) in branch};
var i_ex_im{(l,m) in branch};

# i_or = y_11 v_or + y_12 v_ex
# i_ex = y_21 v_or + y_22 v_ex
subject to  ctr_i_or_re{(l,m) in branch}: i_or_re[l,m] = +y11_re[l,m]*v_re[l]-y11_im[l,m]*v_im[l]+y12_re[l,m]*v_re[m]-y12_im[l,m]*v_im[m];
subject to  ctr_i_or_im{(l,m) in branch}: i_or_im[l,m] = +y11_re[l,m]*v_im[l]+y11_im[l,m]*v_re[l]+y12_re[l,m]*v_im[m]+y12_im[l,m]*v_re[m];

subject to  ctr_i_ex_re{(l,m) in branch}: i_ex_re[l,m] = +y21_re[l,m]*v_re[l]-y21_im[l,m]*v_im[l]+y22_re[l,m]*v_re[m]-y22_im[l,m]*v_im[m];
subject to  ctr_i_ex_im{(l,m) in branch}: i_ex_im[l,m] = +y21_re[l,m]*v_im[l]+y21_im[l,m]*v_re[l]+y22_re[l,m]*v_im[m]+y22_im[l,m]*v_re[m];

var i_re{bus in bus_i};
var i_im{bus in bus_i};


# s_or = v_or i_or*
# s_ex = v_ex i_ex*

subject to ctr_i_re{bus in bus_i}:
  +i_re[bus]*v_re[bus]
  +i_im[bus]*v_im[bus] 
  = 
  +(+Pd[bus]-(if bus in gen_i then p_gen[bus])/100)
  ;
subject to ctr_i_im{bus in bus_i}: 
  +i_im[bus]*v_re[bus] 
  -i_re[bus]*v_im[bus]
  =
  -(+Qd[bus]-(if bus in gen_i then q_gen[bus])/100)
  ;
# (Sd*-S_gen*)/v_or* = (Sd*-S_gen*)/|v_or|² v_or
subject to i_balance_re{bus in bus_i}:
  +i_re[bus]
  +sum{(bus1, bus2) in branch:bus==bus1}100*i_or_re[bus1, bus2]
  +sum{(bus1, bus2) in branch:bus==bus2}100*i_ex_re[bus1, bus2]
  =
  0
  ;
  
  
subject to i_balance_im{bus in bus_i}:
  +i_im[bus]
  +sum{(bus1, bus2) in branch:bus==bus1}100*i_or_im[bus1, bus2]
  +sum{(bus1, bus2) in branch:bus==bus2}100*i_ex_im[bus1, bus2]
  =
  0
  ;


minimize act_cost: sum{gen in gen_i}(c0[gen]+c1[gen]*p_gen[gen]+c2[gen]*p_gen[gen]^2);

subject to ctr_v_min{bus in bus_i}: Vmin[bus]^2 <= v_re[bus]*v_re[bus]+v_im[bus]*v_im[bus];
subject to ctr_v_max{bus in bus_i}: Vmax[bus]^2 >= v_re[bus]*v_re[bus]+v_im[bus]*v_im[bus];

subject to ref_bus_re{bus in bus_i:type[bus]==3}:v_re[bus]>= 0;
subject to ref_bus_im{bus in bus_i:type[bus]==3}:v_im[bus] = 0;

# PQ
#subject to ref_bus_re{bus in bus_i:type[bus]==1}:v_re[bus]>=0;
#subject to ref_bus_im{bus in bus_i:type[bus]==1}:v_im[bus]= 0;
# PV
#subject to ref_bus_re{bus in bus_i:type[bus]==2}:v_re[bus]=Vm[bus]*cos(Va[bus]);
#subject to ref_bus_im{bus in bus_i:type[bus]==2}:v_im[bus]=Vm[bus]*sin(Va[bus]);

