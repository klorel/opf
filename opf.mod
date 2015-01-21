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
  >= Pmin[gen]/100, <= Pmax[gen]/100
  ;
var q_gen{gen in gen_i} 
  >= Qmin[gen]/100, <= Pmax[gen]/100
  ;

# linear equation I=YV
# power flow equation
# s_or = v_or i_or* = v_or(y_11* v_or* + y_12* v_ex*) = y_11* v_or v_or* + y_12* v_or v_ex*
# s_ex = v_ex i_ex* = v_ex(y_21* v_or* + y_22* v_ex*) = y_21* v_ex v_or* + y_22* v_ex v_ex*
# s_or = y_11*(v_or_re^2 +v_or_im^2) + y_12* v_or v_ex*)
# s_ex = y_21* v_ex v_or* + y_22* (v_ex_re^2 +v_ex_im^2)

var v_hProd_re{i in bus_i, j in bus_i} = v_re[i]*v_re[j]+v_im[i]*v_im[j];

var v_hProd_im{i in bus_i, j in bus_i} = 
  if(i == j) then(
    0
  )else(
    v_re[j]*v_im[i]-v_re[i]*v_im[j]
  );
# i_or = y_11 v_or + y_12 v_ex
# i_ex = y_21 v_or + y_22 v_ex

var  i_or_re{(l,m) in branch} = +y11_re[l,m]*v_re[l]-y11_im[l,m]*v_im[l]+y12_re[l,m]*v_re[m]-y12_im[l,m]*v_im[m];
var  i_or_im{(l,m) in branch} = +y11_re[l,m]*v_im[l]+y11_im[l,m]*v_re[l]+y12_re[l,m]*v_im[m]+y12_im[l,m]*v_re[m];
var  i_ex_re{(l,m) in branch} = +y21_re[l,m]*v_re[l]-y21_im[l,m]*v_im[l]+y22_re[l,m]*v_re[m]-y22_im[l,m]*v_im[m];
var  i_ex_im{(l,m) in branch} = +y21_re[l,m]*v_im[l]+y21_im[l,m]*v_re[l]+y22_re[l,m]*v_im[m]+y22_im[l,m]*v_re[m];


# s_or = y_11* v_or v_or* + y_12* v_or v_ex*
# s_ex = y_22* v_ex v_ex* + y_21* v_ex v_or*


var  p_or{(l,m) in branch} = +y11_re[l,m]*v_hProd_re[l,l]+y12_re[l,m]*v_hProd_re[l,m]+y12_im[l,m]*v_hProd_im[l,m];
var  q_or{(l,m) in branch} = -y11_im[l,m]*v_hProd_re[l,l]+y12_re[l,m]*v_hProd_im[l,m]-y12_im[l,m]*v_hProd_re[l,m];
var  p_ex{(l,m) in branch} = +y21_re[l,m]*v_hProd_re[m,l]+y21_im[l,m]*v_hProd_im[m,l]+y22_re[l,m]*v_hProd_re[m,m];
var  q_ex{(l,m) in branch} = +y21_re[l,m]*v_hProd_im[m,l]-y21_im[l,m]*v_hProd_re[m,l]-y22_im[l,m]*v_hProd_re[m,m];

subject to p_balance{bus in bus_i}:
  +Pd[bus]/100
  +sum{(bus, m) in branch} (+y11_re[bus,m]*(v_re[bus]*v_re[bus]+v_im[bus]*v_im[bus])+y12_re[bus,m]*(v_re[bus]*v_re[m]+v_im[bus]*v_im[m])+y12_im[bus,m]*(v_im[bus]*v_re[m]-v_re[bus]*v_im[m]))
  +sum{(l, bus) in branch} (+y21_re[l,bus]*(v_re[l]*v_re[bus]+v_im[l]*v_im[bus])+y21_im[l,bus]*(v_re[l]*v_im[bus]-v_im[l]*v_re[bus])+y22_re[l,bus]*(v_re[bus]*v_re[bus]+v_im[bus]*v_im[bus]))
  =
  (if bus in gen_i then p_gen[bus])
  ;
subject to q_balance{bus in bus_i}:
  +Qd[bus]/100
  +sum{(bus, m) in branch} (-y11_im[bus,m]*(v_re[bus]*v_re[bus]+v_im[bus]*v_im[bus])+y12_re[bus,m]*(v_im[bus]*v_re[m]-v_re[bus]*v_im[m])-y12_im[bus,m]*(v_re[bus]*v_re[m]+v_im[bus]*v_im[m]))
  +sum{(l, bus) in branch} (+y21_re[l,bus]*(v_re[l]*v_im[bus]-v_im[l]*v_re[bus])-y21_im[l,bus]*(v_re[l]*v_re[bus]+v_im[l]*v_im[bus])-y22_im[l,bus]*(v_re[bus]*v_re[bus]+v_im[bus]*v_im[bus]))
  =
  (if bus in gen_i then q_gen[bus])
  ;


minimize act_cost: sum{gen in gen_i}(c0[gen]+c1[gen]*p_gen[gen]+c2[gen]*p_gen[gen]^2);

subject to ctr_v_min{bus in bus_i}: Vmin[bus]^2 <= v_re[bus]^2+v_im[bus]^2;
subject to ctr_v_max{bus in bus_i}: Vmax[bus]^2 >= v_re[bus]^2+v_im[bus]^2;