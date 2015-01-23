
var  p_or{(bus1,bus2) in branch} = 
	+y11_re[bus1, bus2]*v_re[bus1]*v_re[bus1]
	+y11_re[bus1, bus2]*v_im[bus1]*v_im[bus1]
	+y12_re[bus1, bus2]*v_re[bus1]*v_re[bus2]
	+y12_re[bus1, bus2]*v_im[bus1]*v_im[bus2]
	+y12_im[bus1, bus2]*v_im[bus1]*v_re[bus2]
	-y12_im[bus1, bus2]*v_re[bus1]*v_im[bus2];
	
var  p_ex{(bus1,bus2) in branch} = 
	+y22_re[bus1,bus2]*v_re[bus2]*v_re[bus2]
	+y22_re[bus1,bus2]*v_im[bus2]*v_im[bus2]
	+y21_re[bus1,bus2]*v_re[bus1]*v_re[bus2]
	+y21_re[bus1,bus2]*v_im[bus1]*v_im[bus2]
	+y21_im[bus1,bus2]*v_re[bus1]*v_im[bus2]
	-y21_im[bus1,bus2]*v_im[bus1]*v_re[bus2];

var  q_or{(bus1,bus2) in branch} = 
	+y11_re[bus1, bus2]*v_re[bus1]*v_re[bus1]
	+y11_re[bus1, bus2]*v_im[bus1]*v_im[bus1]
	+y12_re[bus1, bus2]*v_re[bus1]*v_re[bus2]
	+y12_re[bus1, bus2]*v_im[bus1]*v_im[bus2]
	+y12_im[bus1, bus2]*v_im[bus1]*v_re[bus2]
	-y12_im[bus1, bus2]*v_re[bus1]*v_im[bus2];
	
var  q_ex{(bus1,bus2) in branch} = 
	+y22_re[bus1,bus2]*v_re[bus2]*v_re[bus2]
	+y22_re[bus1,bus2]*v_im[bus2]*v_im[bus2]
	+y21_re[bus1,bus2]*v_re[bus1]*v_re[bus2]
	+y21_re[bus1,bus2]*v_im[bus1]*v_im[bus2]
	+y21_im[bus1,bus2]*v_re[bus1]*v_im[bus2]
	-y21_im[bus1,bus2]*v_im[bus1]*v_re[bus2];
	
var reduced_cost =
 	-sum{bus in bus_i, (bus1, bus2) in branch:bus==bus1} dual_p[bus]*p_or[bus1,bus2]
	-sum{bus in bus_i, (bus1, bus2) in branch:bus==bus2} dual_p[bus]*p_ex[bus1,bus2]
	-sum{bus in bus_i, (bus1, bus2) in branch:bus==bus1} dual_q[bus]*q_or[bus1,bus2]
	-sum{bus in bus_i, (bus1, bus2) in branch:bus==bus2} dual_q[bus]*q_ex[bus1,bus2]
	-dual_convexity
	#-sum{bus in bus_i} (v_re[bus]^2+v_im[bus]^2)*(dual_vmin[bus]+dual_vmax[bus])
	;
subject to ctr_v_min{bus in bus_i}: v_re[bus]^2+v_re[bus]^2 >= Vmin[bus]^2 ;
subject to ctr_v_max{bus in bus_i}: v_re[bus]^2+v_re[bus]^2 <= Vmax[bus]^2;