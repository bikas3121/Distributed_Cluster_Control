function [J_avg_1, J_avg_2, J_avg_3, J_avg_4, J_tot_1, J_tot_2, J_tot_3, J_tot_4, J_tot, J_avg, delta_J, delta_J_per ] = total_avg_cost( nCk, J_avg_inst_1, J_avg_inst_2, J_avg_inst_3, J_avg_inst_4, J_org_int_1, J_org_int_2, J_org_int_3,J_org_int_4, index_T)

%% Total cum cost
J_avg_1 =nCk(1)*cumsum(J_avg_inst_1(index_T:end,:));
J_avg_2 = nCk(2)*cumsum(J_avg_inst_2(index_T:end,:));
J_avg_3 = nCk(3)*cumsum(J_avg_inst_3(index_T:end,:));
J_avg_4 = nCk(4)*cumsum(J_avg_inst_4(index_T:end,:));


%% Cluster Cost
J_tot_1 = cumsum(J_org_int_1(index_T:end,:));
J_tot_2 = cumsum(J_org_int_2(index_T:end,:));
J_tot_3 = cumsum(J_org_int_3(index_T:end,:));
J_tot_4 = cumsum(J_org_int_4(index_T:end,:));

%% Cost Difference
J_tot  = [J_tot_1(end); J_tot_2(end); J_tot_3(end); J_tot_4(end)];
J_avg = [J_avg_1(end); J_avg_2(end); J_avg_3(end); J_avg_4(end)];

delta_J =abs(J_tot-J_avg);
delta_J_per = (delta_J./J_tot)*100;
end