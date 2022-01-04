function [ J_tot] = total_cost1( J_k_inst_1,J_k_inst_2,J_k_inst_3,J_k_inst_4)



% The function [ J_tot] = total_cost1( J_k_inst_1,J_k_inst_2,J_k_inst_3,J_k_inst_4)
% calculates the cummulative cost from the instataneous cost 


%% Cluster Cost
J_tot_1 = cumsum(J_k_inst_1);
J_tot_2 = cumsum(J_k_inst_2);
J_tot_3 = cumsum(J_k_inst_3);
J_tot_4 = cumsum(J_k_inst_4);

%% Cost Difference
J_tot  = [J_tot_1(end); J_tot_2(end); J_tot_3(end); J_tot_4(end)];
end