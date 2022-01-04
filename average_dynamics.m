function [t_av, y_av, y_av1, y_av2, y_av1_err, y_av2_err, J_avg_inst] = average_dynamics(nx, nCk, Ck, a, b2, Kext, Lext, init, tspan, delta_t, Q_k, Rk)
% This function [t_av, y_av, y_av1, y_av2] = average_dynamics(nx, nCk, a,
% b2, Kext, Lextbar, init, tspan) solves the average dynamics corresponding
% to the equation (34) in the paper. 

%% INPUTS:
% nCk - Arrary of the agents in the cluster,
% nx -Dimension of the system,
% a - state of each individual agent. all agents have same 'a' i.e.,
% network is homogeneous,
% b2 - input matrix for each agent,
% Lextbar - Weighted external Laplacian matrix obtained after
% transformation,
% Kext - External Gain that safisfies the LMI condition,
% init - initial condition,
% tspan - time span to solve ODE,

% OUTPUTS:
% t_av - time overwhich the dynamics is solved
% y_av - state value of the avg network dynamics
% y_av1 - first component of x
% y_av2 - second component of x
%  y_av1_err - error between components of y_av1, i.e, y_av1(1) -
%  y_av1(i)
%  y_av2_err - error between components of y_av2, i.e, y_av2(1) -
%  y_av2(i)
% J_avg_inst -  Instantaneous average cost in the equation 35 in the paper.
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = sum(nCk);  % total number of agents in the network
m = length(nCk);  % total number of clusters.
nk = nCk(Ck);
[Lextbar,  H, W] = weighted_laplacian(nCk, Lext);
A_av = kron(eye(m),a);
B_av = kron(eye(m),b2);
[Kext_k, Kext_minus_k]  = gain_gen (Kext, Ck);
Kbar_ext = blkdiag(Kext_k,Kext_minus_k);


initavg  =kron(W, eye(2))*kron(H,eye(2))'*init;

A_cl_av = (A_av - B_av*Kbar_ext*kron(Lextbar, eye(nx)));

[t_av, y_av] = ode15s(@(t,y) A_cl_av*y, tspan, initavg);

l = size(y_av,2);
y_av1 = [];  % first compnent of x_int
y_av2 = [];    % second compnent of x_int
for i=  2:2:l
    y_av1 = [y_av1,y_av(:,i-1)];
    y_av2 = [y_av2, y_av(:,i)];
end

y_av1_err = y_av1(:,1:m-1)- kron(ones(1,m-1),y_av1(:,m));
y_av2_err = y_av2(:,1:m-1)- kron(ones(1,m-1),y_av2(:,m));

% Q1 =nk*kron(Lextbark, eye(2));


U = [];
for i = 1:m 
    blk = ones(nCk(i),1);
    U = blkdiag(U, blk);
end



if Ck == 1
    ct_m = Lext(1:nCk(Ck),:);
else
    st = sum(nCk(1:Ck-1));
    ct_m = Lext(st+1: st+ nCk(Ck),:);
end

% Q2 = kron(U, eye(2))'*kron(ct_m'*ct_m, Kext'*Rk*Kext)*kron(U, eye(2));

%%  Instantaneous Average Cost -- Undiscounted

[~, t_avg, Y_k, ~, ~, ~, ~, ~] = average_error_dynamics(nx, nCk, Ck, a, b2, Lext, Kext, init, tspan, Rk, delta_t);

inst_avg_cost = zeros(length(t_avg),1);
index = 0;
for i = 1:length(t_avg)
    index = index +1;
    cost = Y_k(i,:)*Q_k*Y_k(i,:)' ; % undiscounted cost
    inst_avg_cost(index)  = delta_t* cost;
end
J_avg_inst = inst_avg_cost;

end
