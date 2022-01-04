function [J_k_inst] = original_cost(nCk, Ck, Lint, Lext, x, Kext, Kint, Rk, delta_t )

% The Function function [J_k_inst] = original_cost(nCk, Ck, Lint, Lext, x,
% Kext, Kint, Rk, delta_t )  calculates the orignal  cost corresponding to
% network dynamics for the cluster Ck. It corresponds to the equation (5)
% in the paper. 

% Inputs:
% nCk - vectors containing the number of agents in the cluster.
% Ck - represents the corresponding cluster.
% a - internal dynamics : state matrix of an agent
% b2 - internal dynamics : input matrix of an agent
% Lint - Internal Laplacian of the network 
% Lext - External Laplacian of the network 
% Kint - The optimal  internal gain obtained from LQR
% Kext - External gain that satisfies the LMI (eq 48 in the paper)
% R_k - Weighing matrix in the original cost
% init & tspan - Initial values and the time span for the solving the odes
% x - the states obtained after solving the network dynamics in equation
% (14) in the paper. 


% Outputs:
% t - time overwhich the dynamics is solved
% x - state value of the network dynamics
% x_1 - first component of x
% x_2 - second component of x
%  x_1_err - error between components of x_1, i.e, x_1(1) -
%  x_1(i)
%  x_2_err - error between components of x2_int, i.e, x_2(1) -
%  x_2(i)
% total_cost -  cummulative cost of the cluster Ck (eq 29 in the paper)


%% Test Parameters

% theta = 30; 
% a =1*[cosd(theta),-sind(theta);sind(theta),cosd(theta)];
% b2=[1;1]; 
% p = 0.3;
% Kext = [-10, 10; -10, 10; -10, 10];
% Kint = [-10, 10; -10, 10; -10, 10];
% nCk = [3 4 5];
% [Lint,  Lext] = laplacianER(nCk, p);
% init =  initial_conditions(nCk);
% 
% Ck = 1;
% 
% % [Lintk, Lextk] = laplacian_cluster(Lint, Lext, Ck, nCk);
%  delta_t = 0.1;
% tspan = 0:delta_t:1;
% [t, x, x_1_err, x_2_err]= network_dynamics (nCk,  a, b2, Lint, Lext, Kext, Kint, init, tspan);
% 
% m = length(nCk);
% % Weight Matrix Rk
% % Rk = [1/nCk(1), 1/nCk(2), 1/nCk(3), 1/nCk(4)];
% R = 0.01;
% Rk = [];
% for i = 1:m
%     Rk = [Rk, R];
% end


%% 
nx = 2;
% n = sum(nCk);
nk = nCk(Ck);
[Lintk, Lextk] = laplacian_cluster(Lint, Lext, Ck, nCk);
%% Internal Laplacian of the Cluster
L_int_k = [];
for i = 1:length(nCk) 
    if Ck == i
        L_int_k = blkdiag(L_int_k, Lintk);
    else
        Lzero = zeros(nCk(i), nCk(i));
        L_int_k = blkdiag(L_int_k, Lzero);
    end
end

%%
Q = L_int_k + Lextk;
len = size(x,1);
J_1 = zeros(len,1);
for i = 1:len
    J_1_i = x(i,:)*kron(Q, eye(nx))*x(i,:)';   % undiscounted cost
  J_1(i) = delta_t*J_1_i;
end


%% External Control 

if Ck == 1
    Lk_red = Lext(1:nCk(Ck),:);
else
    st = sum(nCk(1:Ck-1));
    Lk_red = Lext(st+1: st+ nCk(Ck),:);
end

uk_ext = [];
for i = 1:len
    uk_ext_i = - kron(Lk_red,Kext(Ck,:))*x(i,:)';
    uk_ext=[uk_ext; uk_ext_i'];
end

%% Internal Control


uk_int = [];
for i = 1:len
      uk_int_i1 = kron(L_int_k,eye(2))*x(i,:)';
    uk_int_i = zeros(nk*nx,1);
    indx = 0;
     for j = 1:length(uk_int_i1)
         if uk_int_i1(j) ~= 0
             indx = indx+1;
             uk_int_i(indx) = uk_int_i1(j);
         end
     end
     uk_int_i = (- kron(eye(nk),Kint(Ck,:)))*uk_int_i;
     uk_int = [uk_int; uk_int_i'];
end


uk = uk_int + uk_ext;


J_2  = zeros(len,1);
for i = 1:len
    J_2_i = uk(i,:)*kron(eye(nk),Rk(Ck))*uk(i,:)';  % undiscounted cost
    J_2(i) = delta_t*J_2_i;
end

J_k_inst = J_1 + J_2;

%% Cumulative Cost
% J_k_cum = cumsum(J_k_inst);

end