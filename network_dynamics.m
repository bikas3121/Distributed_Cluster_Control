function [t, x, x_1, x_2, x_1_err, x_2_err]= network_dynamics (nCk,  a, b2, Lint, Lext, Kext, Kint, init, tspan)

% The Function [t, x, x_1, x_2, x_1_err, x_2_err, total_cost] = network_dynamics(nx, nCk, Ck, a, b2,Lint, Lintk, Lextk, Lext, Kint, Kext, Rk, init, tspan)
% calculates the network dynamics and the cost related to the cluster Ck.

% Inputs:
%  nx- dimension of the dynamiccal system
% nCk - vectors containing the number of agents in the cluster.
% Ck - represents the corresponding cluster.
% a - internal dynamics : state matrix of an agent
% b2 - internal dynamics : input matrix of an agent
%  Lintk - Internal Laplacian of the the cluster Ck
%  Lextk - Exernal Laplacian of the the cluster Ck
% Lint - Internal Laplacian of the network 
% Lext - External Laplacian of the network 
% Kint - The optimal  internal gain obtained from LQR
% Kext - External gain that satisfies the LMI (eq 48 in the paper)
% R_k - Weighing matrix in the original cost
% init & tspan - Initial values and the time span for the solving the odes


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


%% Test Paramters
% nx =2;
% nCk = [3 4 5];
% Ck  = 3;
% theta = 30; 
% a = 5*[cosd(theta),-sind(theta);sind(theta),cosd(theta)];
% b2=[1;1];     
% p =0.2;
% [Lint,  Lext] = laplacianER(nCk, p);
% Kext = 1*[-2.1535    6.8305; -2.3408    7.4244; -2.9360    9.3123];
% Kint = 5*[-2.8471    9.2071; -2.4832    7.5260; -2.2966    6.5711];

% load ('initial/init1.mat')
% delta_t = 0.01;
% tf=3;
% tspan = 0:delta_t:tf;

%% Network parameters
n = sum(nCk);
m = length(nCk);
% nk = nCk(Ck);
nx= 2;
% [Lintk, Lextk] = laplacian_cluster(Lext, Ck, nCk);
%% Network Dynamics
org_A = kron(eye(n),a);
org_B = kron(eye(n),b2);

K_org_ext = [];
for i = 1:length(nCk)
    K_blk1 = kron(eye(nCk(i)), Kext(i,:)); 
    K_org_ext = blkdiag( K_org_ext,  K_blk1);
end

K_org_int = [];
for i = 1:length(nCk)
    K_blk2 = kron(eye(nCk(i)), Kint(i,:));
    K_org_int = blkdiag( K_org_int,  K_blk2);
end


%Closed Loop dynamics
A_cl = (org_A - org_B*K_org_int*kron(Lint, eye(nx))-org_B*K_org_ext*kron(Lext,eye(nx)));

% ODE


[t, x] = ode15s( @(t,x) A_cl*x , tspan, init);

% Seprate the components of x into x_1 and x_2
l = size(x,2);
x_1 = [];  %component 1 of the x
x_2 = [];   % component 2 of the x
for i=  2:2:l
    x_1 = [x_1 x(:,i-1)];
    x_2 = [x_2 x(:,i)];
end


x_1_err = x_1(:,1:n-1)- kron(ones(1,n-1),x_1(:,n));
x_2_err = x_2(:,1:n-1)- kron(ones(1,n-1),x_2(:,n));

end































