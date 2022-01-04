% function [] = fast_dynamics(a, b2, Kint, )


nx=2;                   % Agent's dimensio
% nCk = [3 4 5 6];
nCk = [120 140 170 200];
m = length(nCk);
n = sum(nCk);
% nk = nCk(Ck);
theta = 30; 
a =1*[cos(theta),-sin(theta);sin(theta),cos(theta)];
b2=[1;1];     
p  =0.5;
[Lint,  Lext] = laplacianER(nCk, p);
[Lextbar,  H, W] = weighted_laplacian(nCk, Lext);
%% 
load ('Data-21-09-2021/init.mat')
% load ('Data-21-09-2021/Lext.mat')
% load ('Data-21-09-2021/Lint.mat')
% load ('Data-21-09-2021/Kext_opt.mat')
load ('Data-21-09-2021/Kint_opt.mat')

%% Laplacians of the Clusters
 [Lint1, ~] = laplacian_cluster(Lint, Lext, 1, nCk);
 [Lint2, ~] = laplacian_cluster(Lint, Lext, 2, nCk);
 [Lint3, ~] = laplacian_cluster(Lint, Lext, 3, nCk);
 [Lint4, ~] = laplacian_cluster(Lint, Lext, 4, nCk);

%% 
delta_t = 0.01;
tf=15;
tspan = 0:delta_t:tf;

%% Solve Cluster 1:
Ic_1 = eye(nCk(1)-1);
[V1,D1,W1] = eig(Lint1);
V1= V1./(V1(1,1));
Z1 = V1;
Z1(:,1) = [];
init_1 = init(1:2*nCk(1));
int_err_1 = kron(eye(2), Z1)'*init_1;
Cl1 = kron(Ic_1,a) - kron(nCk(1)*Ic_1, b2*Kint(1,:));
[t1, y1] = ode23(@(t1,y1) Cl1*y1, tspan, int_err_1);

figure
plot(t1, y1)
% end