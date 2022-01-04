function [A_k_k, Q_k] = external_cost_matrices(nCk, nx, Ck, a, b2, Lext, init, Kext_k, Kext_minus_k, Rk)
% The function [Q1_ck_ext, Q2_ck_ext] = external_cost_matrices(nCk, nx, Ck,
% Lextbar, Lext, Kext, R_k) calculates the weight matrices used in the
% cost function of the external error dynamics.

%% Inputs:
% nCk - Arrary of the agents in the cluster.
% nx -Dimension of the system
% Ck - Label corresponding to the cluster k.
% Lextbar - Weighted external Laplacian matrix obtained after transformation
% Lext - External Laplacian matrix
% Kext - External Gain that safisfies the LMI condition
% R_k - Weight matrix used in the original cost in 'x' variable. (Rk used
% in the equation (5) in the paper. 

%% Outputs:
% Q1_ck_ext - Weighing matrix Q1^ext
% Q2_ck_ext - Weighing matrix Q2^ext

%% Test Inputs:


% clear
% clearvars
% a = 4*[0,-0.25;0.25,0];
% b2=[1;0];
% nx =2 ;
% Ck = 2;
% nCk = [3 4 5];
% % nk = nCk(Ck);
% % n = length(nCk);
% Kextk = [ 1.8, 1.1];
% init =  initial_conditions(nCk);
% p =0.5;
% [ Lint, Lext] = laplacianER(nCk,p) ; 
% 
% Kext = [1 2; 3 4; 5  6];
% [Kext_k, Kext_minus_k]  = gain_gen (Kext, Ck);
% R_k = eye(1)*(nk^(-1));

%% State Matrices
% Kext = Kext(Ck, :);
nk = nCk(Ck);
n = sum(nCk);  % total number of agents in the network
m = length(nCk);  % total number of clusters.
[Lextbar, H, W] = weighted_laplacian(nCk, Lext);


initavg  =kron(W, eye(2))*kron(H,eye(2))'*init;  % average of the initial values (init), which is the initial conditions for the bar_y variable
% initavg
% Initial Values for the average errro dynamics=
if Ck == 1
    init_Ck = initavg(1:2,:);
    init_Ck_minus = initavg;
    init_Ck_minus(1:2,:) = [];
else
    st  = Ck*nx-1;
    en = st +1;
    init_Ck = initavg(st:en,:);
    init_Ck_minus = initavg;
    init_Ck_minus(st:en,:) = [];
end
% init_Ck_minus
init_avg_err = init_Ck_minus- kron(ones(m-1,1), init_Ck);
% init_avg_err


%% Gain
% In the average error dynamics we are applying the same control Kext
% obtained from the Jihene's algorithm. In case the different gain is
% applied, the Kext should be changed in the fnction
% 'average_error_dynamics'
% Kext_minus_k= kron(eye(m-1), Kext)     ;
% Lextbar

%% Laplacian Graph
Lextbar_minus_k = Lextbar;
Lextbar_minus_k(Ck,:) = [];
Lextbar_minus_k(:,Ck) = [];

Lextbar_k = Lextbar(Ck,:);
Lextbar_k(:,Ck) = [];
Lextbar_k = - Lextbar_k;
%% 

% External Error Dynamics 
A_k = kron(eye(m-1), a) - kron(eye(m-1), b2)*Kext_minus_k*kron(Lextbar_minus_k,eye(nx));
B_k = -kron(ones((m-1),1), b2);

%closedloop error dynamics 
A_k_k = (A_k + B_k*Kext_k*kron(Lextbar_k,eye(nx)));



%% Cost Matrices
%% Q2_ext for  cluster Ck
Lew  = Lextbar;

% Build the Qk for the given ck
Lew_ck = -Lew(Ck,:); %extracts the row from the Laplacian correspoding to the ck-th cluster
Lew_ck(Ck) = [];  % removes the ck-th element form the row
% Q1_ck_ext = nk*kron(diag(Lew_ck),eye(nx));
Q1_ck_ext = kron(diag(Lew_ck),eye(nx));

%% Q2_ext for  cluster Ck
nCk_minus_k = nCk;
nCk_minus_k(Ck) = [];

% Matrix Uk and U_minus_k
U_minus_k = [];
for i = 1:length(nCk_minus_k) 
    Uk = kron(ones(nCk_minus_k(i),1),eye(nx));
    U_minus_k = blkdiag(U_minus_k,Uk);
end

% Extract L_k_red from the external Laplacian Lext corresponding to the
% cluster Ck.

if Ck == 1
    n1 = nCk(1);
    L_k_red = Lext(1:n1,:);
    L_k_red(:,1:n1) = [];
else
        nCkminus = nCk(1:Ck-1);
        snCkminus = sum(nCkminus);
        inds =  snCkminus + 1;
        inde = inds + nCk(Ck) -1;
        L_k_red = Lext(inds:inde,:);
        L_k_red(:,inds:inde) = [];
end
R_k = Rk(Ck);
Q2_ck_ext =U_minus_k'*kron(L_k_red'*L_k_red, Kext_k'*R_k*Kext_k)*U_minus_k;

Q_k = Q1_ck_ext+(1/nk)*Q2_ck_ext;

end







