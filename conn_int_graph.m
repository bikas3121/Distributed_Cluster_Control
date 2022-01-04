function [L_con, L_norm] =  conn_int_graph(Lint, nCk)
% The function [L_con, L_norm] =  conn_int_graph(Lint, nCk) generates the
% connected internal graph. It takes complete Laplacian and the number of
% agents in each cluster as the input and generates the connected internal
% Laplaician matrix. 

%Inputs:
% Lint - Complete Internal Laplacian of the network. It is a block diagonal
% matrix. 
%nCk - vector containing the number of agents in each cluster. 

%Output
% L_con - The connected internal block diagonal matrix.
% L_norm - Norm of the matrix L_con.


%% Test Parameters
% nCk = [6 7 9 8];
% nCk = [120 140 170 200];
%% Generate Lint Sparse
% L_norm = [];
Lint_sparse = [];
for i = 1:length(nCk)
    nk = nCk(i);
    [G]=erdosRenyi(nk,0.1,40);
adj =full (G.Adj);
% make adj symmetric
for  i =1:nk
    for j = 1:nk
        if adj(i,j) == 1
            adj(j,i) = 1;
        end
    end
end

deg = zeros(nk,nk); %degree matrix. 
for i = 1:nk
    for j = 1:nk
        if i == j
            deg(i,j) = sum(adj(i,:));
        end
    end
end
L = deg- adj;
Lint_sparse = blkdiag(Lint_sparse, L);
end

% Generate the connected Laplacian matrix by deleting some connections from
% complete internal Laplacian matrix. 
L_con = Lint- Lint_sparse;



% Calculates the norm of the internal Laplacians
Lintk1 = L_con;
L_norm = [];
for i = 1:length(nCk)
if i == 1
    Lintk = Lintk1(1:nCk(i), 1:nCk(i));
    L_norm_i = norm(Lintk, 2);
    L_norm = [L_norm,  L_norm_i];
else
    st1 = sum(nCk(1:i-1));
    ed1= st1+nCk(i);
    Lintk = Lintk1(st1+1:ed1, st1+1:ed1);
    L_norm_i = norm(Lintk, 2);
    L_norm = [L_norm,  L_norm_i];
end
end
% L_con
% L_norm

%% Connected Laplacian

end