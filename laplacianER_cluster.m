function [Lint] = laplacianER_cluster(nCk, Ck, p)
%Function [Lintk, Lint, Lext] = laplacianER(nCk, p, Ck) generates the
%complete internal Laplacian and the random external Laplacian of the graph
%using Erdos-Renyi Algorithm

% Inputs: 
% nCk - is the vectors containing the number of agents in each
% cluster. For example [3 4 5] means three clusters with 3,4 and 5 agetns
% p  - is the probablity for the Erdos-Reny Algorithm
% Ck  - represent the cluster i.e., Ck = 1 means Lintk in the output will
% be for the cluster 1.


%Outputs
%Lintk - Laplacian of the complete graph corresponding to cluster Ck
%Lint - Laplacian of the overall network while ignoring the external
%connections
% Lext - Random external Laplacian graph generated using ER algorithm.

% clear all
% clear vars
% clc

% 
% p =0.9;
% nCk = [7 4 5];
% 
% Ck = 1;
% m =length(nCk);

nk = nCk(Ck );
[G]=erdosRenyi(nk,p,1);
adjG =full (G.Adj);
adj = adjG;

% make adj symmetric
for  i =1:nk
    for j = 1:nk
        if adj(i,j) == 1
            adj(j,i) = 1;
        end
    end
end

% remove diagonal elements
for i = 1:nk
    adj(i,i) =0;
end
% 
% tf1 = issymmetric(adj)
% adj


% Generate the Degree matrix. 

deg = zeros(nk,nk); %degree matrix. 
for i = 1:nk
    for j = 1:nk
        if i == j
            deg(i,j) = sum(adj(i,:));
        end
    end
end
% deg
% 
Lint = deg - adj;
end