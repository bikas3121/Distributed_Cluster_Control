function [Lint,  Lext] = laplacianER(nCk, p, Kreg)
%Function [Lint, Lext] = laplacianER(nCk, p, Ck) generates the
%complete internal Laplacian and the random external Laplacian of the graph
%using Erdos-Renyi Algorithm

% Inputs: 
% nCk - is the vectors containing the number of agents in each
% cluster. For example [3 4 5] means three clusters with 3,4 and 5 agetns
% p  - is the probablity for the Erdos-Reny Algorithm
% Ck  - represent the cluster i.e., Ck = 1 means Lintk in the output will
% be for the cluster 1.


%Outputs
%Lint - Laplacian of the overall network while ignoring the external
%connections
% Lext - Random external Laplacian graph generated using ER algorithm.

% clear all
% clear vars
% clc

% 
% p =0.3;
% nCk = [3 4 5];


m =length(nCk);
n = sum(nCk);
[G]=erdosRenyi(n,p,Kreg);

adjG =full (G.Adj);

% Remove the connections inside the clusters
adjG1 = ones(n,n);
for i = 1:m
    if i == 1
        inds =   1;
        inde = inds + nCk(i) -1;
         adjG1(inds:inde,inds:inde) = 0;
    else 
        nCkminus = nCk(1:i-1);
        snCkminus = sum(nCkminus);
        inds =  snCkminus + 1;
        inde = inds + nCk(i) -1;
         adjG1(inds:inde,inds:inde) = 0;
    end
end

adj = adjG.*adjG1;

% make adj symmetric
for  i =1:n
    for j = 1:n
        if adj(i,j) == 1
            adj(j,i) = 1;
        end
    end
end

% tf1 = issymmetric(adj)
% adj
%Generate the Degree matrix. 

deg = zeros(n,n); %degree matrix. 
for i = 1:n
    for j = 1:n
        if i == j
            deg(i,j) = sum(adj(i,:));
        end
    end
end
% deg

Lext = deg - adj;

%% Internal Laplacians
Lint = [];
for i = 1:length(nCk)
    nki =  nCk(i);
    Linti = zeros(nki,nki);
    for j = 1:nki
        for k = 1:nki
            if j ==k 
                Linti(j,k) = nki-1;
            else
                Linti(j,k) = -1;
            end
        end
    end
    Lint = blkdiag(Lint, Linti);
end

end