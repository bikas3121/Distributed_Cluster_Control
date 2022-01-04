function [Lint, Lext, Lint_norm] = lap_gen_conn (nCk, pint, pext)
% pext =0.2;
% pint = 0.9
% nCk = [3 4 5];

m =length(nCk);
n = sum(nCk);
[G]=erdosRenyi(n,pext,1);

adjG =full (G.Adj);

%% Remove the connections inside the clusters
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
% 
% tf1 = issymmetric(adj)
% adj

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

%% Internal Block Diagonal Laplacian

Lint = [];
Lint_norm = [];
for i = 1:length(nCk)
    [Lintk] = laplacianER_cluster(nCk, i, pint);
    Lint_norm = [Lint_norm, norm(Lintk,2) ];
    Lint = blkdiag(Lint, Lintk);
end
% Lint
end