function [init] =  initial_conditions(nCk)
% The function [init] =  initial_conditions(nCk, Ck) generates the random
% initial conditions in a evenly spaced manner,

% Inputs:
% nCk - Vectors of the agents in the clusters
% Ck - Cluster label

% Outputs:
% init - vector of the initial state of the agetns.

% nCk = [3 4 5 6 7];
m =length(nCk);

rv = mod(m,2);
if rv == 1
    median_m = ceil(m/2) ;
else
    median_m = m/2;
end
rg = 30;
init = [];
for i = 1:median_m
%         b = rand(1)* rand(1);
%         b = rg*i;
        b = (rg)*i;
    initi = b + (7)*rand(2*nCk(i),1);
    init = [init; -initi];
end
% for j = median_m+1:m
%     a = rand(1)*rand(1);
% %     a = 1*rand(1)*j;
%     initj = -a + ((10+a)-a)*rand(2*nCk(j),1);
%     init = [init; initj];
%     
% end

for j = median_m+1:m
%     a = 50*(j-median_m)
%     b = rg*(j-median_m);
    b = (rg)*(j-median_m);
        initi = b + (7)*rand(2*nCk(j),1);
    init = [init; initi];
end

% init
end