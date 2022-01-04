function [Lextbar, H, W] = weighted_laplacian(nCk, Lext)

% % % % The function [Lextbar] = weighted_laplacian(nCk, Lext) generates the
% % % % weighted Laplacian Lextbar(in paper). 
% % % % Inputs: 
% % % % nCk - vector of the agents in the clusters.
% % % % Lext - External Laplacian
% % % 
% % % % Outputs
% % % % Lextbar - Weighted Laplacian
% % % % H & W - tranformation matrices 


% clear
% clear vars 
% clear all
% clc



% p =1;
% Ck =1 ;
% nCk = [ 3 4 5];
% [Lint, Lext] = laplacianER(nCk,p, Ck) ; 

%%%%%%%%%%%%%
P = [];

for i = 1:length(nCk)
    nCki = nCk(i);
    P1 = ones(nCki,1);
    P = blkdiag(P,P1);
end


H = [];
for i = 1:length(nCk)
    nCki = nCk(i);
    H1 = ones(nCki,1)*sqrt(1/nCki);
    H = blkdiag(H,H1);
end

W = [];

for i = 1:length(nCk)
    nCki = nCk(i);
    W1 = sqrt(1/nCki);
    W = blkdiag(W,W1);
end


Lextbar=W*H'*(Lext)*H*W^(-1);

%% Weighted Laplacian for the cluster.

% Lextbark = Lextbar;
% 
% if Ck == 1
%     Lextbark(2:end,:) = 0;
%     Lextbark2 = Lextbark';
%     Lextbark2(Ck,Ck) = 0;
%     Lextbark = Lextbark + Lextbark2;
% for i = 1:length(Lextbark)
%     if sum((Lextbark(i,:))) <- 0.0000000001
%         Lextbark(i,i) = -sum(Lextbark(i,:));
%     end
% end
% else
%      Lextbark(1:Ck-1,:) =0;
%      Lextbark(Ck+1:end,:)=0;
%      Lextbark2 = Lextbark;
%      Lextbark2(Ck,Ck) =0;
%      Lextbark = Lextbark +Lextbark2';
% for i = 1:length(Lextbark)
%     if sum((Lextbark(i,:))) <- 0.0000000001
%         Lextbark(i,i) = -sum(Lextbark(i,:));
%     end
% end
% end


end