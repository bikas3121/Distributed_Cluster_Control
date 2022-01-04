function [Lintk, Lextk] = laplacian_cluster(Lint, Lext, Ck, nCk)

nk = nCk(Ck);
n = sum(nCk);

%% External Laplacian for Cluster Ck
Lextk = Lext;
if Ck == 1
    Lextk(nCk(Ck)+1:end,:) = 0;
     Lextk(nCk(Ck)+1:end,1:nk) =Lextk(1:nk,nCk(Ck)+1:end)';
     for i  = (nCk(Ck) +1): n
        if sum(Lextk(i,:)) ~=0
            Lextk(i,i) = -sum(Lextk(i,:));
        end
     end
else
    st = sum(nCk(1:Ck-1));
    ed = st+nCk(Ck)+1;
    Lextk(1:st,:) = 0;
    Lextk(ed:end,:) = 0;
     Lextk2 = Lextk';
 for i = 1:length(Lextk2)
        Lextk2(i,i) = 0;
 end
Lextk = Lextk + Lextk2;
    for i = 1:length(Lextk)
        if sum(Lextk(i,:)) ~= 0
            Lextk(i,i)  = -sum(Lextk(i,:));
        end
    end
end


%% Internal Laplacian for Ck
Lintk1 = Lint;
if Ck == 1
    Lintk = Lintk1(1:nCk(Ck), 1:nCk(Ck));
else
    st1 = sum(nCk(1:Ck-1));
    ed1= st1+nCk(Ck);
    Lintk = Lintk1(st1+1:ed1, st1+1:ed1);
end







% 
% Lintk = zeros(nk,nk);
% for i = 1:nk
%     for j = 1:nk
%         if i == j
%         Lintk(i,j) = nk-1;
%         else
%         Lintk(i,j) = -1;
%         end
%     end
% end

end