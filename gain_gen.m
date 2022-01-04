function [Kext_k, Kext_minus_k]  = gain_gen (Kext, Ck)

Kext_k = Kext(Ck,:);
Kext(Ck,:) = [];
Kext_minus_k = [];
for i = 1:length(Kext)
    Kext_minus_k = blkdiag(Kext_minus_k, Kext(i,:));
end
end