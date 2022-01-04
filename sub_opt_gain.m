function [Kopt, eig_Pext_opt, eig_LMI_opt, K_iter] = sub_opt_gain (nCk, nx, Kext, Ck, a, b2, Lext, init, Rk, sigma)
m = length(nCk);
iter = 1;
 j = 1; 

Kext1 = Kext;
[Kext_k, Kext_minus_k]  = gain_gen (Kext1, Ck);
K_iter = [Kext_k];
 eig_Pext_iter = [];
 eig_LMI_iter = [];
 
while iter >0
       K_iter = [K_iter; Kext_k];
       [A_k_k, Q_k] = external_cost_matrices(nCk, nx, Ck, a, b2, Lext, init, Kext_k, Kext_minus_k, Rk);
       nP =(m-1)*nx;
        [valP, valF] = solverLMI(nP, sigma, A_k_k, Q_k);
        eig_Pext = eig(valP)';
        eig_LMI = eig(valF)';
        eig_Pext_iter = [eig_Pext_iter; eig_Pext];
       eig_LMI_iter = [eig_LMI_iter; eig_LMI];
         [sgn_eig] = err_eig_check(A_k_k);
        [sgn] = lmi_conditions(valP, valF);
        if sgn  == 1 && sgn_eig ==1
            Kext_k = j*Kext_k;
            j = j-0.1; 
        else
               break
        end
end
if length(K_iter) == 1
    ed = 1;
else
    ed = length(K_iter) -1;
end
Kval = K_iter(ed,:);
    for k = 1:2
        Kext1(Ck, k) = Kval(1,k);
    end
eig_Pext_opt = eig_Pext_iter(ed, :);
eig_LMI_opt = eig_LMI_iter(ed, :);    

Kopt = Kext1;


end
