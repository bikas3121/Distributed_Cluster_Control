function  [val_count] = initial_gain_check(nCk, nx, a, b2, Lext, init, Kext, Rk, sigma)
m  = length(nCk);
val_count = 0;
for i = 1: length(nCk)
[Kext_k, Kext_minus_k]  = gain_gen (Kext, i);
 [A_k_k, Q_k] = external_cost_matrices(nCk, nx, i, a, b2, Lext, init, Kext_k, Kext_minus_k, Rk);
 nP =(m-1)*nx;
[valP, valF] = solverLMI(nP, sigma, A_k_k, Q_k);
 [sgn] = lmi_conditions(valP, valF);
 if sgn == 1
     val_count = val_count +1;
 else
     val_count = val_count;
 end
end
end