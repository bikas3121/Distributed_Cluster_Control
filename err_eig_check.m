function [sgn_eig] = err_eig_check(A_k_k)

% A_k_k = [ 9.1825  -19.0354   -2.1146    7.2058;
%     7.3523  -12.2052   -2.1146    7.2058;
%    -2.1194    7.2223   11.3091  -26.2825;
%    -2.1194    7.2223    9.4790  -19.4524];

eig_A_k_k = real(eig(A_k_k));
 p = 0;
for i = 1:length(eig_A_k_k)
    if eig_A_k_k(i) <= -0.0000001
        p = p+1;
    end
end

if p == length(eig_A_k_k)
    sgn_eig = 1;
else
    sgn_eig  =0;
end

end