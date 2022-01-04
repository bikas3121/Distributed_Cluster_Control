function [sgni] = lmi_conditions(valP, valF)


eig_Pext = eig(valP)';
eig_LMI = eig(valF)';


for i = 1:size(eig_LMI,2)
    sgn = eig_LMI(i);
    if sgn < 0
        sgni = 1;
    elseif sgn >= 0
        sgni = 0; 
    end
end

if sgni ~=1
    disp('************************************')
    disp('LMI Conditions Not Satisfied')                
    disp('************************************')
end
end
