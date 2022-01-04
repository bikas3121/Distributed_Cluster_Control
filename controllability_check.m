function controllability_check(a,b2)
cab  = ctrb(a,b2);  
if length(a) == rank(cab)
    disp('********************************************')
    disp('Controllable System ')
    disp('********************************************')
else
    disp('********************************************')
    disp('Uncontrollable System ')
    disp('********************************************')
end
pause(2)
end