function [P_eps,Gain_Jih,G] = Gain_Ji_Jihene(a, b2, Lextbar,R)
%Cout Ji basé sur la méthode de Jihene
% global a b2
% global Lextbar 


%%Calcul du gain de Jihene
A=a;               %Dynamique d'un agent
vp_L=eig(Lextbar); %Valeur propre du Laplacien

eig_Lextbar = [vp_L]; % save the eigenvalues in different array
eig_Lextbar = sort(eig_Lextbar,'ascend');
eig_Lextbar(1) = [];

lambdaStar=min(eig_Lextbar) ;%Lambda star
Q=max(eig_Lextbar);         %Lambda circle
% 
B=sqrt(2)*lambdaStar*b2;%Matrice de commande
%  R=1;
[P_eps,~,G]=icare(A,B,Q,R,[],[],[]);%Resolution LQR
Gain_Jih =(1/R)*b2'*P_eps;
% Gain_Jih=inv(R)*B'*P_eps;%Gain Methode de Jihene


% B=sqrt(2)*lambdaStar*b2;%Matrice de commande
% R = 1;
% [P_eps,~,G]=icare(A,B,Q,R,[],[],[]);%Resolution LQR
% Gain_Jih=lambdaStar*inv(R)*b2'*P_eps;%Gain Methode de Jihene



end



