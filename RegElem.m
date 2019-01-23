function [Tnew,Tout,dQ] = RegElem(Tin,Told,dmgas,dt)
%RegElem.m evaluates the new conditions in one element of the
%regenerator

%%%% Variables
% Tin and Told are the inlet termperature of the gas through the
% element (K) and the old temperature of metal in the element (K),
% dmgas is the mass of gas flowing through the element in one 
% time step (kg). 

%%%% Operation
% The new tmeperature of the gas is found. The heat transferred from
% gas to metal in the element in one time step is then found. Using 
% this heat transfer, the new temperature of the metal is found.

global UAr dmmet cpmet cpgas

%%%% TEMPERATURE OUT
A=dmmet*cpmet*Told+dmgas*cpgas*Tin;
B=dmmet*cpmet+dmgas*cpgas;
C=UAr/(dmgas*cpgas*dmmet*cpmet);
Tout=(A/B)-((A/B)-Tin)*exp(-B*C*dt);

%%%% HEAT TRANSFERRED
dQ=dmgas*cpgas*(Tout-Tin);

%%%% NEW TEMPERATURE OF METAL    
Tnew=Told-(dQ/(dmmet*cpmet));

end

