function [TC,Ts,Te,Vs,Ve,Tmet,P,xp,t2,dmev,m,dQi,dQo,dQR,...
                      dU,Uf] = Stage4(t,xp,P,TC,Ts,Te,Vs,Ve,Tmet,Ui)
% Stage4.m updates all temperatures and pressures inside the
% cylinder for the fourth stage of the cycle 

%   Stage 4 - inlet check valve open, piston in reverse half of its
%   revolution

global VelemC Ac Lc dxc nc nr om f rev
global Th Tc Tamb R cvgas cpgas UAe  

%% %%            Find new conditions of the system             %% %%

% kp indicates element occupied by piston
kpold=ceil(xp/dxc);
if kpold==0
    kpold=1;
end

% Find new piston position and time interval for step. Use full
% element volume unless piston occupies first element - in this
% case use partial volume after the piston
if kpold>1
    
    % Determine new piston position
    xp=xp-dxc;
    kpnew=ceil(xp/dxc);
    
    % Find time interval
    t2=(rev/f)-(acos((Lc-2*xp)/Lc)/om);
    dt=t2-t;

    % Determine mass through reg and heater
    dmgas=(P*VelemC)/(R*TC(1));
    
    % Run through heater and regenerator directly
    % Heater
    Tout=Th-(Th-TC(1))*exp((-UAe*dt)/(dmgas*cpgas));
    dQi=dmgas*cpgas*(Tout-TC(1));
    % Regenerator
    for i=nr:-1:1
        [Tmet(i),Tout,dQR(i)]=RegElem(Tout,Tmet(i),dmgas,dt);
    end
    % Cooler
    Tin=Tout;
    Tout=Tamb+(Tin-Tamb)*exp((-UAe*dt)/(dmgas*cpgas));
    dQo=dmgas*cpgas*(Tout-Tin);
    
    dVnew=dmgas*R*Tout/P;
    dxnew=dVnew/Ac;
    
    % Find mass of gas enetering system
    dmev=(dmgas/Tc)*(TC(1)-Tout);
    dVev=dmev*R*Tc/P;
    dxev=dVev/Ac;

    % Update temperatures by left shift
    TC(1:kpnew-1)=TC(2:kpold-1);
    TC(kpnew+1:nc-1)=TC(kpold+1:nc);
    TC(nc)=(Tout*dxnew+Tc*dxev)/dxc;
else
    
%%%% Case where piston occupies final element
    % Set piston position
    xp=0;
    kpnew=1;
    dxs=Vs/Ac;
    
    % Find time interval
    t2=rev/f;
    dt=t2-t;
    
    % Determine mass through reg and heater
    dmgas=(P*Ve)/(R*Te);
    
    % Run through heater and regenerator directly
    % Heater
    Tout=Th-(Th-Te)*exp((-UAe*dt)/(dmgas*cpgas));
    dQi=dmgas*cpgas*(Tout-Te);
    % Regenerator
    for i=nr:-1:1
        [Tmet(i),Tout,dQR(i)]=RegElem(Tout,Tmet(i),dmgas,dt);
    end
    % Cooler
    Tin=Tout;
    Tout=Tamb+(Tin-Tamb)*exp((-UAe*dt)/(dmgas*cpgas));
    dQo=dmgas*cpgas*(Tout-Tin);
    
    dVnew=dmgas*R*Tout/P;
    dxnew=dVnew/Ac;
    
    % Find mass of gas entering system
    dmev=(dmgas/Tc)*(Te-Tout);
    dVev=dmev*R*Tc/P;
    dxev=dVev/Ac;
    
    T2=TC;
    
    % Update temperature by simple mix
    Ts=dxc/(dxs/Ts+(dxc-dxs)/T2(2));
    for i=2:nc-1
        TC(i)=dxc/(dxs/T2(i)+(dxc-dxs)/T2(i+1));
    end
    TC(nc)=dxc/(dxs/T2(nc)+dxnew/Tout+dxev/Tc);
    Ve=0;
    Vs=VelemC;
    
end

%% %%    Check for mass and internal energy conservation       %% %%

% Find system mass
m=(P*VelemC/R)*sum(1./TC(:,1:kpnew-1));
m=m+((P*Ve)/(R*Te))+((P*Vs)/(R*Ts));
m=m+(P*VelemC/R)*sum(1./TC(:,kpnew+1:nc));

% Find change in internal energy
Uf=(P*cvgas/R)*((nc-1)*VelemC+Ve+Vs);
dU=Uf-Ui;
dQR=sum(dQR);

end

