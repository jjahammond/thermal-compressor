function [TC,Ts,Te,Vs,Ve,Tmet,P,xp,t2,dmcond,Tcond,m,dQi,dQo,dQR,...
                      dU,Uf] = Stage2(t,xp,P,TC,Ts,Te,Vs,Ve,Tmet,Ui)
% Stage2.m updates all temperatures and pressures inside the
% cylinder for the second stage of the cycle

%   Stage 2 - outlet check valve open, piston in forwards half of
%   revolution

global om f rev
global VelemC dxc nc nr Ac Lc
global Th Tamb R cvgas cpgas UAr UAe dmmet cpmet

%% %%            Find new conditions of the system             %% %%

% kp indicates element occupied by piston
kpold=ceil(xp/dxc);
if kpold==0
    kpold=1;
end

% Find new piston position and time interval for step. Use full
% element volume unless piston occupies final element - in this
% case use partial volume after the piston
if kpold<nc
    xp=xp+dxc;
    kpnew=ceil(xp/dxc);

    t2=((rev-1)/f)+acos((Lc-2*xp)/Lc)/om;
    dt=t2-t;

%%%% Method of Bisections to find temperature out of regenerator
    dm=(P*VelemC)/(R*TC(nc));
    dmmax=dm;
    dmmin=0;
    H1=dmmax*TC(nc);
    H2=0;
    % Iterate until enthalpies converge
    while abs(H1-H2)>1e-10
        % Set mass to analyse
        dmgas=(dmmax+dmmin)/2;
        % Cooler (use temporary variable through interations)
        Totemp=Tamb+(TC(nc)-Tamb)*exp((-UAe*dt)/(dmgas*cpgas));
        % Regenerator
        for i=1:nr
            A=dmmet*cpmet*Tmet(i)+dmgas*cpgas*Totemp;
            B=dmmet*cpmet+dmgas*cpgas;
            C=UAr/(dmgas*cpgas*dmmet*cpmet);
            Totemp=(A/B)-((A/B)-Totemp)*exp(-B*C*dt);   
        end
        % Heater
        Totemp=Th-(Th-Totemp)*exp((-UAe*dt)/(dmgas*cpgas));

        % Enthalpy of gas through reg and heater
        H2=dmgas*Totemp;

        % Conditions to update mass analysed
        if H2>H1
            dmmax=dmgas;
        elseif H2<H1
            dmmin=dmgas;
        end
    end
    
    % Run through heater and regenerator again storing results
    Tout=Tamb+(TC(nc)-Tamb)*exp((-UAe*dt)/(dmgas*cpgas));
    dQo=dmgas*cpgas*(Tout-TC(nc));
    for i=1:nr
        [Tmet(i),Tout,dQR(i)]=RegElem(Tout,Tmet(i),dmgas,dt);
    end
    Tin=Tout;
    Tout=Th-(Th-Tin)*exp((-UAe*dt)/(dmgas*cpgas));
    dQi=dmgas*cpgas*(Tout-Tin);
    
    % Find mass and temp of gas leaving system
    dmcond=dm-dmgas;
    Tcond=TC(nc);

    % Update temperatures by right shift
    TC(kpnew+1:nc)=TC(kpold+1:nc-1);
    TC(2:kpnew-1)=TC(1:kpold-1);
    TC(1)=Tout;
else
    
%%%% Case where piston occupies final element
    xp=Lc;
    kpnew=nc;
    dxh=Vs/Ac;
    
    t2=((rev-1)/f)+acos((Lc-2*xp)/Lc)/om;
    dt=t2-t;
    
%%%% Method of Bisections to find temperature out of regenerator    
    dm=(P*Vs)/(R*Ts);
    dmmax=dm;
    dmmin=0;
    H1=dmmax*Ts;
    H2=0;
    % Iterate until enthalpies converge
    while abs(H1-H2)>1e-10
        % Set mass to analyse
        dmgas=(dmmax+dmmin)/2;
        % Cooler (use temporary variable through interations)
        Totemp=Tamb+(Ts-Tamb)*exp((-UAe*dt)/(dmgas*cpgas));
        % Regenerator
        for i=1:nr
            A=dmmet*cpmet*Tmet(i)+dmgas*cpgas*Totemp;
            B=dmmet*cpmet+dmgas*cpgas;
            C=UAr/(dmgas*cpgas*dmmet*cpmet);
            Totemp=(A/B)-((A/B)-Totemp)*exp(-B*C*dt);   
        end
        % Heater
        Totemp=Th-(Th-Totemp)*exp((-UAe*dt)/(dmgas*cpgas));

        % Enthalpy of gas through reg and heater
        H2=dmgas*Totemp;

        % Conditions to update mass analysed
        if H2>H1
            dmmax=dmgas;
        elseif H2<H1
            dmmin=dmgas;
        end
    end
    
    % Run through heater and regenerator again storing results
    Tout=Tamb+(Ts-Tamb)*exp((-UAe*dt)/(dmgas*cpgas));
    dQo=dmgas*cpgas*(Tout-Ts);
    for i=1:nr
        [Tmet(i),Tout,dQR(i)]=RegElem(Tout,Tmet(i),dmgas,dt);
    end
    Tin=Tout;
    Tout=Th-(Th-Tin)*exp((-UAe*dt)/(dmgas*cpgas));
    dQi=dmgas*cpgas*(Tout-Tin);
    
    % Find mass and temp of gas leaving system
    dmcond=dm-dmgas;
    Tcond=Ts;
    T2=TC;
    
    % Update temperatures by simple mix
    TC(1)=dxc/(dxh/Th+(dxc-dxh)/T2(1));
    for i=2:nc-1
        TC(i)=dxc/(dxh/T2(i-1)+(dxc-dxh)/T2(i));
    end
    Te=dxc/(dxh/T2(nc-1)+(dxc-dxh)/Te);
    Ve=VelemC;
    Vs=0;
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

