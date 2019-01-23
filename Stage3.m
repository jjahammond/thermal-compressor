function [TC,Ts,Te,Vs,Ve,Tmet,P2,xp,t2,mtot,dQi,dQo,dQR,...
                dU,Uf] = Stage3(t,dtmax,xp,P,TC,Ts,Te,Vs,Ve,Tmet,Ui)
% Stage3.m updates all temperatures and pressures inside the
% cylinder for the third stage of the cycle 

%   Stage 3 - both check valves closed, piston in reverse half of
%   revolution

global Vtot VelemC Ac Lc dxc nc nr om rev f
global Th Tamb R y cvgas cpgas UAr UAe dmmet cpmet

%% %%             Find new pressure of the system              %% %%

% kp indicates element occupied by piston
kpold=ceil(xp/dxc);     
if kpold==0
    kpold=1;
end
if kpold>nc
    kpold=nc;
end

% Find mass of gas in first element and new pressure when passed 
% through the heater. Use full element volume unless piston occupies
% final element - in this case use partial volume after the piston
if kpold>1
    dmgas=(P*VelemC)/(R*TC(1));     % Mass of heated gas

%%%% Method of Bisections to find temperature out of regenerator
    dtmin=0;
    xp1=0;
    xp2=1;
    % Iterate until displacements converge
    while abs(xp2-xp1)>1e-10
        % Set initial time interval to analyse
        dt=(dtmin+dtmax)/2;
        % Heater (use temporary variable through interations)
        Totemp=Th-(Th-TC(1))*exp((-UAe*dt)/(dmgas*cpgas));
        % Regenerator
        for i=nr:-1:1
            A=dmmet*cpmet*Tmet(i)+dmgas*cpgas*Totemp;
            B=dmmet*cpmet+dmgas*cpgas;
            C=UAr/(dmgas*cpgas*dmmet*cpmet);
            Totemp=(A/B)-((A/B)-Totemp)*exp(-B*C*dt);
        end
        % Cooler
        Totemp=Tamb+(Totemp-Tamb)*exp((-UAe*dt)/(dmgas*cpgas));

        % Determine new pressure
        P2=fzero(@(P2) Vtot-(Vtot-VelemC)*(P2/P)^(-1/y)-...
                                (dmgas*R*Totemp/P2),P);

        % Displacement resulting from compression
        xp1=(P2/P)^(-1/y)*((kpold-2)*VelemC+Ve)/Ac;

        % Displacement resulting from time interval alone
        t2=t+dt;
        xp2=Lc*(1-cos(om*t2))/2;

        % Conditions to update interval analysed
        if xp1<xp2
            dtmin=dt;
        elseif xp1>xp2
            dtmax=dt;
        end
    end
    
    % Run through heater and regenerator again storing results
    Tout=Th-(Th-TC(1))*exp((-UAe*dt)/(dmgas*cpgas));
    dQi=dmgas*cpgas*(Tout-TC(1));
    for i=nr:-1:1
        [Tmet(i),Tout,dQR(i)]=RegElem(Tout,Tmet(i),dmgas,dt);
    end
    Tin=Tout;
    Tout=Tamb+(Tin-Tamb)*exp((-UAe*dt)/(dmgas*cpgas));
    dQo=dmgas*cpgas*(Tout-Tin);
    
else
    
%%%% Case where piston occupies final element
    dmgas=(P*Ve)/(R*Te);
    
    % Find time interval
    t2=rev/f;
    dt=t2-t;
    
    % Run through heater and regenerator directly
    Tout=Th-(Th-Te)*exp((-UAe*dt)/(dmgas*cpgas));
    dQi=dmgas*cpgas*(Tout-Te);
    for i=nr:-1:1
        [Tmet(i),Tout,dQR(i)]=RegElem(Tout,Tmet(i),dmgas,dt);
    end
    Tin=Tout;
    Tout=Tamb+(Tin-Tamb)*exp((-UAe*dt)/(dmgas*cpgas));
    dQo=dmgas*cpgas*(Tout-Tin);

    % Determine new system pressure
    P2=fzero(@(P2) Vtot-(Vtot-Ve)*(P2/P)^(-1/y)-...
                            (dmgas*R*Tout/P2),P);
end


%% %%             Update temperatures and volumes              %% %%

% Temperatures of variable volumes
T2(2:kpold-1)=TC(2:kpold-1)*((P2/P)^((y-1)/y));
T2e=Te*((P2/P)^((y-1)/y));      % Partial element before piston
T2s=Ts*((P2/P)^((y-1)/y));      % Partial element after piston
T2(kpold+1:nc)=TC(kpold+1:nc)*((P2/P)^((y-1)/y));

% Variable volumes
V2=VelemC*(P2/P)^(-1/y);        % General variable volume
V2e=Ve*(P2/P)^(-1/y);           % Partial volume before piston
V2s=Vs*(P2/P)^(-1/y);           % Partial volume after piston
dVnew=dmgas*R*Tout/P2;          % New volume of heated gas


%% %%       Find position of all elements within cylinder      %% %%

% Variable volume lengths
dxVar=V2/Ac;        % General variable volume length
dx2end=V2e/Ac;      % Partial volume before piston length
dx2start=V2s/Ac;    % Partial volume after piston length
dxnew=dVnew/Ac;

x2=zeros(2,nc);     % Position Matrices: Row 1 gives x-coordinate
x2start=zeros(2,1);

% x-coordinates and their corresponding element no. of all points: 
% Variable volume elements up to element preceeding occupied
for n=2:kpold-1                
    x2(1,n)=dxVar*(n-1);
    x2(2,n)=ceil(x2(1,n)/dxc);
end

% New piston position (and end of Vend volume)
xp=dxVar*(kpold-2)+dx2end;
if kpold==1
    xp=0;
end

x2start(1)=xp+dx2start;     % Vstart volume end position
x2start(2)=ceil(x2start(1)/dxc);

% Variable volume elements up to cylinder end
for n=kpold+1:nc             
    x2(1,n)=x2start(1)+dxVar*(n-kpold);
    x2(2,n)=ceil(x2(1,n)/dxc);
end

% Prepare for mixing procedure
% Transfer all temperatures to vector Temps
Temps(2:kpold-1)=T2(2:kpold-1);
Temps(kpold)=T2e;
Temps(kpold+1)=T2s;
Temps(kpold+2:nc+1)=T2(kpold+1:nc);
Temps(nc+2)=Tout;

% Transfer all positions and elemenet no. to matrix Pos
Pos(:,2:kpold-1)=x2(:,2:kpold-1);
Pos(1,kpold)=xp;
Pos(2,kpold)=ceil(xp/dxc);
Pos(:,kpold+1)=x2start;
Pos(:,kpold+2:nc+1)=x2(:,kpold+1:nc);
Pos(1,nc+2)=Pos(1,nc+1)+dxnew;
Pos(2,nc+2)=nc+1;

%% %%  Mix variable volumes back into initial spatial volumes  %% %%

kpnew=ceil(xp/dxc); % Indicates newly occupied element
if kpnew==0
    kpnew=1;
end

% For each cylinder element
for j=1:nc
    n=2;
    % Incrememnt n to first gas element occupying j
    while Pos(2,n)<j
        n=n+1;
    end
    % If n is only element in j - Update temperature
    if Pos(2,n)>j
        TC(j)=Temps(n);
    % Otherwise sum all contriutions of gas to element j
    else
        TC(j)=(Pos(1,n)-dxc*(j-1))/Temps(n);
        n=n+1;
        while Pos(2,n)<j+1
            TC(j)=TC(j)+((Pos(1,n)-Pos(1,n-1))/Temps(n));
            n=n+1;
        end
    % Finish by mixing contributions
        TC(j)=dxc/(TC(j)+(dxc*j-Pos(1,n-1))/Temps(n));
    end
end

% Special case for partial element preceding piston
j=kpnew;
if kpold>1
    n=2;
    while Pos(2,n)<j
        n=n+1;
    end
    if n==kpold
        Te=Temps(n);
    else
        Te=(Pos(1,n)-dxc*(j-1))/Temps(n);
        for n=n+1:kpold
            Te=Te+((Pos(1,n)-Pos(1,n-1))/Temps(n));
        end
        Te=(xp-dxc*(j-1))/Te;
    end
end

% Special case for partial element following piston
n=kpold+1;
if Pos(2,n)>j;
    Ts=Temps(n);
else
    Ts=(Pos(1,n)-xp)/Temps(n);
    n=n+1;
    while Pos(2,n)<j+1
        Ts=Ts+((Pos(1,n)-Pos(1,n-1))/Temps(n));
        n=n+1;
    end
    Ts=(dxc*j-xp)/(Ts+(dxc*j-Pos(1,n-1))/Temps(n));
end


%% %%                   Update partial volumes                 %% %%

Ve=(xp-(kpnew-1)*dxc)*Ac;
Vs=(kpnew*dxc-xp)*Ac;


%% %%    Check for mass and internal energy conservation       %% %%

% Find system mass
mtot=(P2*VelemC/R)*sum(1./TC(:,1:kpnew-1));
mtot=mtot+((P2*Ve)/(R*Te))+((P2*Vs)/(R*Ts));
mtot=mtot+(P2*VelemC/R)*sum(1./TC(:,kpnew+1:nc));

% Find change in internal energy
Uf=(P2*cvgas/R)*((nc-1)*VelemC+Ve+Vs);
dU=Uf-Ui;
dQR=sum(dQR);

end

