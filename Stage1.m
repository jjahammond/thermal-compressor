function [TC,Ts,Te,Vs,Ve,Tmet,P2,xp,t2,m,dQi,dQo,dQR,...
                  dU,Uf]=Stage1(t,dtmax,xp,P,TC,Ts,Te,Vs,Ve,Tmet,Ui)
% Stage1.m updates all temperatures and pressures inside the
% cylinder for the first stage of the cycle

%   Stage 1 - both check valves closed, piston in forwards half 
%   of its revolution

global om f rev
global Vtot VelemC Lc Ac dxc nc nr
global Th Tamb R y cvgas cpgas UAr UAe dmmet cpmet

%% %%             Find new pressure of the system              %% %%

% kpold indicates original element occupied by piston
kpold=ceil(xp/dxc);
if kpold==0
    kpold=1;
end

% Find mass of gas in final element and new pressure when passed 
% through the heater. Use full element volume unless piston occupies
% final element - in this case use partial volume after the piston
if kpold<nc
    dmgas=(P*VelemC)/(R*TC(nc));    % Mass of heated gas
    
%%%% Method of Bisections to find temperature out of regenerator
    dtmin=0;
    xp1=0;
    xp2=1;
    % Iterate until displacements converge
    while abs(xp2-xp1)>1e-10
        % Set initial time interval to analyse
        dt=(dtmin+dtmax)/2;
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
        
        % Determine new pressure
        P2=fzero(@(P2) Vtot-(Vtot-VelemC)*(P2/P)^(-1/y)-...
                                (dmgas*R*Totemp/P2),P);
        Vtemp=dmgas*R*Totemp/P2;  
        
        % Displacement resulting from compression
        xp1=(Vtemp+((P2/P)^(-1/y))*((kpold-1)*VelemC+Ve))/Ac;

        % Displacement resulting from time interval alone
        t2=t+dt;
        xp2=Lc*(1-cos(om*t2))/2;

        % Conditions to update interval analysed
        if xp1>xp2
            dtmin=dt;
        elseif xp1<xp2
            dtmax=dt;
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
    
else

%%%% Case where piston occupies final element
    dmgas=(P*Vs)/(R*Ts);
    
    % Find time interval
    t2=(2*rev-1)/(2*f);
    dt=t2-t;
    
    % Run through heater and regenerator directly
    Tout=Tamb+(Ts-Tamb)*exp((-UAe*dt)/(dmgas*cpgas));
    dQo=dmgas*cpgas*(Tout-Ts);
    for i=1:nr
        [Tmet(i),Tout,dQR(i)]=RegElem(Tout,Tmet(i),dmgas,dt);
    end
    Tin=Tout;
    Tout=Th-(Th-Tin)*exp((-UAe*dt)/(dmgas*cpgas));
    dQi=dmgas*cpgas*(Tout-Tin);

    % Determine new system pressure
    P2=fzero(@(P2) Vtot-(Vtot-Vs)*(P2/P)^(-1/y)-...
                            (dmgas*R*Tout/P2),P);
end


%% %%            Update temperatures and volumes               %% %%

% Temperatures of variable volumes  / K
T2(1:kpold-1)=TC(1:kpold-1)*((P2/P)^((y-1)/y));
T2(kpold+1:nc-1)=TC(kpold+1:nc-1)*((P2/P)^((y-1)/y));
T2e=Te*((P2/P)^((y-1)/y));      % Partial element before piston
T2s=Ts*((P2/P)^((y-1)/y));      % Partial element after piston

% Variable volumes  / m^3
V2=VelemC*(P2/P)^(-1/y);        % General variable volume
V2e=Ve*(P2/P)^(-1/y);           % Partial volume before piston
V2s=Vs*(P2/P)^(-1/y);           % Partial volume after piston
dVnew=dmgas*R*Tout/P2;          % New volume of heated gas


%% %%      Find position of all elements within cylinder       %% %%

% Variable volume lengths
dxVar=V2/Ac;        % General variable volume length
dx2e=V2e/Ac;        % Partial volume before piston length
dx2s=V2s/Ac;        % Partial volume after piston length

xh=zeros(2,1);      % Position Matrices: Row 1 gives x-coordinate
x2=zeros(2,nc-1);   %                    Row 2 gives element no.
x2s=zeros(2,1);

% x-coordinates and their corresponding element no. of all points: 
xh(1)=dVnew/Ac;     % Heated element end point
xh(2)=ceil(xh(1)/dxc);

% Variable volume elements up to element preceeding occupied
for n=1:kpold-1        
    x2(1,n)=xh(1)+dxVar*n;
    x2(2,n)=ceil(x2(1,n)/dxc);
end

% New piston position (and end of Ve volume)
xp=xh(1)+dxVar*(kpold-1)+dx2e;
if kpold==nc
    xp=Lc;
end

x2s(1)=xp+dx2s;     % Vstart volume end position
x2s(2)=ceil(x2s(1)/dxc);

% Variable volume elements up to cylinder end
for n=kpold+1:nc-1                  
    x2(1,n)=xh(1,1)+dx2e+dx2s+dxVar*(n-1);
    x2(2,n)=ceil(x2(1,n)/dxc);
end

% Prepare for mixing procedure
% Transfer all temperatures to vector Temps
Temps(1)=Tout;
Temps(2:kpold)=T2(1:kpold-1);
Temps(kpold+1)=T2e;
Temps(kpold+2)=T2s;
Temps(kpold+3:nc+1)=T2(kpold+1:nc-1);

% Transfer positions and element no. to matrix Pos
Pos(:,1)=xh;
Pos(:,2:kpold)=x2(:,1:kpold-1);
Pos(1,kpold+1)=xp;
Pos(2,kpold+1)=ceil(xp/dxc);
Pos(:,kpold+2)=x2s;
Pos(:,kpold+3:nc+1)=x2(:,kpold+1:nc-1);
Pos(2,nc+1)=nc+1;


%% %%  Mix variable volumes back into initial spatial volumes  %% %%

kpnew=ceil(xp/dxc); % Indicates newly occupied element
if kpnew==nc+1
    kpnew=nc;
end

% For each cylinder element
for j=1:nc
    n=1;
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
            TC(j)=TC(j)+(Pos(1,n)-Pos(1,n-1))/Temps(n);
            n=n+1;
        end
    % Finish by mixing contributions
        TC(j)=dxc/(TC(j)+(dxc*j-Pos(1,n-1))/Temps(n));
    end
end

% Special case for partial element preceeding piston 
j=kpnew;
n=1;
while Pos(2,n)<j
    n=n+1;
end
if n==kpold+1
    Te=Temps(n);
else
    Te=(Pos(1,n)-dxc*(j-1))/Temps(n);
    for n=n+1:kpold+1
        Te=Te+(Pos(1,n)-Pos(1,n-1))/Temps(n);
    end
    Te=(xp-dxc*(j-1))/Te;
end

% Special case for partial element following piston
if kpold<nc
    n=kpold+2;
    if Pos(2,n)>j
        Ts=Temps(n);
    else
        Ts=(Pos(1,n)-xp)/Temps(n);
        n=n+1;
        while Pos(2,n)<j+1
            Ts=Ts+(Pos(1,n)-Pos(1,n-1))/Temps(n);
            n=n+1;
        end
        Ts=(dxc*j-xp)/(Ts+(dxc*j-Pos(1,n-1))/Temps(n));
    end
end


%% %%                   Update partial volumes                 %% %%

Ve=(xp-(kpnew-1)*dxc)*Ac;
Vs=(kpnew*dxc-xp)*Ac;


%% %%     Check for mass and internal energy conservation      %% %%

% Find system mass
m=(P2*VelemC/R)*sum(1./TC(:,1:kpnew-1));
m=m+((P2*Ve)/(R*Te))+((P2*Vs)/(R*Ts));
m=m+(P2*VelemC/R)*sum(1./TC(:,kpnew+1:nc));

% Find change in internal energy
Uf=(P2*cvgas/R)*((nc-1)*VelemC+Ve+Vs);
dU=Uf-Ui;
dQR=sum(dQR);

end

