function [qideal,qactual,qratio,mdot,a] = Comp(Pratio,Pi,UAhe,...
                                               UAreg,nreg,Lcyl,rcyl)
% Comp.m: Finds the parameters inside a thermal compresser as the 
% piston moves back and forth along a cylinder

% Model based on taking a constant volume of gas out of the end of
% the cylinder with each step


%% %%                     Set-up Parameters                    %% %%

    % Piston
    global f
    f=10;               % Piston frequency                     / Hz
    nrev=10;            % No. revolutions in simulation

    % Cylinder
    global Lc nc
    Lc=Lcyl;            % Cylinder length                      / m
    rc=rcyl;            % Cylinder radius                      / m
    nc=250;             % Number of elements in cylinder

    % Regenerator
    global nr UAr UAe
    rr=0.02;            % Regenerator radius                   / m
    nr=nreg;            % Number of elements in regenerator
    dxr=0.001;          % Regenerator element length           / m
    p=0.99;             % Volume fraction of metal in the reg.
    UAr=UAreg;          % Reg Thermal 'conductance'            / W/K
    UAe=UAhe;           % HE Thermal 'conductance'             / W/K

    % Temperatures
    global Th Tc Tamb
    Tamb=293;           % Temperature of cooler (ambient)      / K
    Th=873;             % Temperature of heater                / K
    Tc=268;             % Temp of gas from evaporator          / K

    % Pressures
    Pinit=2000000;      % Initial system pressure              / Pa
    Pchki=Pi;           % Inlet check valve pressure           / Pa
    Pchko=Pchki*Pratio; % Outlet check valve pressure          / Pa


%% %%                    Physical Constants                    %% %%

    % Cylinder
    global Ac dxc VelemC Vtot
    Ac=pi*rc*rc;                % Cylinder c-s area         / m^2
    dxc=Lc/nc;                  % Cylinder element length   / m
    VelemC=Ac*dxc;              % Cylinder element volume   / m^3
    Vtot=(Ac*Lc);               % Cylinder total volume     / m^3

    % Regenerator
    Ar=pi*rr*rr;                % Regenerator c-s area      / m^2
    VelemR=Ar*dxr;              % Regenerator element volume/ m^3

    % Piston
    global om
    om=2*pi*f;                  % Piston angular frequency  / rads/s
    dtmax1=acos((Lc-8*dxc)/Lc)/om; % Max step through reg - s1   / s
    dtmax3=acos((Lc-2*dxc)/Lc)/om; % Max step through reg - s3   / s


%% %%             Constant thermodynamic properties            %% %%

    % System wide
    global R y cvgas cpgas

    R=8314.5/44;   % Gas constant for CO2                   / J/kgK
    cpgas=844;     % Specific heat at constant pressure     / J/kgK
    cvgas=cpgas-R; % Specific heat at constant volume       / J/kgK
    y=cpgas/cvgas; % Gamma - ratio of specific heats (cp/cv)

    % Regenerator
    global cpmet dmmet
    cpmet=106;          % Specific heat capacity of metal   / J/kgK 
    rho=7870;           % Density of metal                  / kg/m3
    dmmet=p*VelemR*rho; % Metal mass in regenerator element / kg


%% %%                   Initial Conditions                     %% %%

    % System wide
    t(1)=0;             % Initial time                      / s
    xp(1)=0;            % Initial piston position           / m
    P(1)=Pinit;         % Initial pressure of the system    / Pa

    % Masses
    mtot(1)=(P(1)*Vtot)/(R*Tamb);   % Initial system mass   / kg
    dmcd(1)=0;              % Mass of gas to condenser      / kg
    dmev(1)=0;              % Mass of gas from evaporator   / kg

    % Energies
    dQin(1)=0;              % Heat in with each step        / J
    dQout(1)=0;
    dQR(1)=0;
    dU(1)=0;                % Change in internal energy     / J
    Ui=mtot(1)*cvgas*Tamb;  % Initial system internal energy/ J

    % Temperatures        
    TC=Tamb*ones(1,nc);     % Cylinder temperatures         / K
    TRmet=Tamb*ones(1,nr);  % Regenerator metal temps       / K
    TRgas=zeros(1,nr);      % Regenerator gas temps         / K
    Tcd(1)=0;

    % Partial volumes and temperatures for occupied element
    Vs=VelemC;  % Partial element volume after piston       / m^3
    Ts=Tamb;    % Temperature in the start variable element / K
    Ve=0;       % Partial element volume before piston      / m^3
    Te=Tamb;    % Temperature in the end variable element   / K
    
    % Initialise necessary variables
    mh=0;       % Average maximum mass
    ml=0;       % Average minimum mass
    qdot=0;     % Heat in per unit time
    n=0;        % Step counter
    a=0;        % Condition for whether cycle is stable
                % a=0 - Both check valves remain closed
                % a=1 - Only outlet check valve opens
                % a=2 - Only inlet check valve opens
                % a=3 - Both open, system working

    i=1;        % Inititalise step index to 1
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%                    Start simulation                      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    global rev
    for rev=1:nrev
    
        % Display current revolution number
        disp('revolution')
        disp(rev)
        
    
%% %% Stage 1: Check valves closed - Forward half of revolution%% %%
    
        % Run stage 1 while conditions are satisfied
        while P(i)<Pchko && xp(i)<Lc
            [TC(i+1,:),Ts,Te,Vs,Ve,TRmet(i+1,:),P(i+1),xp(i+1),...
                t(i+1),mtot(i+1),dQin(i+1),dQout(i+1),dQR(i+1),...
                dU(i+1),Ui] = Stage1(t(i),dtmax1,xp(i),P(i),...
                TC(i,:),Ts,Te,Vs,Ve,TRmet(i,:),Ui);
            % Increment step index
            i=i+1;
            % If in final revolution perform stability check
            if rev==nrev
                mh=mh+mtot(i);
                qdot=qdot+dQin(i);
                n=n+1;
            end
        end
        if rev==nrev
            mh=mh/n;
            n=0;
        end
        
    
%% %% Stage 2: Outlet valve open - Forward half of revolution  %% %%
    
        % Run stage 2 while conditions are satisfied
        while xp(i)<Lc
            [TC(i+1,:),Ts,Te,Vs,Ve,TRmet(i+1,:),P(i+1),xp(i+1),...
                t(i+1),dmcd(i+1),Tcd(i+1),mtot(i+1),dQin(i+1),...
                dQout(i+1),dQR(i+1),dU(i+1),Ui] = Stage2(t(i),...
                xp(i),P(i),TC(i,:),Ts,Te,Vs,Ve,TRmet(i,:),Ui);
            % Increment step index
            i=i+1;
            % If in final revolution perform stability check
            if rev==nrev
                a=1;
                qdot=qdot+dQin(i);
            end
        end
        

%% %% Stage 3: Check valves closed - Reverse half of revolution%% %%    
    
        % Run stage 3 while conditions are satisfied
        while P(i)>Pchki && xp(i)>0
            [TC(i+1,:),Ts,Te,Vs,Ve,TRmet(i+1,:),P(i+1),xp(i+1),...
                t(i+1),mtot(i+1),dQin(i+1),dQout(i+1),dQR(i+1),...
                dU(i+1),Ui] = Stage3(t(i),dtmax3,xp(i),P(i),...
                TC(i,:),Ts,Te,Vs,Ve,TRmet(i,:),Ui);
            % Increment step index
            i=i+1;
            % If in final revolution perform stability check
            if rev==nrev
                ml=ml+mtot(i);
                qdot=qdot+dQin(i);
                n=n+1;
            end
        end
        if rev==nrev
            ml=ml/n;
            mdot=mh-ml;
        end
        

%% %% Stage 4: Inlet valve open - Reverse half of revolution   %% %%    
    
        % Run stage 4 while conditions are satisfied
        while xp(i)>0
            [TC(i+1,:),Ts,Te,Vs,Ve,TRmet(i+1,:),P(i+1),xp(i+1),...
                t(i+1),dmev(i),mtot(i+1),dQin(i+1),dQout(i+1),...
                dQR(i+1),dU(i+1),Ui] = Stage4(t(i),xp(i),P(i),...
                TC(i,:),Ts,Te,Vs,Ve,TRmet(i,:),Ui);
            % Increment step index
            i=i+1;
            % If in final revolution perform stability check
            if rev==nrev
                if a==0 || a==2
                    a=2;
                elseif a==1 || a==3
                    a=3;
                end
                qdot=qdot+dQin(i);
            end
        end
        

    end
    
    
%% %%               Post Simulation Processing                 %% %%

    % Determine the ideal heat for the analysis
    qideal=cpgas*Tc*((Pchko/Pchki)^((y-1)/y)-1)*(Th/(Th-Tamb))/1000;
    % Determine the actual heat for the analysis
    qactual=(qdot/mdot)/1000;
    % Ratio of actual heat to ideal heat - efficiency measure
    qratio=qactual/qideal;
    
    
end

