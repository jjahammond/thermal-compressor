% CompAnalysis.m sets up and executes analyses for varying input
% conditions to the compressor model

% The output is a measure of performance variable - qratio
% This gives the ratio of the actual heat required per unit mass
% of gas passing through the system to the ideal heat required
% if the compression were isentropic. This indicator tends to 
% infinity if the system stops producing the require pressure
% ratio as the mass passing through the system with each revolution
% falls to zero.


%% %% Choose parameters to analyse and set others as constant  %% %%

%Pratio=1.6;    % Pressure Ratio
Pi=3000000;     % Inlet check valve openning pressure   / Pa
UAhe=500;       % Thermal 'conductance' of exchanger    / W/K
UAreg=500;      % Thermal 'conductance' of regenerator  / W/K
%nreg=150;      % No. of regenerator elements
Lcyl=0.275;     % Cylinder length                       / m
rcyl=0.035;     % Cylinder radius                       / m


for i=1:16
    % Set relationship for parameter A
    nreg(i)=25*i;
    for j=1:40
        % Set relationship for parameter B
        Pratio(j)=1+0.05*j;
        % Display current analysis
        disp(i);
        disp(j);
        % Execute Comp.m model for current analysis
        [qideal(i,j),qactual(i,j),qratio(i,j),mdot(i,j),...
         a(i,j)] = Comp(Pratio(j),Pi,UAhe,UAreg,nreg(i),Lcyl,rcyl);
        % Clear display
        clc
    end
end