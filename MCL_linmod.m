%% MCL_init
% This script initialises the MCL model from Simulink
% and transfers it into the state-space linearised model  
% and into the transfer function matrix 

clear all;
clc;

%% init all parameters
MCL_init;

%% Trim model

% Calculate operation point specifications
op = operspec('MCL');

% Limit input
op.inputs(1).Description='U_GP1';
set(op.inputs(1), 'Min', 0);
set(op.inputs(1), 'Max', 24);
op.inputs(1).u=10;
op.inputs(2).Description='U_VCA1';
set(op.inputs(2), 'Min', 0);
set(op.inputs(2), 'Max', 0);
op.inputs(2).u=7;
op.inputs(3).Description='U_GP12';
set(op.inputs(3), 'Min', 0);
set(op.inputs(3), 'Max', 24);
op.inputs(3).u=5;
op.inputs(4).Description='U_GP2';
set(op.inputs(4), 'Min', 0);
set(op.inputs(4), 'Max', 24);
op.inputs(4).u=2;
op.inputs(5).Description='U_VCA2';
set(op.inputs(5), 'Min', -24);
set(op.inputs(5), 'Max', 0);
op.inputs(5).u=1;

% Limit outputs/set output conditions
op.outputs(1).Description='x_VCA1';
set(op.outputs(1), 'Min', -0.0127);   % limit x_VCA[m]
set(op.outputs(1), 'Max',  10);   % limit x_VCA[m]
op.outputs(2).Description='P_pc1';
set(op.outputs(2), 'Min', 40);       % limit P_pc1 [mmHG]
set(op.outputs(2), 'Max', 200);       % limit P_pc1 [mmHG]

op.outputs(3).Description='x_VCA2';
set(op.outputs(3), 'Min', -0.0127);   % limit x_VCA2[m]
set(op.outputs(3), 'Max',  0.0127);   % limit x_VCA2[m]
op.outputs(4).Description='P_pc2';
set(op.outputs(4), 'Min', -20);       % limit P_pc2 [mmHG]
set(op.outputs(4), 'Max', 100);       % limit P_pc2 [mmHG]

op.outputs(5).Description='F_GP1';
set(op.outputs(5), 'Min', 0);        % limit F_GP1 [mmHG]
set(op.outputs(5), 'Max', 100);        % limit F_GP1 [mmHG]
 op.outputs(5).y=3;                   % F_GP1 [L/min]

% Speed of rotation
op.states(4).Known = 1; %Duke
op.states(2).Known = 1; %Duke
op.states(13).Known = 1; %Duke
op.states(4).SteadyState = 1; %Duke
op.states(2).SteadyState = 1; %Duke
op.states(13).SteadyState = 1; %Duke
op.states(4).x=700*(1/RADS2TURNMIN);  % w_GP1 %Dipankar
op.states(2).x=500*(1/RADS2TURNMIN);  % w_GP12 %Dipankar
op.states(13).x=500*(1/RADS2TURNMIN); % w_GP2 %Dipankar

% Limit states/set state conditions
set(op.states(1), 'Min', -7.5);   % Limit i_GP12 [A]
set(op.states(1), 'Max', 7.5);    % Limit i_GP12 [A]
set(op.states(3), 'Min', -7.5);   % Limit i_GP1 [A]
set(op.states(3), 'Max', 7.5);    % Limit i_GP1 [A]
set(op.states(12), 'Min', -7.5);  % Limit i_GP2 [A]
set(op.states(12), 'Max', 7.5);   % Limit i_GP2 [A]
set(op.states(11), 'Min', -7.5);  % Limit i_VCA1 [A]
set(op.states(11), 'Max', 7.5);   % Limit i_VCA1[A]
set(op.states(20), 'Min', -7.5);  % Limit i_VCA2 [A]
set(op.states(20), 'Max', 7.5);   % Limit i_VCA2[A]
set(op.states(5), 'Min', -20);    % Limit P_pc1 [mmHG]?
set(op.states(5), 'Max', 200);    % Limit P_pc1 [mmHG]?
set(op.states(14), 'Min', -20);   % Limit P_pc2 [mmHG]?
set(op.states(14), 'Max', 200);   % Limit P_pc2 [mmHG]?
set(op.states(10), 'Min', 0);     % limit x_VCA1 [m]
set(op.states(10), 'Max', 0.0137) % limit x_VCA1 [m]
set(op.states(19), 'Min', 0);     % limit x_VCA2 [m]
set(op.states(19), 'Max', 0.0137) % limit x_VCA2 [m]

%opt = findopOptions; %Duke
opt = linoptions; %Dipankar
opt.DisplayReport='on';
%opt.OptimizerType='graddescent-elim'; %Duke
opt.OptimizerType='lsqnonlin'; %Dipankar
opt.OptimizationOptions.Algorithm= 'trust-region-reflective'; %Dipankar
opt.OptimizationOptions.DiffMaxChange = 0.1;
opt.OptimizationOptions.MaxIter = 5000;
opt.OptimizationOptions.MaxFunEvals = 1000;
opt.OptimizationOptions.TolFun = 1.0e-004;
opt.OptimizationOptions.TolX = 1.0e-004;

%% Trim model
[op_point, op_report] = findop('MCL', op, opt);

% Get trim results
U  = [op_point.inputs(1).u; op_point.inputs(2).u; op_point.inputs(3).u; op_point.inputs(4).u; op_point.inputs(5).u];
DX = [op_report.states(1).dx; op_report.states(2).dx; op_report.states(3).dx; op_report.states(4).dx; op_report.states(5).dx; op_report.states(6).dx; op_report.states(7).dx; op_report.states(8).dx; op_report.states(9).dx; op_report.states(10).dx; op_report.states(11).dx; op_report.states(12).dx; op_report.states(13).dx; op_report.states(14).dx; op_report.states(15).dx; op_report.states(16).dx; op_report.states(17).dx; op_report.states(18).dx; op_report.states(19).dx; op_report.states(20).dx ];
X  = [op_report.states(1).x; op_report.states(2).x; op_report.states(3).x; op_report.states(4).x; op_report.states(5).x; op_report.states(6).x; op_report.states(7).x; op_report.states(8).x; op_report.states(9).x; op_report.states(10).x; op_report.states(11).x; op_report.states(12).x; op_report.states(13).x; op_report.states(14).x; op_report.states(15).x; op_report.states(16).x; op_report.states(17).x; op_report.states(18).x; op_report.states(19).x; op_report.states(20).x];
Y  = [op_report.outputs(1).y; op_report.outputs(2).y; op_report.outputs(3).y; op_report.outputs(4).y; op_report.outputs(5).y];

%% init the MCL model from Simulink

% Load state-space data from MCL simulink model
[Al, Bl, Cl, Dl]=linmod('MCL', X, U);

% State-space linearised model
Sys=(ss(Al, Bl, Cl, Dl));
set(Sys,'inputname',[{'U\_GP1'} {'U\_VCA1'} {'U\_GP12'} {'U\_GP2'} {'U\_VCA2'}], 'outputname',[{'x\_VCA1'} {'P\_pc1'} {'x\_VCA2'} {'P\_pc2'} {'F\_GP1'}]);

%% Task 3.1.1. The stable operating point
fprintf('\nTask 3.1.1. The stable operating point');
format shortG; %Duke
display(U);
display(Y)
%display(X([4, 2, 13],1)*RADS2TURNMIN); %Duke

%% Task 3.1.2. Analysis of the linearised model
fprintf('Task 3.1.2. Analysis of the linearised model\n');
fprintf('The dimension of the system matrix is: %g.\n', length(Al));
[b,a] = ss2tf(Al,Bl,Cl,Dl,5); %Have not remove 0 elements of b
fprintf('The order of the tf matrix of the system is: %g.\n', max(max(size(a)),max(size(b)))-1);

polesSys = pole(Sys);
fprintf('Number of unstable poles: %g.\n', size(polesSys(real(polesSys)>0,:),1));

zerosSys = tzero(Sys);
fprintf('Number of RHP zeros: %g.\n\n', size(zerosSys(real(zerosSys)>0,:),1));

%% Task 3.1.3. Controllability and observability
fprintf('Task 3.1.3. Controllability and observability\n');
tol = 5E-08;
% Controllability
[Actrb,Bctrb,Cctrb,Tctrb,kctrb] = ctrbf(Al,Bl,Cl,tol);
% Observability
[Aobsv,Bobsv,Cobsv,Tobsv,kobsv] = obsvf(Al,Bl,Cl,tol);

fprintf('The number of states are not controllable: %g.\n', length(Al)-sum(kctrb));
fprintf('The number of states are not observable: %g.\n\n', length(Al)-sum(kobsv));

%% Task 3.2.1. Hankel Singular Values
fprintf('Task 3.2.1. Hankel Singular Values\n\n');
hsv = hsvd(Sys);
figure(1);
h = hsvplot(Sys);

%% Task 3.2.2. Reduce the model
fprintf('Task 3.2.2. Reduce the model\n');
[redSys3] = balred(Sys,3);
[redSys5] = balred(Sys,5);
[redSys10] = balred(Sys,10);
[redSys11] = balred(Sys,11);
[redSys12] = balred(Sys,12);

t = 0:0.1:10;
s0 = step(Sys,t);
s3 = step(redSys3,t);
s10 = step(redSys10,t);
s11 = step(redSys11,t);
s12 = step(redSys12,t);

figure(2);
i=4; j=4;
%plot(s0(:,i,j),'DisplayName','original system');
plot(s0(:,i,j)-s3(:,i,j),'DisplayName','redSys3'); hold on;
plot(s0(:,i,j)-s10(:,i,j),'DisplayName','redSys10');
plot(s0(:,i,j)-s11(:,i,j),'DisplayName','redSys11');
plot(s0(:,i,j)-s12(:,i,j),'DisplayName','redSys12');
legend;
xlabel('time'); ylabel('response');

[b_red,a_red] = ss2tf(redSys12.a,redSys12.b,redSys12.c,redSys12.d,5);
fprintf('The order of the reduced model: %g.\n\n', max(max(size(a_red)),max(size(b_red)))-1);

%% Task 3.2.3. Minimal realisation
fprintf('Task 3.2.3. Minimal realisation');
minSys = minreal(redSys3, tol);
display(minSys);

%% Task 3.2.4. RGA
fprintf('Task 3.2.4. Static Relative Gain Array');
minG = zeros(5);
for i=1:5
    [minA, minB] = ss2tf(minSys.a, minSys.b, minSys.c, minSys.d, i);
    minG(:,i) = minA(:,1) / minB(1);
end
sRGA = minG.*pinv(minG).';
rsRGA = round(sRGA,2);
display(rsRGA);

%% Task 3.2.5. Scaling
fprintf('Task 3.2.5. Scaling');
Du = diag([24, 24, 24, 24, 24]);
De = eye(5);
Dd = eye(5);

scaled_G = pinv(De) * minG * Du;
display(scaled_G);

%% Task 3.3.1. Sources for uncertainty
fprintf('Task 3.3.1. Sources for uncertainty\n\n');

%% Task 3.3.2. Weighting functions for uncertainty
fprintf('Task 3.3.2. Weighting functions for uncertainty\n\n');

%% Task 3.3.3. Visualize uncertain system
fprintf('Task 3.3.3. Visualize uncertain system\n\n');

%% Task 3.4. H-inf-controller synthesis
fprintf('Task 3.4. H-inf-controller synthesis\n');

%%
fprintf('\nThank you for reviewing this report!\n');
