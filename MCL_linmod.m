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
op.states(4).Known = 1; 
op.states(2).Known = 1; 
op.states(13).Known = 1;
op.states(4).x=700*(1/(RADS2TURNMIN));  % w_GP1
op.states(2).x=500*(1/(RADS2TURNMIN));  % w_GP12
op.states(13).x=500*(1/(RADS2TURNMIN)); % w_GP2


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

%opt = findopOptions('DisplayReport','on','OptimizerType','graddescent');
opt = linoptions;
opt.DisplayReport='on';
opt.OptimizerType='lsqnonlin';
opt.OptimizationOptions.Algorithm= 'trust-region-reflective';
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

format shortG;
display(X([4, 2, 13],1)*RADS2TURNMIN);

%% init the MCL model from Simulink
op_point_new = load('MCL_op.mat');
X = cell2mat(get(op_point_new.op.States,'x'));
U = cell2mat(get(op_point_new.op.Input,'u'));

% Load state-space data from MCL simulink model
[Al, Bl, Cl, Dl]=linmod('MCL', X, U);

% State-space linearised model
Sys=(ss(Al, Bl, Cl, Dl));
set(Sys,'inputname',[{'U\_GP1'} {'U\_VCA1'} {'U\_GP12'} {'U\_GP2'} {'U\_VCA2'}], 'outputname',[{'x\_VCA1'} {'P\_pc1'} {'x\_VCA2'} {'P\_pc2'} {'F\_GP1'}]);

%%%%%%%%%%%%%%%%%%%


%% Task 3.1.1. The stable operating point
fprintf('\nTask 3.1.1. The stable operating point');

format shortG;
display(U);
display(Y);
display(DX);
if sum(DX) ~= 0
    fprintf('We do not have a STABLE operating point (DX not equal 0).\n');
end


%% Task 3.1.2 Analysis of the Linearised Model
fprintf('\nTask 3.1.2. Analysis of the linearised model\n\n');
fprintf('The dimension of the tranfer function matrix matrix is: %g.\n', length(Al));

TFSys = tf(Sys);
OSys=order(Sys);

%[b,a] = ss2tf(Al,Bl,Cl,Dl,5); %Tranfer Function of the System
fprintf('The order of the System is: %g.\n',OSys);

P = pole(Sys); %Poles of the Tranfer Function
fprintf('Number of unstable poles of the System is: %g. \n', size(P(real(P)>0,:),1));
Z = tzero(Sys);
fprintf('Number of RHP zeros of the system is: %g. \n\n', size(Z(real(Z)>0,:),1));


%% Task 3.1.3 Controllabililty & Observability
fprintf('Task 3.1.3. Controllability & Observability\n\n');
tol = 5E-8;

%Controlabillity Check
[AbarC,BbarC,CbarC,TC,kC] = ctrbf(Al,Bl,Cl,tol);

fprintf("\nNumber of Un-Controllable states in the system are: %d", length(Al)-sum(kC));

%Observability Check
[AbarO,BbarO,CbarO,TO,kO] = obsvf(Al,Bl,Cl,tol);

fprintf("\nNumber of Un-Observable states in the system are: %d\n\n", length(Al)-sum(kO));


%% Task 3.2.1 Hankel Singular Values of the Dynamic System
fprintf('Task 3.2.1. Hankel Singular Values of the System\n\n');

hsv = hsvd(Sys);
figure(1);
h = hsvplot(Sys);

%% Task 3.2.2 Model Reduction
fprintf('Task 3.2.2. Model Reduction\n\n');

[SysS, SysU] = stabsep(ssbal(Sys));
hsvs = hsvd(SysS);

figure(2);
hsvplot(SysS);

% Reduce the stable model part
redOrder = 5;
[redSysS] = balred(SysS,redOrder);

fprintf('The order of the reduced stable model : %g.\n\n', order(redSysS));
figure(3);
step(SysS,'r',redSysS,'g--')
legend;
xlabel('time'); ylabel('response');

% Combine the reduced stable part and the unstable part
orderRedSysU = length(SysU.A);
redA = [redSysS.A,                    zeros(redOrder,orderRedSysU);
        zeros(orderRedSysU,redOrder), SysU.A];
redB = [redSysS.B; SysU.B];
redC = [redSysS.C, SysU.C];
redD = redSysS.D + SysU.D;
redSys = ss(redA, redB, redC, redD);
figure(4);
step(Sys,'r',redSys,'g--')
legend;
xlabel('time'); ylabel('response');

%% Task 3.2.3 Minimal Realisation
fprintf('Task 3.2.3. Minimal Realisation of the System\n\n');

minSys = minreal(redSysS, tol);
display(minSys); %minimal realization of the reduced model 


%% Task 3.2.4 RGA
fprintf('\nTask 3.1.2. RGA evaluation\n');

G = dcgain(minSys);
ARga = (G).*pinv(G).'; %RGA computation from SS model
display(round(ARga,2));

TFrSysr = tf(minSys);
%G1 = dcgain(TFrSysr);
%ARga1 = (G1).*pinv(G1).'; %RGA computation from TF of SS model
%display(ARga1);

%% Task 3.2.5 Scaling
fprintf('\nTask 3.2.5. Scaling\n');
Du = diag([24,24,24,24,24]);
Dy = eye(5);
Dd = eye(5);

TFrSysr_scaled = pinv(Dy)*TFrSysr*Du;
%display(TFrSysr_scaled);


%% Task 3.3.2. Weighting functions for uncertainty
fprintf('Task 3.3.2. Weighting functions for uncertainty\n\n');

%G = minreal(ss(TFrSysr_scaled)).NominalValue;
%Wxvca12 = makeweight(.1,5,1);
%Wppc12 = makeweight(.3,2.5,1);
%wppfgpq = makeweight(.1,1.5,1);
%Delta1 = ultidyn('Delta1',[1 1]);
%Delta2 = ultidyn('Delta2',[1 1]);
%Delta3 = ultidyn('Delta3',[1 1]);
%G = H*blkdiag(1+W1*Delta1,1+W2*Delta2,1+W1*Delta1,1+W2*Delta2, 1+w1*Delta3);

Gs = minreal(ss(TFrSysr_scaled));
G = uss(Gs);
[Li,phi,w0] = bode(Gs-G.NominalValue/G.NominalValue);
