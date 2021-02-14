%% MCL_init
% This script initialises the MCL model from Simulink
% and transfers it into the state-space linearised model  
% and into the transfer function matrix 

% init all parameters
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
set(op.outputs(5), 'Min', 0);         % limit F_GP1 [mmHG]
set(op.outputs(5), 'Max', 100);       % limit F_GP1 [mmHG]
 op.outputs(5).y=3;                   % F_GP1 [L/min]

% Speed of rotation
op.states(4).Known = 1;
op.states(2).Known = 1;
op.states(13).Known = 1;
op.states(4).x=700*(1/RADS2TURNMIN);  % w_GP1
op.states(2).x=500*(1/RADS2TURNMIN);  % w_GP12
op.states(13).x=500*(1/RADS2TURNMIN); % w_GP2

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

opt = findopOptions;
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

%% init the MCL model from Simulink
op_point_new = load('MCL_op.mat');
X = cell2mat(get(op_point_new.op.States,'x'));
U = cell2mat(get(op_point_new.op.Input,'u'));

% Load state-space data from MCL simulink model
[Al, Bl, Cl, Dl]=linmod('MCL', X, U);

% State-space linearised model
Sys=(ss(Al, Bl, Cl, Dl));
set(Sys,'inputname',[{'U\_GP1'} {'U\_VCA1'} {'U\_GP12'} {'U\_GP2'} {'U\_VCA2'}], 'outputname',[{'x\_VCA1'} {'P\_pc1'} {'x\_VCA2'} {'P\_pc2'} {'F\_GP1'}]);

%% Task 3.1.1. The stable operating point
fprintf('\nTask 3.1.1. The stable operating point');
format shortG;
display(U);
display(Y)
display(X([4, 2, 13],1)*RADS2TURNMIN);
display(DX);
if sum(DX) ~= 0
    fprintf('We do not have a STABLE operating point (DX not equal 0).\n');
end

%% Task 3.1.2. Analysis of the linearised model
fprintf('Task 3.1.2. Analysis of the linearised model\n');
fprintf('The dimension of the system matrix is: %g.\n', length(Al));
TFSys = tf(Sys);
fprintf('The order of the tf matrix of the system is: %g.\n', order(Sys));

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
[SysS, SysU] = stabsep(ssbal(Sys));
hsvs = hsvd(SysS);

figure(2);
hsvplot(SysS);

%% Reduce the stable model part
redOrder = 5;
[redSysS] = balred(SysS,redOrder);

fprintf('The order of the reduced stable model : %g.\n\n', order(redSysS));
figure(3);
step(SysS,'r',redSysS,'g--')
legend;
xlabel('time'); ylabel('response');
title('Step response of stable and reduced stable models');

%% Combine the reduced stable part and the unstable part
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
title('Step response of initial and reduced models');
%% Task 3.2.3. Minimal realisation
fprintf('Task 3.2.3. Minimal realisation\n');
minSys = minreal(redSysS, tol);
display(minSys);

%% Task 3.2.4. RGA
fprintf('Task 3.2.4. Static Relative Gain Array');
minsG = dcgain(minSys);
sRGA = minsG.*pinv(minsG).';
rsRGA = round(sRGA,3);
display(rsRGA);

%% Task 3.2.5. Scaling
fprintf('Task 3.2.5. Scaling');
Du = diag([24, 24, 24, 24, 24]);
De = diag([0.0137, 200, 0.0137, 200, 1]);
Dd = eye(5);

minG = tf(minSys);
scaled_G = pinv(De) * minG * Du;
display(scaled_G);

%% Task 3.3.2. Weighting functions for uncertainty
fprintf('Task 3.3.2. Weighting functions for uncertainty\n\n');

Delta = ultidyn('Delta', [5 5], 'Bound', 1);
Wxvca = makeweight(.1, 1, 5);
Wppc  = makeweight(.3, 1, 2.5);
Wfgp1 = makeweight(.1, 1, 1.5);

WA = diag([Wxvca.A,Wppc.A,Wxvca.A,Wppc.A,Wfgp1.A]);
WB = diag([Wxvca.B,Wppc.B,Wxvca.B,Wppc.B,Wfgp1.B]);
WC = diag([Wxvca.C,Wppc.C,Wxvca.C,Wppc.C,Wfgp1.C]);
WD = diag([Wxvca.D,Wppc.D,Wxvca.D,Wppc.D,Wfgp1.D]);
W = ss(WA,WB,WC,WD);

scaledSys = minreal(ss(scaled_G));
scaledSys.InputName = {'U_GP1', 'U_VCA1', 'U_GP12', 'U_GP2', 'U_VCA2'};
scaledSys.OutputName = {'x_VCA1', 'P_PC1', 'x_VCA2', 'P_PC2', 'F_GP1'};

Gp = (eye(5)+W*Delta)*scaledSys;

Gp.OutputName = {'xm_VCA1', 'Pm_PC1', 'xm_VCA2', 'Pm_PC2', 'Fm_GP1'};
%display(Gp);

%% Task 3.3.3. Visualize uncertain system
fprintf('Task 3.3.3. Visualize uncertain system\n\n');
opt = stepDataOptions('StepAmplitude',0.001);

figure(5);
step(scaledSys,'r.',Gp,'b--',7,opt)
xlabel('time'); ylabel('response'); legend;
title('Step response with amplitude of 0.001');

figure(6);
step(scaledSys,'r.',Gp,'b--',7)
xlabel('time'); ylabel('response'); legend;
title('Step response with amplitude of 1');

%% Task 3.4. H-inf-controller synthesis
fprintf('Task 3.4. H-inf-controller synthesis\n');

% Define the objective values for H-inf-controller to minimize
%{
syms x
sigma = 0.0137 / 3;
f1 = 1 - 1/(sigma*sqrt(2*pi))*exp(-0.5*x^2/sigma^2)/72;
f2 = 0.5-tanh(x)/2;
Wxvca1, Wxvca2 = laplace(f1);
Wfgp1 = laplace(f2);
%}

s = tf('s');
Wppc1 = s;
Wppc2 = s;
Wxvca1 = ss(0); Wxvca2 = ss(0); Wfgp1 = ss(0);
%Wxvca1 = 1/s - (52252789446435703*2^(1/2)*21298081631339^(1/2)*pi^(1/2)*exp((6155145591456971*s^2)/590295810358705651712)*erfc((17*2^(1/2)*21298081631339^(1/2)*s)/34359738368))/87042659012253300578844672;
%Wxvca2 = 1/s - (52252789446435703*2^(1/2)*21298081631339^(1/2)*pi^(1/2)*exp((6155145591456971*s^2)/590295810358705651712)*erfc((17*2^(1/2)*21298081631339^(1/2)*s)/34359738368))/87042659012253300578844672;
%Wfgp1 = psi(s/4)/4 - psi(s/4 + 1/2)/4 + 1/s;
Wxvca1.u = 'xm_VCA1'; Wxvca1.y = 'e1';
Wxvca2.u = 'xm_VCA2'; Wxvca2.y = 'e2';
Wppc1.u = 'Pm_PC1'; Wppc1.y = 'e3';
Wppc2.u = 'Pm_PC2'; Wppc2.y = 'e4';
Wfgp1.u = 'Fm_GP1'; Wfgp1.y = 'e5';

ICinputs = {'U_GP1';'U_VCA1';'U_GP12';'U_GP2';'U_VCA2'}; %'F_VAD'
ICoutputs = {'e1';'e2';'e3';'e4';'e5';'xm_VCA1';'Pm_PC1';'xm_VCA2';'Pm_PC2';'Fm_GP1'};
mcl = connect(Gp, Wxvca1, Wxvca2, Wppc1, Wppc2, Wfgp1, ICinputs, ICoutputs);
ncont = 5; % 5 control signals
nmeas = 5; % 5 measurement signals
[K,~,gamma] = hinfsyn(mcl,nmeas,ncont);

display(K)
%% The End
fprintf('\nThank you for reviewing this report!\n');
