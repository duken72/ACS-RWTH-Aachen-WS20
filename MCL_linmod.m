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


%% Task 3.2.5 Scaling
fprintf('\nTask 3.2.5. Scaling\n');
Du = diag([24,24,24,24,24]);
Dy = diag([0.0137,200,0.0137,200,100]);
Dd = eye(5);

TFrSysr_scaled = pinv(Dy)*TFrSysr*Du;
%display(TFrSysr_scaled);


%% Task 3.3.2. Weighting functions for uncertainty
fprintf('Task 3.3.2. Weighting functions for uncertainty\n\n');

G = minreal(ss(TFrSysr_scaled));
%Wxvca12 = makeweight(.1,5,1);
%Wppc12 = makeweight(.3,2.5,1);
%wppfgpq = makeweight(.1,1.5,1);
%Delta1 = ultidyn('Delta1',[1 1]);
%Delta2 = ultidyn('Delta2',[1 1]);
%Delta3 = ultidyn('Delta3',[1 1]);
%G = H*blkdiag(1+W1*Delta1,1+W2*Delta2,1+W1*Delta1,1+W2*Delta2, 1+w1*Delta3);

Delta = ultidyn('Delta',[5 5],'Bound',1);
Wxvca12 = makeweight(.1,1,5);
Wppc12 = makeweight(.3,1,2.5);
wppfgpq = makeweight(.1,1,1.5); 
W1 = diag([Wxvca12.A,Wppc12.A,Wxvca12.A,Wppc12.A,wppfgpq.A]);
W2 = diag([Wxvca12.B,Wppc12.B,Wxvca12.B,Wppc12.B,wppfgpq.B]);
W3 = diag([Wxvca12.C,Wppc12.C,Wxvca12.C,Wppc12.C,wppfgpq.C]);
W4 = diag([Wxvca12.D,Wppc12.D,Wxvca12.D,Wppc12.D,wppfgpq.D]);
W = ss(W1,W2,W3,W4);
E = Delta*W;
% Delta = ultidyn('Delta',[5 1],'Bound',1);
% s = tf('s');
% R01 = 0.1; R02 = 0.3; R03 = 0.1;
% Rinf1 = 5; Rinf2 = 2.5; Rinf3 = 1.5;
% tau1 = 1; tau2 = 2; tau3 = 1;
% 
% w1 = (tau1*s + R01)/(1+s*tau1/Rinf1);
% w2 = (tau2*s + R02)/(1+s*tau2/Rinf2);
% w3 = (tau3*s + R03)/(1+s*tau3/Rinf3);
% 
% W = [w1,w2,w1,w2,w3];
% E = Delta*W;
% for i=1:5
%     E(i,i+1:end)=0;
%     E(i,1:i-1)=0;
% end
Gp = G*(eye(5)+E);
set(Gp,'inputname',[{'U_GP1'} {'U_VCA1'} {'U_GP12'} {'U_GP2'} {'U_VCA2'}]);
set(Gp,'outputname',[{'xm_VCA1'} {'Pm_pc1'} {'xm_VCA2'} {'Pm_pc2'} {'Fm_GP1'}]);
bodemag(Gp);
%%

opt = stepDataOptions('StepAmplitude',0.001);
figure(6);
step(Gp,opt)
legend;
xlabel('time'); ylabel('response');

opt1 = stepDataOptions('StepAmplitude',1);
figure(7);
step(Gp,opt1)
legend;
xlabel('time'); ylabel('response');

%%
% for i=1:5
% [YppC,tppC] = step(Gpp(2,i));
% sserror(i)=abs(1-YppC(end));
% end
% 
% [M,Delta1] = lftdata(Gpp);
% [K,CL,gamma] = hinfsyn(M);
% Tinf = feedback(Gpp*K,1);
% figure(8);
% step(Tinf);

s = tf('s');
Wppc1 = 1/s;
Wppc2 = 1/s;
r1 = 0.01372; Wxvca2 = 0.01372; Wfgp1 = ss(0);
Sum1 = sumblk('e1 = r1 - x_VCA1',1);
Sum2 = sumblk('e3 = r1 - x_VCA2',1);
%Wxvca1 = 1/s - (52252789446435703*2^(1/2)*21298081631339^(1/2)*pi^(1/2)*exp((6155145591456971*s^2)/590295810358705651712)*erfc((17*2^(1/2)*21298081631339^(1/2)*s)/34359738368))/87042659012253300578844672;
%Wxvca2 = 1/s - (52252789446435703*2^(1/2)*21298081631339^(1/2)*pi^(1/2)*exp((6155145591456971*s^2)/590295810358705651712)*erfc((17*2^(1/2)*21298081631339^(1/2)*s)/34359738368))/87042659012253300578844672;
%Wfgp1 = psi(s/4)/4 - psi(s/4 + 1/2)/4 + 1/s;
%Wxvca1.u = 'xm_VCA1'; Wxvca1.y = 'e1';
%Wxvca2.u = 'xm_VCA2'; Wxvca2.y = 'e2';
Wppc1.u = 'Pm_pc1'; Wppc1.y = 'e2';
Wppc2.u = 'Pm_pc2'; Wppc2.y = 'e4';
Wfgp1.u = 'Fm_GP1'; Wfgp1.y = 'e5';

ICinputs = {'U_GP1';'U_VCA1';'U_GP12';'U_GP2';'U_VCA2'}; %'F_VAD'
ICoutputs = {'xm_VCA1';'Pm_pc1';'xm_VCA2';'Pm_pc2';'Fm_GP1'};
mcl = connect(Gp, Sum1, Sum2, Wppc1, Wppc2, Wfgp1, ICinputs, ICoutputs);
ncont = 5; % 5 control signals
nmeas = 5; % 5 measurement signals
[K,~,gamma] = hinfsyn(mcl,nmeas,ncont);

K.u={'xm_VCA1';'Pm_pc1';'xm_VCA2';'Pm_pc2';'Fm_GP1'};
K.y={'U_GP1';'U_VCA1';'U_GP12';'U_GP2';'U_VCA2'};

%Closed Loop

CL = connect(Gp.nominal,K,ICinputs,ICoutputs );

figure(8);
step(CL);
%% The End
fprintf('\nThank you for reviewing this report!\n');
