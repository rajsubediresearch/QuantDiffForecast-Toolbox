% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cadfilename1, caddisease, datatype, dist1, numstartpoints, B, model, params, vars, windowsize1, tstart1, tend1, printscreen1] = options_fit_GRM_jalisco_measles_dist1_3

global method1

cadfilename1='JALISCO_2025-08-18_2026-06-29-trimmed';
caddisease='measles';
datatype='cases';

method1=3;   % MLE Neg Binomial (VAR = mean + alpha*mean)
dist1=3;
switch method1
    case 1
        dist1=1;
    case 3
        dist1=3;
    case 4
        dist1=4;
    case 5
        dist1=5;
end

numstartpoints=50;
B=300;

model.fc=@GRM;
model.name='GRM model';

params.num=4; % number of model parameters
params.phenom=1; % enable data-driven seed + bound scaling (phenomenological auto-config)
params.label={'r','p','a','K'};
params.LB=[0    0    0     1000];
params.UB=[10   1.2  10    10000];
params.initial=[0.79 0.84 1.3 7300];
params.fixed=[0 0 0 0];
params.fixI0=1;
params.composite='';
params.composite_name='';
params.extra0=[];

vars.label={'C'};
vars.initial=1;
vars.fit_index=1;
vars.fit_diff=1;

windowsize1=40;
tstart1=1;
tend1=1;
printscreen1=1;

end
