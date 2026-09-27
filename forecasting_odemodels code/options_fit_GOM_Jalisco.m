
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cadfilename1,caddisease,datatype, dist1, numstartpoints,B, model,params,vars,windowsize1,tstart1,tend1,printscreen1]=options_fit_GOM_Jalisco

% OPTIONS_FIT_GOM_JALISCO  Configure QuantDiffForecast fitting of the Gompertz growth model
% to weekly measles case data from Jalisco (2025-08-18 to 2026-06-29).
%
% Overview
%   Returns all configuration needed by Run_Fit_ODEModel to calibrate the
%   Gompertz growth model and quantify uncertainty (parametric bootstrap).
%   The resulting output is then plotted with plotFit_ODEModel.
%
% Usage
%   [cadfilename1, caddisease, datatype, dist1, numstartpoints, B, ...
%    model, params, vars, windowsize1, tstart1, tend1, printscreen1] = ...
%    options_fit_GOM_Jalisco;
%
% Example run (fit first, then plot with the SAME tstart/tend/window):
%   Run_Fit_ODEModel(@options_fit_GOM_Jalisco, 1, 1, 20);
%   plotFit_ODEModel(@options_fit_GOM_Jalisco, 1, 1, 20);
%
% Model parameters
%   r = growth rate, a = exponential decay rate of growth
%
% NOTE: fit option files return 13 outputs (no getperformance or
% forecastingperiod). Run_Fit_ODEModel sets forecastingperiod=0 internally.
%
% See options_fit_SEIR_flu1918.m for full documentation of every field.
%
% See also
%   Run_Fit_ODEModel, plotFit_ODEModel, options_forecast_GOM_Jalisco


% <============================================================================>
% <=================== Declare global variables =======================================>
% <============================================================================>

global method1 % Parameter estimation method

% <============================================================================>
% <================================ Datasets properties =======================>
% <============================================================================>
% Located in the input folder, the time series data file is a text file with extension *.txt. 
% The first column corresponds to time index: 0,1,2, ... and the second
% column corresponds to the observed time series data.

cadfilename1='JALISCO_2025-08-18_2026-06-29-trimmed'; % Name of the data file (without extension)

caddisease='measles'; % string indicating the name of the disease related to the time series data

datatype='cases'; % string indicating the nature of the data (cases, deaths, hospitalizations, etc)

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

method1=3; % Type of estimation method

% Nonlinear least squares (LSQ)=0,
% MLE Poisson=1,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

dist1=3; % Define dist1 which is the type of error structure. See below:

%dist1=0; % Normal distribution to model error structure (method1=0)
%dist1=1; % Poisson error structure (method1=0 OR method1=1)
%dist1=2; % Neg. binomial error structure where var = factor1*mean where
                  % factor1 is empirically estimated from the time series
                  % data (method1=0)
%dist1=3; % MLE (Neg Binomial) with VAR=mean+alpha*mean  (method1=3)
%dist1=4; % MLE (Neg Binomial) with VAR=mean+alpha*mean^2 (method1=4)
%dist1=5; % MLE (Neg Binomial)with VAR=mean+alpha*mean^d (method1=5)

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

numstartpoints=100; % Number of initial guesses for optimization procedure using MultiStart

B=300; % number of bootstrap realizations to characterize parameter uncertainty

% <==============================================================================>
% <============================== ODE model =====================================>
% <==============================================================================>

model.fc=@GOM; % name of the model function
model.name='Gompertz model';   % string indicating the name of the ODE model (used in output file names)

params.num=2; % number of model parameters
params.label={'r','a'};  % list of symbols to refer to the model parameters
params.LB=[0.05 0.05]; % lower bound values of the parameter estimates
params.UB=[5 1]; % upper bound values of the parameter estimates
params.initial=[0.5 0.5]; % initial parameter values/guesses
params.fixed=[0 0]; % Boolean vector: 1 = fixed at params.initial, 0 = estimated
params.fixI0=1; % 1 = initial value of the fitting variable fixed to the first observation; 0 = estimated
params.composite=[];  % no composite parameter for this phenomenological model
params.composite_name=''; % name of the composite parameter (must exist, even if empty; plotFit_ODEModel prints it)
params.extra0=[]; % used to pass any extra parameters (e.g., data, static variables) to the model function

vars.num=1; % number of variables comprising the ODE model
vars.label={'C'}; % list of symbols to refer to the variables included in the model
vars.initial=1; % initial condition (first case count in the series; overridden when params.fixI0=1)
vars.fit_index=1; % index of the model's variable that will be fit to the observed time series data
vars.fit_diff=1; % 1 = fit the derivative of C (incidence) to the data

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

windowsize1=20;  % moving window size (calibration period, weeks)

tstart1=1; % time point for the start of rolling window analysis

tend1=1;  % time point for the end of the rolling window analysis

printscreen1=1; % print outputs/figures to screen (1) or not (0)
