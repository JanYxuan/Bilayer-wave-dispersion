% ========================================================================
% Description:
%   This is the main script configured and executed by the user.
%   the material and geometric parameters of the bilayer system, along with
%   the frequency points of interestcan be specified by the user;
%   the corresponding phase velocities are then computed and output.
%
% Usage:
%   Users can configure the relevant parameters in the Parameter Settings section. 
%   The wave dispersion (i.e. frequency - phase velocity relation) will be calculated.
%
% Author: Yuxuan Jiang, @ Wellman Center for Photomedicine
% Created date: 2025-12
% GitHub: https://github.com/JanYxuan/Bilayer-wave-dispersion.git
% Reference:
%   [1] X Feng, GY Li, Y Jiang, O Shortt-Nguyen, SH Yun. Optical coherence
%   elastography measures mechanical tension in the lens and capsule. Acta
%   Biomaterialia, 199:252-261, 2025.

clc,clear,close all

%% ------- Parameter Settings ---------
% ---- parameters for the plate (Layer #1)---
mu=500e3; % [Pa], shear modulus of layer 1
sigma=40e3; % [Pa], stress of layer 1 (caution: sigma<1.5*mu)
rho=1e3; % [kg/m^3], mass density of layer 1
h=55e-6; % [m], wall thickness

% --- parameters for the substrate (Layer #2)---
mu_s=5e3; % [Pa], shear modulus of substrate
sigma_s=0.4e3; % [Pa], stress of substrate (caution: sigma_s<1.5*mu_s)
rho_s=1e3; % [kg/m^3], mass density of substrate

% ---- frequency range ----
freq_list=0.1:500:10001; % [Hz], frequency

%% ---- Solving the wave dispersion ------
paras = [mu,sigma,mu_s,sigma_s,rho,rho_s,h];
phase_vel=SolveStressedBilayerWaveDispersion(paras,freq_list);

figure(1)
plot(freq_list/1e3,phase_vel)
xlabel('Frequency (kHz)')
ylabel('Phase Velocity (m/s)')
