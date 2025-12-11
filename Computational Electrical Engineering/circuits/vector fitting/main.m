%% VECTOR FITTING ANALYSIS
clc;clear;close all;

% addpath(genpath("vector_fitting/")) % folder and all subfolders % UNCOMMENT
% TO LOAD THE LIBRARY 
load dataY_cap.mat
bigH = reshape(YY, [1, 1, length(freqs)]);
s = 2*pi*1j*freqs; 

% options of vf routine
opts=[]; opts.N = 2; opts.poletype = 'logcmplx'; opts.asymp = 2;
opts.stable=1; opts.passive_DE=1; poles = []; opts.cmplx_ss = 1;
[SERY,rmserr,bigYfit,opts2]=VFdriver(bigH,s,poles,opts);

A = full(SERY.A);
B = full(SERY.B);
C = full(SERY.C);
D = full(SERY.D);

[filename,Port_node_name,Pole_node_name,ground_node_name]= fun_VFmodel2Netlist_MultiPort(A,B,C,D,'netlist.m');

[lib_filename] = generate_LTspice_netlist_MultiPort(size(B,2), filename,'MyVFmodel',Port_node_name,ground_node_name);