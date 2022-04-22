%% Root locus
 
clear;
clc;
close all;
 
%% parametri 
 
R = 4.3 ;  % ohm
L = 0.027750000000000 ; %H
 
F_load = tf([1],[L R]) ; % TF del carico
 
%% regolatore PI
 
BW_des=10000/10; % Hz
 
kp = 2*pi*BW_des*L
 
% m = 1.05 ;
ki_div_kp = (R/L)  % per la cancellazione zp
 
F_pi = tf([1 ki_div_kp],[1 0])  % TF PI senza kp
 
%% fdt open loop
 
G = F_load*F_pi
F = feedback(G,1)
%pzmap(G)
rlocus(G)
%bode(G)