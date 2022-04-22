%% SPEED OBSERVER
% machine parameters
 
J = 0.00179
P = 4
 
%% Speed observer parameters
p1 = 2*pi*300
p2 = 2*pi*30
p3 = 2*pi*3
 
%% Transfer function
 
K_p = (p1*p2 + p2*p3 + p1*p3)*J/(P/2)
K_i = p1*p2*p3*J/(P/2)
K_d = (p1 + p2 + p3)/(P/2)
 
Ki_div_kp = K_i/K_p
 
G = tf([P*K_d*J P*K_p K_i*P],[2*J P*J*K_d P*K_p P*K_i])
pzmap(G)
%bode(G)
%rlocus(G)
%% SPEED REGULATOR
% Transfer function
 
tau_i = 1/(2*pi*200)
K_p_omega = J/(2*tau_i)
K_i_omega = J/(8*(tau_i)^2)
K_i_omega/K_p_omega 
F = tf([K_p_omega K_i_omega],[tau_i*J J 0 0])
W = tf([K_p_omega/K_i_omega 1],[J*tau_i/K_i_omega J/K_i_omega K_p_omega/K_i_omega  1])
%rlocus(F)
%pzmap(W)
%step(G)
%figure
bode(W)
 
% pzplot(G)
%% Current transfer function
%G = tf([1],[tau_i 1])
%bode(G)