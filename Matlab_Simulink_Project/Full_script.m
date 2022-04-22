%% Studio in alta frequenza dell azionamento di un motore sincrono a magneti permanenti: impatto sul sistema di controllo
%% High frequency analysis of the electrical drive of a PMSM: impact on the control system

% % % % % % % % % % % % % %
%                         %
%   author:               %
%   johnbonham1           %
%                         %  
% % % % % % % % % % % % % %
%% Machine Parameters

% LF parameters
lambdamax = 0.272 ;                 % peak value of PM flux [Wb]
p = 4 ;                             % n. of poles couples
Ld = 0.027 ;                        % Leakage Inductance [H]
Lm = 0.0005 ;                       % Magnetizing Inductance[H]
J =   1000000000  ;                 % setted to high value if a constant speed load is used
R_s = 4.3 ;                         % Resistance of the stator windings [ohm]
kt = (3/2)*(p)*lambdamax ;          % torque/q-current coefficient 

% HF parameters
R_g = 150 ;                         % high frequency resistance [ohm]
C_g = 3.92e-9 ;                     % high frequency capacitance [F]

%% Matrix builder and computations

% stator inductance abc frame
L = [Ld+Lm,    -Lm*(1/2), -Lm*(1/2)  ;      
    -Lm*(1/2), Ld+Lm,     -Lm*(1/2)  ;
    -Lm*(1/2), -Lm*(1/2),  Ld+Lm   ] ;

L_inv = inv(L) ;

% stator resistance
R = [R_s, 0,   0    ;
     0,   R_s, 0    ;
     0,   0,   R_s] ;

L_s = Ld + (3/2)* Lm ;

% ground resistance
Rg = [R_g, 0,   0     ;
      0  , R_g, 0     ;
      0  , 0  , R_g]  ;

Rg_inv = inv(Rg) ;

% ground capacitance
Cg = [C_g, 0,   0     ;
      0  , C_g, 0     ;
      0  , 0  , C_g]  ;

Cg_inv = inv(Cg) ;

%% PI Regulator
 
BW_des=200 ;                          % Desired Bandwidht [Hz]
 
kp = 2*pi*BW_des*L_s ;
kidivkp = (R_s/L_s) ;                 % pole-zero cancellation
ki = kidivkp*kp ;

%% Speed Observer

% desired frequency of the speed observer poles
p1 = 2*pi*300 ;
p2 = 2*pi*30 ;
p3 = 2*pi*3 ;
 
kp_so = (p1*p2 + p2*p3 + p1*p3)*J/p ;
ki_so = p1*p2*p3*J/(p) ;
kd_so = (p1 + p2 + p3)/(p) ;
 
kidivkp_so = ki_so/kp_so ;
 

%% Operative Conditions
% In this section will be defined the operative conditions of the Simulink
% simulation

Tload = 2 ;                         % Load Torque [N*m]

omega_rpm = 5000 ;                  % constant-speed load speed [rpm]
omega_rads = omega_rpm*(2*pi)/60 ;  % constant-speed load speed [rad/s]               
omega_rads_ae = omega_rads*p ;      % electrical speed
f_ae = omega_rads_ae*7/(2*pi) ;     % electrical frequency

Iq = 5 ;                            % q-current [A]
Id = 0 ;                            % d-current [A]

T_des = 10 ;                        % torque reference

%% Analysis with harmonics of the emf
% this section modifies parameters of the Simulink simulation

harmonics = 0 ;                     % if you want to add the 5th e 7th harmonics to e, set this variable as 1
fifth_module = lambdamax/10 ;       % peak value of the 5th harmonic
seventh_module = lambdamax/20 ;     % peak value of the 7th harmonic

%% HF analysis
% In this section transfer functions will be calcuated

% High Frequency (HF) transfer function iq''(vq)
numerator1 = [(C_g^2 * R_g), C_g, C_g^2 * R_g * omega_rads_ae^2] ;
denominator1 = [(C_g^2 * R_g^2),(2*C_g*R_g),1+C_g^2 * R_g^2 * omega_rads_ae^2] ;
iq_hf = tf(numerator1,denominator1) ;

% Low Frequency (LF) transfer function iq'(vq)
numerator5 = [1] ;
denominator5 = [L_s R_s] ;
iq_lf = tf(numerator5,denominator5) ;


% Total stator current transfer function iq(vq)
iq = iq_lf + iq_hf ;

%% Cloed Loop Analysis

G = iq ;                            % process transfer function
R = tf([kp,ki],[1 0]) ;             % canonical controller transfer function

L = G*R ;                           % open loop tranfer function
F = L/(1+L) ;                       % cloed loop tranfer function

%% nyquist
% In this section the nyquist analysis is performed

% plot of the circles that define the BW
n = 1000;                           
theta = linspace(0, 2*pi, n);       
x = sqrt(2)*cos(theta);
y = sqrt(2)*sin(theta);

%plot(x+1, y,'red'); %// Unit circle
%plot(x-2, y,'red'); %// Unit circle

%% Analysis of the open-loop transfer function F(s)

% PI regulator
BW_des_prova=200 ;                          % [Hz]
kp_prova = 2*pi*BW_des_prova*L_s ;
kidivkp_prova = (R_s/L_s) ;                 % pole-zero cancellation
ki_prova = kidivkp_prova*kp_prova ;

% F(s) coefficients, explicit calculus
% numerator
n3 = (C_g^(2))*R_g*L_s.*ki_prova ;
n2 = ((C_g*R_g)^2 + (C_g^(2))*R_g*R_s+C_g*L_s).*ki_prova ;
n1 = (2*C_g*R_g+C_g*R_s+(C_g^(2))*R_g*(omega_rads_ae^(2))*R_s).*ki_prova ;
n0 = (1+(C_g*R_g*omega_rads_ae)^(2)+(C_g^(2))*R_g*(omega_rads_ae)^(2)*R_s).*ki_prova ;
% denominator
d3 = R_s*((C_g*R_g)^(2)) + n3 ;
d2 = R_s*2*C_g*R_g + n2 ;
d1 = R_s*(1+((C_g*R_g*omega_rads_ae)^(2))) + n1 ;
d0 = n0 ;

% open-loop transfer function
L_alg = tf([n3 n2 n1 n0],[d3-n3 d2-n2 d1-n1 d0-n0]) ;

% closed-loop transfer function
F_alg = tf([n3 n2 n1 n0],[d3 d2 d1 d0]) ;

%% rlocus analysis
% In this section an rlocus analysis is performed

% The L(s) tf is manipulated with the aim of considering the desired BW as the variable gain 
L_alg_BW = L_alg/(ki_prova/(2*pi*R_s)) ;

% rlocus(L_alg_BW) ;
% [r,k] = rlocus(L_alg_BW) ;                          % deviation calculus
% pole1 = -r(1,1:47)'/(2*pi) ;                        % deviation calculus  
% bw = k(1,1:47)' ;                                  
% plot(bw,-pole1)

%% explicit calculus of zeros of L(s)

qn = n0./n3 - (n1.*n2)./(3.*(n3.^(2))) + (2.*(n2.^(3)))/(27.*(n3.^(3))) ;               
pn = n1./n3 - (n2.^(2))/(3*((n3).^(2))) ;

deltan = (qn.^(2))./4 + (pn.^(3))./27 ;

un = -(abs(-qn./2 + (deltan).^(1/2))).^(1/3) ;
vn = -(abs(-qn./2 - (deltan).^(1/2))).^(1/3) ;
 
yn1 = un + vn ;
xn1 = -n2./(3*n3) + yn1 ;                                                           % z1 

yn2 = un.*(-(1/2)+1i.*(sqrt(3)/2)) + vn.*(-(1/2)-1i.*(sqrt(3)/2)) ;
xn2 = (-n2./(3.*n3) + yn2 ) ;                                                       % z2

xn3 = real(xn2) - i*imag(xn2) ;                                                     % z3    

% complex zeros dampening and omega res 
smorzamenton = -real(xn2)./imag(xn2).*1./(sqrt(1+(real(xn2)./imag(xn2)).^(2))) ;
omegan = real(xn2)./smorzamenton ;

%% explicit calculus of poles of L(s)

qd = d0./d3 - (d1.*d2)./(3*(d3.^(2))) + (2*(d2.^(3)))/(27*(d3.^(3))) ;
pd = d1./d3 - (d2.^(2))./(3*(d3).^(2)) ;

deltad = (qd.^(2))./4 + (pd.^(3))./27 ;

ud = -(abs(-qd./2 + (deltad).^(1/2))).^(1/3) ;
vd = -(abs(-qd./2 - (deltad).^(1/2))).^(1/3) ;

yd1 = ud + vd ;

xd1 = -d2./(3*d3) + yd1 ;                                                           % d1

yd2 = ud.*(-(1/2)+1i.*(sqrt(3)/2)) + vd.*(-(1/2)-1i.*(sqrt(3)/2)) ;

xd2 = -d2./(3.*d3) + yd2 ;                                                          % d2

% complex poles dampening and omega res
smorzamentod = abs(-real(xd2)./imag(xd2).*1./(sqrt(1+(real(xd2)./imag(xd2)).^(2)))) ;
smorzamentod(isnan(smorzamentod))=1 ;
omegad = -real(xd2)./smorzamentod ;

%% Canonical Regulator

% PI design G(s)
BW_des_provanc = 10e+03;
kp_provanc = 2*pi*BW_des_provanc*L_s ;
kidivkp_provanc = (R_s/L_s) ;                 
ki_provanc = kidivkp_provanc*kp_provanc ;

% numerator of L(s)*G(s)
n3nc = (C_g^(2))*R_g*L_s.*ki_provanc ;
n2nc = ((C_g*R_g)^2 + (C_g^(2))*R_g*R_s+C_g*L_s).*ki_provanc ;
n1nc = (2*C_g*R_g+C_g*R_s+(C_g^(2))*R_g*(omega_rads_ae^(2))*R_s).*ki_provanc ;
n0nc = (1+(C_g*R_g*omega_rads_ae)^(2)+(C_g^(2))*R_g*(omega_rads_ae)^(2)*R_s).*ki_provanc ;
% denominator of L(s)*G(s)
d3nc = R_s*((C_g*R_g)^(2)) + n3nc ;
d2nc = R_s*2*C_g*R_g + n2nc ;
d1nc = R_s*(1+((C_g*R_g*omega_rads_ae)^(2))) + n1nc ;
d0nc = n0nc ;

% open-loop L(s)*G(s) transfer function
L_algnc = tf([n3nc n2nc n1nc n0nc],[d3nc-n3nc d2nc-n2nc d1nc-n1nc d0nc-n0nc]) ;

% closed-loop F(s) of L(s)*G(s) transfer function
F_algnc = tf([n3nc n2nc n1nc n0nc],[d3nc d2nc d1nc d0nc]) ;

%% PROJECT 1

BW_des_prova = 10e+03/(0.0001148);                          % linear scaling of BW

% PI design G(s)
kp_prova = 2*pi*BW_des_prova*L_s ;
kidivkp_prova = (R_s/L_s) ;                 
ki_prova = kidivkp_prova*kp_prova ;

% numerator of L(s)*G(s)
n3 = (C_g^(2))*R_g*L_s.*ki_prova ;
n2 = ((C_g*R_g)^2 + (C_g^(2))*R_g*R_s+C_g*L_s).*ki_prova ;
n1 = (2*C_g*R_g+C_g*R_s+(C_g^(2))*R_g*(omega_rads_ae^(2))*R_s).*ki_prova ;
n0 = (1+(C_g*R_g*omega_rads_ae)^(2)+(C_g^(2))*R_g*(omega_rads_ae)^(2)*R_s).*ki_prova ;
% denominator of L(s)*G(s)
d3 = R_s*((C_g*R_g)^(2)) + n3 ;
d2 = R_s*2*C_g*R_g + n2 ;
d1 = R_s*(1+((C_g*R_g*omega_rads_ae)^(2))) + n1 ;
d0 = n0 ;

% open-loop L(s)*G(s) transfer function
L_alg = tf([n3 n2 n1 n0],[d3-n3 d2-n2 d1-n1 d0-n0]) ;
% closed-loop F(s) of L(s)*G(s) transfer function
F_alg = tf([n3 n2 n1 n0],[d3 d2 d1 d0]) ;

% complex poles for the pole-zero cancelation G1(s)
R_canc = tf ([1],[1 -xn2-xn3 xn2*xn3]) ;

% regulator additional pole G2(s)
f_zero = 30e+03 ;                                         % pole frequency [Hz]
R_ant2 = tf ([1 f_zero*(2*pi)],[1]) ;

% open-loop transfer function L'(s)=L(s)*G(s)*G1(s)*G2(s)
L_ant = L_alg*R_canc*R_ant2 ;
F_ant = feedback(L_ant,1) ;

%% PROJECT 2

BW_des_prova = 10e+03;                                                  

% PI design G(s)
kp_prova = 2*pi*BW_des_prova*L_s ;
kidivkp_prova = (R_s/L_s) ;                 
ki_prova = kidivkp_prova*kp_prova ;

% numerator of L(s)*G(s)
n3 = (C_g^(2))*R_g*L_s.*ki_prova ;
n2 = ((C_g*R_g)^2 + (C_g^(2))*R_g*R_s+C_g*L_s).*ki_prova ;
n1 = (2*C_g*R_g+C_g*R_s+(C_g^(2))*R_g*(omega_rads_ae^(2))*R_s).*ki_prova ;
n0 = (1+(C_g*R_g*omega_rads_ae)^(2)+(C_g^(2))*R_g*(omega_rads_ae)^(2)*R_s).*ki_prova ;
% denominator of L(s)*G(s)
d3 = R_s*((C_g*R_g)^(2)) + n3 ;
d2 = R_s*2*C_g*R_g + n2 ;
d1 = R_s*(1+((C_g*R_g*omega_rads_ae)^(2))) + n1 ;
d0 = n0 ;

% open-loop L(s)*G(s) transfer function
L_alg = tf([n3 n2 n1 n0],[d3-n3 d2-n2 d1-n1 d0-n0]) ;
% closed-loop F(s) of L(s)*G(s) transfer function
F_alg = tf([n3 n2 n1 n0],[d3 d2 d1 d0]) ;

% PI regulator transfer function G1(s)
k_p_reg = 1.7268 ;
k_i_reg = 17508.4592 ;
alpha = k_i_reg/k_p_reg ;
T = 1/alpha ;
omega_rit = 1/(T*sqrt(alpha)) ;
R_ant2 = tf ([k_p_reg k_i_reg],[1 0]) ;
R_ant3 = 1 ;

% open-loop transfer function L'(s)=L(s)*G(s)*G1(s)*G2(s)
L_ant = L_alg*R_ant2*R_ant3 ;
F_ant = feedback(L_ant,1) ;

%% PROJECT 3

BW_des_prova = 10+03;                                                  

% PI design G(s)
kp_prova = 2*pi*BW_des_prova*L_s ;
kidivkp_prova = (R_s/L_s) ;                 
ki_prova = kidivkp_prova*kp_prova ;

% numerator of L(s)*G(s)
n3 = (C_g^(2))*R_g*L_s.*ki_prova ;
n2 = ((C_g*R_g)^2 + (C_g^(2))*R_g*R_s+C_g*L_s).*ki_prova ;
n1 = (2*C_g*R_g+C_g*R_s+(C_g^(2))*R_g*(omega_rads_ae^(2))*R_s).*ki_prova ;
n0 = (1+(C_g*R_g*omega_rads_ae)^(2)+(C_g^(2))*R_g*(omega_rads_ae)^(2)*R_s).*ki_prova ;
% denominator of L(s)*G(s)
d3 = R_s*((C_g*R_g)^(2)) + n3 ;
d2 = R_s*2*C_g*R_g + n2 ;
d1 = R_s*(1+((C_g*R_g*omega_rads_ae)^(2))) + n1 ;
d0 = n0 ;

% open-loop L(s)*G(s) transfer function
L_alg = tf([n3 n2 n1 n0],[d3-n3 d2-n2 d1-n1 d0-n0]) ;
F_alg = tf([n3 n2 n1 n0],[d3 d2 d1 d0]) ;

% PID regulator transfer function G1(s)
k_p_reg = 1.7268 ;
k_d_reg = 1.1702e-05 ;
k_i_reg = 63704.81 ;
alpha = k_i_reg/k_p_reg ;
T = 1/alpha ;
omega_rit = 1/(T*sqrt(alpha)) ;
R_ant2 = tf ([k_d_reg k_p_reg k_i_reg],[1 0]) ;
R_ant3 = 1 ;

% open-loop transfer function L'(s)=L(s)*G(s)*G1(s)*G2(s)
L_ant = L_alg*R_ant2*R_ant3 ;
F_ant = feedback(L_ant,1) ;

%% step-response analysis

% comparison between the canonical regulator and the project 1 regulator
% step(F_ant)
% hold on
% step(F_algnc)
% hold off
% legend('Progetto 1','Schema di controllo canonico')

%% pzmap analysis

% comparison between the canonical regulator and the project 1 regulator
% pzmap(F_ant)
% hold on
% pzmap(F_alg)
% hold off
% legend('Progetto 1','Schema di controllo canonico')

%% bode analysis

% comparison between the canonical regulator and the project 1 regulator
% bode(L_ant)
% hold on
% bode(L_algnc)
% hold off
% legend('Progetto 1','Schema di controllo canonico')

%% nyquist analysis

% comparison between the canonical regulator and the project 1 regulator
% nyquist(L_ant)
% hold on
% nyquist(L_algnc)
% n = 1000;                           
% theta = linspace(0, 2*pi, n);
% x = sqrt(2)*cos(theta);
% y = sqrt(2)*sin(theta);
% plot(x+1, y,'red'); 
% plot(x-2, y,'red'); 
% hold off
% legend('Progetto 1','Schema di controllo canonico') ;

%% rlocus analysis

% L_alg_BW = L_alg/((ki_prova)/(2*pi*R_s)) ;                              % manipolo L_alg per far si che il gain sia BW_des
% rlocus(L_alg_BW*R_canc*R_ant2) ;
% 
% % estraggo i dati per il calcolo della deviazione
% [r,k] = rlocus(L_alg_BW*R_ant*R_ant2*R_canc) ;                          % calcolo della deviazione
% pole2 = -r(1,1:50)'/(2*pi) ;                                            % calcolo della deviazione  
% bw2 = k(1,1:50)' ;                                                      % calcolo della deviazione da interpolare

%% influence of the parameter estimation error on the project 1 regulator performance

% overestimation of 10%
xn2_1 = (-n2./(3.*n3) + yn2 )*1.1 ;
xn3_1 = real(xn2_1) - i*imag(xn2_1) ;

R_canc_1 = tf ([1],[1 -xn2_1-xn3_1 xn2_1*xn3_1]) ;

L_ant_1 = L_alg*R_canc_1*R_ant2 ;
F_ant_1 = feedback(L_ant_1,1) ;

% overestimation of 30%
xn2_2 = (-n2./(3.*n3) + yn2 )*1.3 ;
xn3_2 = real(xn2_2) - i*imag(xn2_2) ;

R_canc_2 = tf ([1],[1 -xn2_2-xn3_2 xn2_2*xn3_2]) ;

L_ant_2 = L_alg*R_canc_2*R_ant2 ;
F_ant_2 = feedback(L_ant_2,1) ;

% underestimation of 10%
xn2_4 = (-n2./(3.*n3) + yn2 )*0.9 ;
xn3_4 = real(xn2_4) - i*imag(xn2_4) ;

R_canc_4 = tf ([1],[1 -xn2_4-xn3_4 xn2_4*xn3_4]) ;

L_ant_4 = L_alg*R_canc*R_ant2 ;
F_ant_4 = feedback(L_ant_4,1) ;

% underestimation of 30%
xn2_3 = (-n2./(3.*n3) + yn2 )*0.7 ;
xn3_3 = real(xn2_3) - i*imag(xn2_3) ;

R_canc_3 = tf ([1],[1 -xn2_3-xn3_3 xn2_3*xn3_3]) ;

L_ant_3 = L_alg*R_canc_3*R_ant2 ;
F_ant_3 = feedback(L_ant_3,1) ;

% comparison
% step(F_ant_1)
% hold on
% step(F_ant_2)
% step(F_ant_4)
% step(F_ant_3)
% legend('sovrastima del 10%','sovrastima del 30%','sottostima del 10%','sottostima del 30%')
% hold off




