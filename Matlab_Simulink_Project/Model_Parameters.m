%% 

clear
clc
close all
%machine inertia (0.00179) 

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

L = [Ld+Lm,    -Lm*(1/2), -Lm*(1/2)  ;
    -Lm*(1/2), Ld+Lm,     -Lm*(1/2)  ;
    -Lm*(1/2), -Lm*(1/2),  Ld+Lm   ] ;

L_inv = inv(L) ;

R = [R_s, 0,   0    ;
     0,   R_s, 0    ;
     0,   0,   R_s] ;

L_s = Ld + (3/2)* Lm ;

Rg = [R_g, 0,   0     ;
      0  , R_g, 0     ;
      0  , 0  , R_g]  ;

Rg_inv = inv(Rg) ;

Cg = [C_g, 0,   0     ;
      0  , C_g, 0     ;
      0  , 0  , C_g]  ;

Cg_inv = inv(Cg) ;

%% PI Regulator
 
BW_des=200 ;                          % Desired Bandwidht [Hz]
 
kp = 2*pi*BW_des*L_s ;
kidivkp = (R_s/L_s) ;                 % pole-zero cancellation
ki = kidivkp*kp ;
% F_motor = tf([1],[L_s R_s]) ;       % motor tranfer function
% F_PI = tf([1 kidivkp],[1 0])        % PI tranfer fucntion for root locus
% G = F_load*F_pi
% F = feedback(G,1)
% pzmap(G)
% rlocus(G)
% bode(G)

%% Speed Observer

% desired frequency of the speed observer poles
p1 = 2*pi*300 ;
p2 = 2*pi*30 ;
p3 = 2*pi*3 ;
 
kp_so = (p1*p2 + p2*p3 + p1*p3)*J/p ;
ki_so = p1*p2*p3*J/(p) ;
kd_so = (p1 + p2 + p3)/(p) ;
 
kidivkp_so = ki_so/kp_so ;
 
%G = tf([P*K_d*J P*K_p K_i*P],[2*J P*J*K_d P*K_p P*K_i])
%pzmap(G)
%bode(G)
%rlocus(G)

%% Operative Conditions

Tload = 2 ;                         % Load Torque [N*m]

omega_rpm = 5000 ;                  % constan-speed load speed [rpm]
omega_rads = omega_rpm*(2*pi)/60 ;                  
omega_rads_ae = omega_rads*p ; 
f_ae = omega_rads_ae*7/(2*pi) ;

Iq = 5 ;                            % [A]
Id = 0 ;                            % [A]

T_des = 10 ;

harmonics = 0 ;                     % if you want to add the 5th e 7th harmonics to e, set this variable as 1
fifth_module = lambdamax/10 ;       % peak value of the 5th harmonic
seventh_module = lambdamax/20 ;     % peak value of the 7th harmonic

%% HF analysis

% iq''(vq)
numerator1 = [(C_g^2 * R_g), C_g, C_g^2 * R_g * omega_rads_ae^2];
denominator1 = [(C_g^2 * R_g^2),(2*C_g*R_g),1+C_g^2 * R_g^2 * omega_rads_ae^2];
iq_hf = tf(numerator1,denominator1);
%pzmap(iq_hf)
%bode(iq_hf)

% iq'(vq)
numerator5 = [1] ;
denominator5 = [L_s R_s] ;
iq_lf = tf(numerator5,denominator5) ;
%pzmap(iq_lf)
%bode(iq_lf)

% iq(vq)
iq = iq_lf + iq_hf ;
%pzmap(iq)
%bode(iq) ;
%step(iq)

%% Cloed Loop Analysis

G = iq ;                            % process transfer function
R = tf([kp,ki],[1 0]) ;             % controller transfer function

L = G*R ;                           % open loop tranfer function
F = L/(1+L) ;                       % cloed loop tranfer function

%% zero-pole map
%pzmap(F)

%% bode
%bode(F)

%% nyquist
%nyquist(L)
n = 1000;                           %// Define number of points on circle
theta = linspace(0, 2*pi, n);
x = sqrt(2)*cos(theta);
y = sqrt(2)*sin(theta);
%hold on;
%plot(x+1, y,'red'); %// Unit circle
%plot(x-2, y,'red'); %// Unit circle
%hold off;

%% step
%step(F)

%% studio F(s)

BW_des_prova=200 ;                          %linspace(200,20000,10000) Desired Bandwidht [Hz]
%% 
% taratura PI
kp_prova = 2*pi*BW_des_prova*L_s ;
kidivkp_prova = (R_s/L_s) ;                 % pole-zero cancellation
ki_prova = kidivkp_prova*kp_prova ;

% coefficienti F(s)
n3 = (C_g^(2))*R_g*L_s.*ki_prova ;
n2 = ((C_g*R_g)^2 + (C_g^(2))*R_g*R_s+C_g*L_s).*ki_prova ;
n1 = (2*C_g*R_g+C_g*R_s+(C_g^(2))*R_g*(omega_rads_ae^(2))*R_s).*ki_prova ;
n0 = (1+(C_g*R_g*omega_rads_ae)^(2)+(C_g^(2))*R_g*(omega_rads_ae)^(2)*R_s).*ki_prova ;
d3 = R_s*((C_g*R_g)^(2)) + n3 ;
d2 = R_s*2*C_g*R_g + n2 ;
d1 = R_s*(1+((C_g*R_g*omega_rads_ae)^(2))) + n1 ;
d0 = n0 ;

L_alg = tf([n3 n2 n1 n0],[d3-n3 d2-n2 d1-n1 d0-n0]) ;
F_alg = tf([n3 n2 n1 n0],[d3 d2 d1 d0]) ;

% %% rlocus analysis
% L_alg_BW = L_alg/(ki_prova/(2*pi*R_s)) ;
% rlocus(L_alg_BW) ;
% [r,k] = rlocus(L_alg_BW) ;                          % calcolo della deviazione
% pole1 = -r(1,1:47)'/(2*pi) ;                        % calcolo della deviazione  
% bw = k(1,1:47)' ;                                   % calcolo della deviazione da interpolare con un polinomio del 9 ordine
% plot(bw,-pole1)


%% pzmap 
%pzmap(F_alg)
%% step
%step(F_alg)
% %%
% qn = n0./n3 - (n1.*n2)./(3.*(n3.^(2))) + (2.*(n2.^(3)))/(27.*(n3.^(3))) ;
% pn = n1./n3 - (n2.^(2))/(3*((n3).^(2))) ;
% 
% deltan = (qn.^(2))./4 + (pn.^(3))./27 ;
% 
% un = -(abs(-qn./2 + (deltan).^(1/2))).^(1/3) ;
% vn = -(abs(-qn./2 - (deltan).^(1/2))).^(1/3) ;
%  
% yn1 = un + vn ;
% 
% xn1 = -n2./(3*n3) + yn1 ;
% 
% yn2 = un.*(-(1/2)+1i.*(sqrt(3)/2)) + vn.*(-(1/2)-1i.*(sqrt(3)/2)) ;
% 
% xn2 = -n2./(3.*n3) + yn2 ;
% 
% smorzamenton = -real(xn2)./imag(xn2).*1./(sqrt(1+(real(xn2)./imag(xn2)).^(2))) ;
% omegan = real(xn2)./smorzamenton ;
% 
% qd = d0./d3 - (d1.*d2)./(3*(d3.^(2))) + (2*(d2.^(3)))/(27*(d3.^(3))) ;
% pd = d1./d3 - (d2.^(2))./(3*(d3).^(2)) ;
% 
% deltad = (qd.^(2))./4 + (pd.^(3))./27 ;
% 
% ud = -(abs(-qd./2 + (deltad).^(1/2))).^(1/3) ;
% vd = -(abs(-qd./2 - (deltad).^(1/2))).^(1/3) ;
% 
% yd1 = ud + vd ;
% 
% xd1 = -d2./(3*d3) + yd1 ;
% 
% yd2 = ud.*(-(1/2)+1i.*(sqrt(3)/2)) + vd.*(-(1/2)-1i.*(sqrt(3)/2)) ;
% 
% xd2 = -d2./(3.*d3) + yd2 ;
% 
% smorzamentod = abs(-real(xd2)./imag(xd2).*1./(sqrt(1+(real(xd2)./imag(xd2)).^(2)))) ;
% smorzamentod(isnan(smorzamentod))=1 ;
% omegad = -real(xd2)./smorzamentod ;
% 
% % plot(BW_des_prova/1000,xd1)
% %%
% plot(BW_des_prova/1000,xd1/(2*pi*1000))
% hold on
% plot(BW_des_prova/1000,xn1/(2*pi*1000))
% title('Real Pole and Zero')
% xlabel('Frequency [kHz]')
% ylabel('Frequency [kHz]')
% legend('pole','zero')
