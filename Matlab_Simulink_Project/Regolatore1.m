%% PROGETTO REGOLATORE 1

%% calcolo zeri numeratore

qn = n0./n3 - (n1.*n2)./(3.*(n3.^(2))) + (2.*(n2.^(3)))/(27.*(n3.^(3))) ;               
pn = n1./n3 - (n2.^(2))/(3*((n3).^(2))) ;

deltan = (qn.^(2))./4 + (pn.^(3))./27 ;

un = -(abs(-qn./2 + (deltan).^(1/2))).^(1/3) ;
vn = -(abs(-qn./2 - (deltan).^(1/2))).^(1/3) ;
 
yn1 = un + vn ;
xn1 = -n2./(3*n3) + yn1 ;                                                           % calcolo z1

yn2 = un.*(-(1/2)+1i.*(sqrt(3)/2)) + vn.*(-(1/2)-1i.*(sqrt(3)/2)) ;
xn2 = (-n2./(3.*n3) + yn2 ) ;                                                       % calcolo z2

xn3 = real(xn2) - i*imag(xn2) ;                                                     % calcolo z3    

% smorzamenton = -real(xn2)./imag(xn2).*1./(sqrt(1+(real(xn2)./imag(xn2)).^(2))) ;
% omegan = real(xn2)./smorzamenton ;
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

%% regolatore canonico

% taratura regolatore
BW_des_provanc = 10e+03;
kp_provanc = 2*pi*BW_des_provanc*L_s ;
kidivkp_provanc = (R_s/L_s) ;                 
ki_provanc = kidivkp_provanc*kp_provanc ;

% calcolo parametri
n3nc = (C_g^(2))*R_g*L_s.*ki_provanc ;
n2nc = ((C_g*R_g)^2 + (C_g^(2))*R_g*R_s+C_g*L_s).*ki_provanc ;
n1nc = (2*C_g*R_g+C_g*R_s+(C_g^(2))*R_g*(omega_rads_ae^(2))*R_s).*ki_provanc ;
n0nc = (1+(C_g*R_g*omega_rads_ae)^(2)+(C_g^(2))*R_g*(omega_rads_ae)^(2)*R_s).*ki_provanc ;
d3nc = R_s*((C_g*R_g)^(2)) + n3nc ;
d2nc = R_s*2*C_g*R_g + n2nc ;
d1nc = R_s*(1+((C_g*R_g*omega_rads_ae)^(2))) + n1nc ;
d0nc = n0nc ;

%funzioni di trasferimento
L_algnc = tf([n3nc n2nc n1nc n0nc],[d3nc-n3nc d2nc-n2nc d1nc-n1nc d0nc-n0nc]) ;
F_algnc = tf([n3nc n2nc n1nc n0nc],[d3nc d2nc d1nc d0nc]) ;

%% regolatore progetto 1

BW_des_prova = 10e+03/(0.0001148);                                                  % scalatura BW_des con approssimazione lineare

%taratura regolatore
kp_prova = 2*pi*BW_des_prova*L_s ;
kidivkp_prova = (R_s/L_s) ;                 
ki_prova = kidivkp_prova*kp_prova ;

% calcolo coefficienti
n3 = (C_g^(2))*R_g*L_s.*ki_prova ;
n2 = ((C_g*R_g)^2 + (C_g^(2))*R_g*R_s+C_g*L_s).*ki_prova ;
n1 = (2*C_g*R_g+C_g*R_s+(C_g^(2))*R_g*(omega_rads_ae^(2))*R_s).*ki_prova ;
n0 = (1+(C_g*R_g*omega_rads_ae)^(2)+(C_g^(2))*R_g*(omega_rads_ae)^(2)*R_s).*ki_prova ;
d3 = R_s*((C_g*R_g)^(2)) + n3 ;
d2 = R_s*2*C_g*R_g + n2 ;
d1 = R_s*(1+((C_g*R_g*omega_rads_ae)^(2))) + n1 ;
d0 = n0 ;

% funzioni di trasferimento senza regolatore
L_alg = tf([n3 n2 n1 n0],[d3-n3 d2-n2 d1-n1 d0-n0]) ;
F_alg = tf([n3 n2 n1 n0],[d3 d2 d1 d0]) ;

% zero regolatore
f_zero = 30e+03 ;                                                                   % frequenza zero del regolatore
R_ant2 = tf ([1 f_zero*(2*pi)],[1]) ;

% cancellazione polo-zero
R_canc = tf ([1],[1 -xn2-xn3 xn2*xn3]) ;


% funzioni di trasferimento con regolatore
L_ant = L_alg*R_canc*R_ant2 ;
F_ant = feedback(L_ant,1) ;

%% step plot
% confronto tras regolatore canonico e progetto 1
step(F_ant)
hold on
step(F_algnc)
hold off
legend('Progetto 1','Schema di controllo canonico')

%% pzmap plot
% confronto tras regolatore canonico e progetto 1
pzmap(F_ant)
hold on
pzmap(F_alg)
hold off
legend('Progetto 1','Schema di controllo canonico')

%% bode plot
% confronto tras regolatore canonico e progetto 1
bode(L_ant)
hold on
bode(L_algnc)
hold off
legend('Progetto 1','Schema di controllo canonico')

%% nyquist
% confronto tras regolatore canonico e progetto 1
nyquist(L_ant)
hold on
nyquist(L_algnc)
n = 1000;                           
theta = linspace(0, 2*pi, n);
x = sqrt(2)*cos(theta);
y = sqrt(2)*sin(theta);
plot(x+1, y,'red'); 
plot(x-2, y,'red'); 
hold off
legend('Progetto 1','Schema di controllo canonico') ;

%% rlocus analysis

L_alg_BW = L_alg/((ki_prova)/(2*pi*R_s)) ;                              % manipolo L_alg per far si che il gain sia BW_des
rlocus(L_alg_BW*R_canc*R_ant2) ;

% estraggo i dati per il calcolo della deviazione
[r,k] = rlocus(L_alg_BW*R_ant*R_ant2*R_canc) ;                          % calcolo della deviazione
pole2 = -r(1,1:50)'/(2*pi) ;                                            % calcolo della deviazione  
bw2 = k(1,1:50)' ;                                                      % calcolo della deviazione da interpolare

%% grafico al variare dell'errore

% sovrastima del 10%
xn2_1 = (-n2./(3.*n3) + yn2 )*1.1 ;
xn3_1 = real(xn2_1) - i*imag(xn2_1) ;

R_canc_1 = tf ([1],[1 -xn2_1-xn3_1 xn2_1*xn3_1]) ;

L_ant_1 = L_alg*R_ant*R_canc_1*R_ant2 ;
F_ant_1 = feedback(L_ant_1,1) ;

% sovrastima del 30%
xn2_2 = (-n2./(3.*n3) + yn2 )*1.3 ;
xn3_2 = real(xn2_2) - i*imag(xn2_2) ;

R_canc_2 = tf ([1],[1 -xn2_2-xn3_2 xn2_2*xn3_2]) ;

L_ant_2 = L_alg*R_ant*R_canc_2*R_ant2 ;
F_ant_2 = feedback(L_ant_2,1) ;

% sottostima del 10%
xn2_4 = (-n2./(3.*n3) + yn2 )*0.9 ;
xn3_4 = real(xn2_4) - i*imag(xn2_4) ;

R_canc_4 = tf ([1],[1 -xn2_4-xn3_4 xn2_4*xn3_4]) ;

L_ant_4 = L_alg*R_ant*R_canc_4*R_ant2 ;
F_ant_4 = feedback(L_ant_4,1) ;

% sottostima del 30%
xn2_3 = (-n2./(3.*n3) + yn2 )*0.7 ;
xn3_3 = real(xn2_3) - i*imag(xn2_3) ;

R_canc_3 = tf ([1],[1 -xn2_3-xn3_3 xn2_3*xn3_3]) ;

L_ant_3 = L_alg*R_ant*R_canc_3*R_ant2 ;
F_ant_3 = feedback(L_ant_3,1) ;

% plot del confronto
step(F_ant_1)
hold on
step(F_ant_2)
step(F_ant_4)
step(F_ant_3)
legend('sovrastima del 10%','sovrastima del 30%','sottostima del 10%','sottostima del 30%')
hold off












