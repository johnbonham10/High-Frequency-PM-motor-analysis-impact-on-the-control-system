%% PROGETTO REGOLATORE 3

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

BW_des_prova = 10+03;                                                  % scalatura BW_des con approssimazione lineare

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

k_p_reg = 1.7268
k_d_reg = 1.1702e-05
k_i_reg = 63704.81
alpha = k_i_reg/k_p_reg
T = 1/alpha
omega_rit = 1/(T*sqrt(alpha))
R_ant2 = tf ([k_d_reg k_p_reg k_i_reg],[1 0])
R_ant3 = 1 ;

bode(R_ant2)

% funzioni di trasferimento con regolatore
L_ant = L_alg*R_ant2*R_ant3 ;
F_ant = feedback(L_ant,1) ;

%% step plot
% confronto tras regolatore canonico e progetto 1
step(F_ant)
hold on
step(F_algnc)
hold off
legend('Progetto 2','Schema di controllo canonico')

%% pzmap plot
% confronto tras regolatore canonico e progetto 1
pzmap(F_ant)
hold on
pzmap(F_alg)
hold off
legend('Progetto 2','Schema di controllo canonico')

%% bode plot
% confronto tras regolatore canonico e progetto 1
bode(L_ant)
hold on
bode(L_algnc)
hold off
legend('Progetto 2','Schema di controllo canonico')

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
legend('Progetto 2','Schema di controllo canonico') ;

%% rlocus analysis

L_alg_BW = L_alg*R_ant2*R_ant3/((ki_prova)/(2*pi*R_s)) ;                              % manipolo L_alg per far si che il gain sia BW_des
rlocus(L_alg_BW) ;

% estraggo i dati per il calcolo della deviazione
[r,k] = rlocus(L_alg_BW) ;                          % calcolo della deviazione
pole2 = -r(1,1:50)'/(2*pi) ;                                            % calcolo della deviazione  
bw2 = k(1,1:50)' ;                                                      % calcolo della deviazione da interpolare





