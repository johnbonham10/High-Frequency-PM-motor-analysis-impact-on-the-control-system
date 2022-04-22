%% HF Motor Model

% Parametric analysis

Rg_min = 0.1;
Rg_max = 400;
Cg_min = 1e-9;
Cg_max = 25e-9;

Lq = L_s;
Rs = 4.3;

R_g = 200;
C_g = 25e-9;

omega_ae = omega_rads_ae

CgRg = linspace(Rg_min * Cg_min,Rg_max * Cg_max,100000);

% %% omega zero
% omega_zero = sqrt(1 + omega_ae^2 * CgRg.^2)./ (CgRg);
% omega_zero1 = sqrt(1 + omega_ae^2 * (Cg_max *Rg_max)^2)./ (Cg_max*Rg_max)
% omega_zero2 = sqrt(1 + omega_ae2^2 * CgRg.^2)./ (CgRg);
% 
% f_zero = omega_zero / (2 * pi) / 1000000; %MHz
% f_zero1 = omega_zero1 / (2 * pi) / 1000 %MHz
% f_zero2 = omega_zero2 / (2 * pi) / 1000000; %MHz
% plot(CgRg,f_zero);
% grid on
% title("Frequency f_o")
% xlabel("Cg * Rg")
% ylabel("Frequency [MHz]")
% ylim([0,5])
% zeta = 1 ./ (1 + omega_ae^2 * (Cg_max*Rg_max)^2)^(1/2)
%% zeta
zeta = 1 ./ (1 + omega_ae^2 .* (CgRg).^2).^(1/2);
plot(CgRg,zeta);

grid on
title("xi")
xlabel("Cg * Rg")
ylabel("Dampening")
ylim([0,2])
%% test1
Rg = R_g;
Cg = C_g;

omega_zero_x = sqrt(1 + omega_ae^2 * (Cg * Rg)^2)/ (Cg * Rg);
f_zero_x = omega_zero_x / (2 * pi) / 1000
%% iq''(vq)
numerator1 = [(Cg^2 * Rg), Cg, Cg^2 * Rg * omega_ae^2]
denominator1 = [(Cg^2 * Rg^2),(2*Cg*Rg),1+Cg^2 * Rg^2 * omega_ae^2];
sys1 = tf(numerator1,denominator1)
pzmap(sys1)
%bode(sys1)
%step(sys1)
%% iq''(vd)
numerator2 = [2 * Cg^2 * Rg * omega_ae, Cg * omega_ae];
denominator2 = [(Cg^2 * Rg^2),(2*Cg*Rg),1+Cg^2 * Rg^2 * omega_ae^2];
sys2 = tf(numerator2,denominator2)
pzmap(sys2)
%% id''(vq)
numerator3 = [-Cg * omega_ae];
denominator3 = [(Cg^2 * Rg^2),(2*Cg*Rg),1+Cg^2 * Rg^2 * omega_ae^2];
sys3 = tf(numerator3,denominator3)
pzmap(sys3)
step(sys3)
%% id''(vd)
numerator4 = [(Cg^2 * Rg), Cg, Cg^2 * Rg * omega_ae^2];
denominator4 = [(Cg^2 * Rg^2),(2*Cg*Rg),1+Cg^2 * Rg^2 * omega_ae^2];
sys4 = tf(numerator4,denominator4)
pzmap(sys4)
%% plot all pzmap
pzmap(sys1,sys2,sys3,sys4)
%% iq'(vq)
numerator5 = [1]
denominator5 = [Lq Rs]
sys5 = tf(numerator5,denominator5)
pzmap(sys5)
%step(sys5)
%% iq(vq)
iq = sys1 + sys5
pzmap(iq)
%bode(iq)
%step(iq)
%% id'(vd)
numerator5 = [1]
denominator5 = [Lq Rs]
sys6 = tf(numerator5,denominator5)
pzmap(sys6)
%% id(vd)
id = sys4 + sys6
pzmap(id)
bode(iq)
%step(iq)