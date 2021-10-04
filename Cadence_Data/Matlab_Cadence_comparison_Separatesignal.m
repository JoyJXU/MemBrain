
clc;
clear;
close all;
%% load cadence data
load zi_modified_sep 
load zj_modified_sep
load pi_modified_sep
load pj_modified_sep
load pij_modified_sep
load logpi_modified_sep
load logpj_modified_sep
load logpij_modified_sep
load wij_modified_sep

%% BCPNN parameter
kp = 1/500;
kz = 1/11;
kfti = 1/11;
eps = 0.01;

%% BCPNN input
points=5000;
si=zeros(1,points);
si(100:4:250)=1;
sj=zeros(1,points);
sj(1100:4:1250)=1;
%% BCPNN
zi(1)=0;
pi(1)=0;
zj(1)=0;
pj(1)=0;
pij(1)=0;
for i= 2: points
    pi(i) = pi(i-1) * (1 - kp) + zi(i-1)*kp;
    zi(i) = zi(i-1) * (1 - kz) + si(i-1)*kz; 
    pj(i) = pj(i-1) *( 1 - kp) + zj(i-1)*kp;
    zj(i) = zj(i-1) * (1 - kz) + sj(i-1)*kfti;
    pij(i) = pij(i-1) * (1 - kp) + zi(i-1)*zj(i-1)*kp;
    wij(i) = log((pij(i-1)+eps^2)/((pi(i-1)+eps)*(pj(i-1)+eps)));
    bj=log(pj+eps);
end

%% memristor model parameters
D = 1;
w_init = 0; 
Roff = 2e5;
Ron = 2000;
alpha_on = 1;
alpha_off = 1;
v_on=-0.02;
v_off=0.02;
delta_t=0.5e-3;
k_on= -28;
k_off = 21;
P_coeff = 1;
%% cadence circuit parameter
current=50e-9;
Ron_current=Ron*current;
%% Data processing parameter
trans_ratio_z_New = D/(Roff-Ron)/current;
trans_ratio_pi_New =1/19.5e5/current;
trans_ratio_pj_New =1/19.5e5/current; 
trans_ratio_pij_New=1/10e6/current;
rv_ratio_p_New = 0.23; 
rv_devia_p_New = -0.023; 
rv_ratio_pij_New = 1; 
rv_devia_pij_New = -0.023;
%% Initial value of the memristor model variable
New_X_i(1)=w_init*D;
New_X_j(1)=w_init*D;
New_PX_i(1)=0;
New_PX_j(1)=0;
New_PX_ij(1)=0;
W_ij(1) = 0;
%% Zi memristor
% memristor input
a=si;
a(a==0)=v_on*(-1e3/11/k_on*2+1);
a(a==1)=v_off*(1e3/11/k_off*2+1);
i_ = 50e-9.*ones(1,points);%sampling signal
New_rand_zi = [a;i_];
New_h_i= New_rand_zi(:)';
for p=2:2*points
    if  (New_h_i(p) > 0) && (New_h_i(p) > v_off)
        New_X_dot_i=k_off*(New_h_i(p)/v_off-1)^alpha_off;
        New_X_i(p)=New_X_i(p-1)+delta_t*New_X_dot_i.*(1-New_X_i(p-1)/D)^(P_coeff);
    elseif (New_h_i(p) <= 0) && (New_h_i(p) < v_on)
        New_X_dot_i=k_on*((New_h_i(p))/v_on-1)^alpha_on;
        New_X_i(p)=New_X_i(p-1)+delta_t*New_X_dot_i.*(New_X_i(p-1)/D)^(P_coeff);
    else
        New_X_i(p)=New_X_i(p-1);
        New_X_dot_i=0;
    end  
end
New_R_i=Roff.*New_X_i./D+Ron.*(1-New_X_i./D);
New_R_i=New_R_i*current;%The data sampled by the current in Cadence

% value to Zi
New_Z_i = (New_R_i(2:2:2*points)-Ron_current).*trans_ratio_z_New;
New_Z_i_shift = [0,New_Z_i(1:4999)];%memristor data right shift

%error analysis
New_Z_error_i = abs(New_Z_i_shift - zi);
New_Z_max_error_i = max(New_Z_error_i);
New_Z_mean_error_i = mean(New_Z_error_i);
corr_zi = corrcoef(zi,New_Z_i_shift);
New_Z_corr_i = corr_zi(2,1);
rrmse_zi= sqrt((sum((New_Z_i_shift - zi).^2)/(sum(zi.^2)))/5000);
%% Pi memristor
% memristor input
New_newR_i =New_Z_i*rv_ratio_p_New + rv_devia_p_New;
New_rand_pi = [New_newR_i;i_];
New_P_input_i = New_rand_pi(:)';

for p=2:2*points
    if  (New_P_input_i(p) > 0) && (New_P_input_i(p) > v_off)
        New_PX_dot_i=k_off*(New_P_input_i(p)/v_off-1)^alpha_off;
        New_PX_i(p)=New_PX_i(p-1)+delta_t*New_PX_dot_i.*(1-New_PX_i(p-1)/D)^(P_coeff);
    elseif (New_P_input_i(p) <= 0) && (New_P_input_i(p) < v_on)
        New_PX_dot_i=k_on*((New_P_input_i(p))/v_on-1)^alpha_on;
        New_PX_i(p)=New_PX_i(p-1)+delta_t*New_PX_dot_i.*(New_PX_i(p-1)/D)^(P_coeff);
    else
        New_PX_i(p)=New_PX_i(p-1);
        New_PX_dot_i=0;
    end
end 
New_PR_i=Roff.*New_PX_i./D+Ron.*(1-New_PX_i./D);
New_PR_i=New_PR_i*current;%The data sampled by the current in Cadence

% value to Pi
New_P_i = (New_PR_i(2:2:2*points)-Ron_current).*trans_ratio_pi_New;
New_P_i_shift = [0,New_P_i(1:4999)];%memristor right shift

% error analysis
New_P_error_i = abs(New_P_i_shift - pi);
New_P_max_error_i = max(New_P_error_i);
New_P_mean_error_i = mean(New_P_error_i);
corr_pi = corrcoef(pi,New_P_i_shift);
New_P_corr_i = corr_pi(2,1);
rrmse_pi= sqrt((sum((New_P_i_shift - pi).^2)/(sum(zi.^2)))/5000);


%% Zj memristor
% memristor input
b=sj;
b(b==0)=v_on*(-1e3/11/k_on*2+1);
b(b==1)=v_off*(1e3/11/k_off*2+1);
New_rand_zj = [b;i_];
New_h_j = New_rand_zj(:)';

for p=2:2*points
    if  (New_h_j(p) > 0) && (New_h_j(p) > v_off)
        New_X_dot_j=k_off*(New_h_j(p)/v_off-1)^alpha_off;
        New_X_j(p)=New_X_j(p-1)+delta_t*New_X_dot_j.*(1-New_X_j(p-1)/D)^(P_coeff);
    elseif (New_h_j(p) <= 0) && (New_h_j(p) < v_on)
        New_X_dot_j=k_on*((New_h_j(p))/v_on-1)^alpha_on;
        New_X_j(p)=New_X_j(p-1)+delta_t*New_X_dot_j.*(New_X_j(p-1)/D)^(P_coeff);
    else
        New_X_j(p)=New_X_j(p-1);
        New_X_dot_j=0;
    end  
end
New_R_j=Roff.*New_X_j./D+Ron.*(1-New_X_j./D);
New_R_j=New_R_j*current;%The data sampled by the current in Cadence

% value to Zj
New_Z_j = (New_R_j(2:2:2*points)-Ron_current).*trans_ratio_z_New;
New_Z_j_shift = [0,New_Z_j(1:4999)];%memristor right shift
%error analysis
New_Z_error_j = abs(New_Z_j_shift - zj);
New_Z_max_error_j= max(New_Z_error_j);
New_Z_mean_error_j = mean(New_Z_error_j);
corr_zj = corrcoef(zj,New_Z_j_shift);
New_Z_corr_j = corr_zj(2,1);
rrmse_zj= sqrt((sum((New_Z_j_shift - zj).^2)/(sum(zi.^2)))/5000);
%% Pj memristor
% memristor input
New_newR_j =New_Z_j*rv_ratio_p_New + rv_devia_p_New;
New_rand_pj = [New_newR_j;i_];
New_P_input_j = New_rand_pj(:)';

for p=2:2*points
    if  (New_P_input_j(p) > 0) && (New_P_input_j(p) > v_off)
        New_PX_dot_j=k_off*(New_P_input_j(p)/v_off-1)^alpha_off;
        New_PX_j(p)=New_PX_j(p-1)+delta_t*New_PX_dot_j.*(1-New_PX_j(p-1)/D)^(P_coeff);
    elseif (New_P_input_j(p) <= 0) && (New_P_input_j(p) < v_on)
        New_PX_dot_j=k_on*((New_P_input_j(p))/v_on-1)^alpha_on;
        New_PX_j(p)=New_PX_j(p-1)+delta_t*New_PX_dot_j.*(New_PX_j(p-1)/D)^(P_coeff);
    else
        New_PX_j(p)=New_PX_j(p-1);
        New_PX_dot_j=0;
    end
end 
New_PR_j=Roff.*New_PX_j./D+Ron.*(1-New_PX_j./D);
New_PR_j=New_PR_j*current;%The data sampled by the current in Cadence

% value to Pj
New_P_j = (New_PR_j(2:2:2*points)-Ron_current).*trans_ratio_pj_New;
New_P_j_shift = [0,New_P_j(1:4999)];%memristor right shift
%error analysis
New_P_error_j = abs(New_P_j_shift - pj);
New_P_max_error_j = max(New_P_error_j);
New_P_mean_error_j = mean(New_P_error_j);
corr_pj = corrcoef(pj,New_P_j_shift);
New_P_corr_j = corr_pj(2,1);
rrmse_pj= sqrt((sum((New_P_j_shift - pj).^2)/(sum(pj.^2)))/5000);
%% bj memristor
b_j=log(New_P_j+eps);
b_j_shift = [[-4.6017],b_j(1:4999)];%memristor right shift

%error analysis
New_b_error_j = abs(b_j_shift - bj);
New_b_max_error_j = max(New_b_error_j);
New_b_mean_error_j = mean(New_b_error_j);
corr_bj = corrcoef(b_j_shift,bj);
New_b_corr_j = corr_bj(2,1);
rrmse_bj= sqrt((sum((b_j_shift - bj).^2)/(sum(bj.^2)))/5000);
%% pij memristor
% memristor input
zij=New_Z_i.*New_Z_j;
New_Z_ij=(zij*rv_ratio_pij_New + rv_devia_pij_New);
New_rand_ij = [New_Z_ij;i_];
New_P_input_ij = New_rand_ij(:)';
for p=2:2*points
    if  (New_P_input_ij(p) > 0) && (New_P_input_ij(p) > v_off)
        New_PX_dot_ij=k_off*(New_P_input_ij(p)/v_off-1)^alpha_off;
        New_PX_ij(p)=New_PX_ij(p-1)+delta_t*New_PX_dot_ij.*(1-New_PX_ij(p-1)/D)^(P_coeff);
    elseif (New_P_input_ij(p) <= 0) && (New_P_input_ij(p) < v_on)
        New_PX_dot_ij=k_on*((New_P_input_ij(p))/v_on-1)^alpha_on;
        New_PX_ij(p)=New_PX_ij(p-1)+delta_t*New_PX_dot_ij.*(New_PX_ij(p-1)/D)^(P_coeff);
    else
        New_PX_ij(p)=New_PX_ij(p-1);
        New_PX_dot_ij=0;
    end
end 
New_PR_ij=Roff.*New_PX_ij./D+Ron.*(1-New_PX_ij./D);
New_PR_ij=New_PR_ij*current;%The data sampled by the current in Cadence

% value to Pij
New_P_ij = (New_PR_ij(2:2:2*points)-Ron_current).*trans_ratio_pij_New;
New_P_ij_shift = [0,New_P_ij(1:4999)];%memristor right shift

%error analysis
New_P_error_ij = abs(New_P_ij_shift - pij);
New_P_max_error_ij = max(New_P_error_ij);
New_P_mean_error_ij= mean(New_P_error_ij);
corr_pij = corrcoef(pij,New_P_ij_shift);
New_P_corr_ij = corr_pij(2,1);
rrmse_pij= sqrt((sum((New_P_ij_shift - pij).^2)/(sum(pij.^2)))/5000);
%% Wij memristor
for i=2:points

        W_ij(i) = log((New_P_ij(i-1)+eps^2)/((New_P_i(i-1)+eps).*(New_P_j(i-1)+eps)));

end

% value to wij
W_ij_shift=[0,W_ij(1:4999)];%memristor right shift

%error analysis
New_W_error_ij = abs(W_ij_shift - wij);
New_W_max_error_ij = max(New_W_error_ij);
New_W_mean_error_ij = mean(New_W_error_ij);
corr_pij = corrcoef(wij,W_ij_shift);
New_W_corr_ij = corr_pij(2,1);
rrmse_wij= sqrt((sum((W_ij_shift - wij).^2)/(sum(wij.^2)))/5000);

%% matlab comparison

% matlab data error analysis
New_mean_error_z_avg_i=mean(New_Z_mean_error_i);
New_max_error_z_avg_i=mean(New_Z_max_error_i);
New_corr_z_avg_i=mean(New_Z_corr_i);
rmse_zi=sqrt(mean((New_Z_i_shift - zi).^2));

New_mean_error_p_avg_i=mean(New_P_mean_error_i);
New_max_error_p_avg_i=mean(New_P_max_error_i);
New_corr_p_avg_i=mean(New_P_corr_i);
rmse_pi=sqrt(mean((New_P_i_shift - pi).^2));

New_mean_error_z_avg_j=mean(New_Z_mean_error_j);
New_max_error_z_avg_j=mean(New_Z_max_error_j);
New_corr_z_avg_j=mean(New_Z_corr_j);
rmse_zj=sqrt(mean((New_Z_j_shift - zj).^2));

New_mean_error_p_avg_j=mean(New_P_mean_error_j);
New_max_error_p_avg_j=mean(New_P_max_error_j);
New_corr_p_avg_j=mean(New_P_corr_j);
rmse_pj=sqrt(mean((New_P_j_shift - pj).^2));

New_mean_error_p_avg_ij=mean(New_P_mean_error_ij);
New_max_error_p_avg_ij=mean(New_P_max_error_ij);
New_corr_p_avg_ij=mean(New_P_corr_ij);
rmse_pij=sqrt(mean((New_P_ij_shift - pij).^2));

New_mean_error_w_avg_ij=mean(New_W_mean_error_ij);
New_max_error_w_avg_ij=mean(New_W_max_error_ij);
New_corr_w_avg_ij=mean(New_W_corr_ij);
rmse_wij=sqrt(mean((W_ij_shift - wij).^2));

New_mean_error_b_avg_j=mean(New_b_mean_error_j);
New_max_error_b_avg_j=mean(New_b_max_error_j);
New_corr_b_avg_j=mean(New_b_corr_j);
rmse_bj=sqrt(mean((b_j_shift - bj).^2));

fprintf('M Zi mean error: %f, Zi max error: %f, Zi corr: %f, Zi rrmse: %f, Zi rmse: %f\n',New_Z_mean_error_i,New_Z_max_error_i,New_Z_corr_i,rrmse_zi,rmse_zi);
fprintf('M Pi mean error: %f, Pi max error: %f, Pi corr: %f, Pi rrmse: %f, Pi rmse: %f\n',New_P_mean_error_i,New_P_max_error_i,New_P_corr_i,rrmse_pi,rmse_pi);

fprintf('M Zj mean error: %f, Zj max error: %f, Zj corr: %f, Zj rrmse: %f, Zj rmse: %f\n',New_Z_mean_error_j,New_Z_max_error_j,New_Z_corr_j,rrmse_zj,rmse_zj);
fprintf('M Pj mean error: %f, Pj max error: %f, Pj corr: %f, Pj rrmse: %f, Pj rmse: %f\n',New_P_mean_error_j,New_P_max_error_j,New_P_corr_j,rrmse_pj,rmse_pj);

fprintf('M Pij mean error: %f, Pij max error: %f, Pij corr: %f, Pij rrmse: %f, Pij rmse: %f\n',New_P_mean_error_ij,New_P_max_error_ij,New_P_corr_ij,rrmse_pij,rmse_pij);
fprintf('M wij mean error: %f, wij max error: %f, wij corr: %f, wij rrmse: %f, wij rmse: %f\n',New_W_mean_error_ij,New_W_max_error_ij,New_W_corr_ij,rrmse_wij,rmse_wij);
fprintf('M bj mean error: %f, bj max error: %f, bj corr: %f, bj rrmse: %f, bj rmse: %f\n',New_Z_max_error_i,New_b_max_error_j,New_b_corr_j,rrmse_bj,rmse_bj);

% color setting
blue1 = [101/255 158/255 206/255];
blue2 = [46/255 117/255 181/255];
green1 = [169/255 209/255 142/255];
green2 = [83/255 129/255 53/255];
orange1 = [244/255 177/255 131/255];
orange2 = [197/255 90/255 17/255];

% Matlab comparison of separative signal
figure(1)
subplot(511),
x=1:points;
plot(x*0.001,si,'Color',blue1, 'linewidth',0.8);hold on
plot(x*0.001,sj,'Color',green1, 'linewidth',0.8);
title('Matlab comparison of separative signal')

subplot(512),
plot(x*0.001,zi,'Color',blue2, 'linewidth',1.5);hold on
plot(x*0.001,New_Z_i_shift,'Color',blue1, 'linewidth',1.5);hold on
plot(x*0.001,zj,'Color',green2, 'linewidth',1.5);hold on
plot(x*0.001,New_Z_j_shift,'Color',green1, 'linewidth',1.5);

subplot(513);
plot(x*0.001,pi,'Color',blue2, 'linewidth',1.5);hold on
plot(x*0.001,New_P_i_shift,'Color',blue1, 'linewidth',1.5);hold on
plot(x*0.001,pj,'Color',green2, 'linewidth',1.5);hold on
plot(x*0.001,New_P_j_shift,'Color',green1, 'linewidth',1.5);hold on
plot(x*0.001,pij,'Color',orange2, 'linewidth',1.5);hold on
plot(x*0.001,New_P_ij_shift,'Color',orange1, 'linewidth',1.5);

subplot(514);
plot(x*0.001,wij,'Color',orange2, 'linewidth',1.5);hold on
plot(x*0.001,W_ij_shift,'Color',orange1, 'linewidth',1.5);

subplot(515);
plot(x*0.001,bj,'Color',green2, 'linewidth',1.5);hold on
plot(x*0.001,b_j_shift,'Color',green1, 'linewidth',1.5);


%% cadence comparison

% cadence data error analysis
Ca_Z_error_i = abs(zi_modified_sep - zi);
Ca_max_error_z_i = max(Ca_Z_error_i);
Ca_mean_error_z_i = mean(Ca_Z_error_i);
Ca_corr_zi = corrcoef(zi,zi_modified_sep);
Ca_corr_z_i = Ca_corr_zi(2,1);
Ca_rrmse_zi=sqrt((sum((zi_modified_sep - zi).^2)/(sum(zi.^2)))/5000);
Ca_rmse_zi=sqrt(mean((zi_modified_sep - zi).^2));

Ca_Z_error_j = abs(zj_modified_sep - zj);
Ca_max_error_z_j = max(Ca_Z_error_j);
Ca_mean_error_z_j = mean(Ca_Z_error_j);
Ca_corr_zj = corrcoef(zj,zj_modified_sep);
Ca_corr_z_j = Ca_corr_zj(2,1);
Ca_rrmse_zj=sqrt((sum((zj_modified_sep - zj).^2)/(sum(zj.^2)))/5000);
Ca_rmse_zj=sqrt(mean((zj_modified_sep - zj).^2));

Ca_P_error_i = abs(pi_modified_sep - pi);
Ca_max_error_p_i = max(Ca_P_error_i);
Ca_mean_error_p_i = mean(Ca_P_error_i);
Ca_corr_pi = corrcoef(pi,pi_modified_sep);
Ca_corr_p_i = Ca_corr_pi(2,1);
Ca_rrmse_pi=sqrt((sum((pi_modified_sep - pi).^2)/(sum(pi.^2)))/5000);
Ca_rmse_pi=sqrt(mean((pi_modified_sep - pi).^2));

Ca_P_error_j = abs(pj_modified_sep - pj);
Ca_max_error_p_j = max(Ca_P_error_j);
Ca_mean_error_p_j = mean(Ca_P_error_j);
Ca_corr_pj = corrcoef(pj,pj_modified_sep);
Ca_corr_p_j = Ca_corr_pj(2,1);
Ca_rrmse_pj=sqrt((sum((pj_modified_sep - pj).^2)/(sum(pj.^2)))/5000);
Ca_rmse_pj=sqrt(mean((pj_modified_sep - pj).^2));

Ca_P_error_ij = abs(pij_modified_sep - pij);
Ca_max_error_p_ij = max(Ca_P_error_ij);
Ca_mean_error_p_ij = mean(Ca_P_error_ij);
Ca_corr_pij = corrcoef(pij,pij_modified_sep);
Ca_corr_p_ij = Ca_corr_pij(2,1);
Ca_rrmse_pij= sqrt((sum((pij_modified_sep - pij).^2)/(sum(pij.^2)))/5000);
Ca_rmse_pij=sqrt(mean((pij_modified_sep - pij).^2));

Ca_W_error_ij = abs(wij_modified_sep - wij);
Ca_max_error_w_ij = max(Ca_W_error_ij);
Ca_mean_error_w_ij = mean(Ca_W_error_ij);
Ca_corr_wij = corrcoef(wij,wij_modified_sep);
Ca_corr_w_ij = Ca_corr_wij(2,1);
Ca_rrmse_wij= sqrt((sum((wij_modified_sep - wij).^2)/(sum(wij.^2)))/5000);
Ca_rmse_wij=sqrt(mean((wij_modified_sep- wij).^2));

logpj_modified_sep = logpj_modified_sep-4.60517;
Ca_b_error_j = abs(logpj_modified_sep - bj);
Ca_max_error_b_j = max(Ca_b_error_j);
Ca_mean_error_b_j = mean(Ca_b_error_j);
Ca_corr_bj = corrcoef(bj,logpj_modified_sep);
Ca_corr_b_j = Ca_corr_bj(2,1);
Ca_rrmse_bj= sqrt((sum((logpj_modified_sep - bj).^2)/(sum(bj.^2)))/5000);
Ca_rmse_bj=sqrt(mean((logpj_modified_sep - bj).^2));


fprintf('Ca Zi mean error: %f, Zi max error: %f, Zi corr: %f, Zi rrmse: %f, Zi rmse: %f\n',Ca_mean_error_z_i,Ca_max_error_z_i,Ca_corr_z_i,Ca_rrmse_zi,Ca_rmse_zi);
fprintf('Ca Pi mean error: %f, Pi max error: %f, Pi corr: %f, Pi rrmse: %f, Pi rmse: %f\n',Ca_mean_error_p_i,Ca_max_error_p_i,Ca_corr_p_i,Ca_rrmse_pi,Ca_rmse_pi);

fprintf('Ca Zj mean error: %f, Zj max error: %f, Zj corr: %f, Zj rrmse: %f, Zj rmse: %f\n',Ca_mean_error_z_j,Ca_max_error_z_j,Ca_corr_z_j,Ca_rrmse_zj,Ca_rmse_zj);
fprintf('Ca Pj mean error: %f, Pj max error: %f, Pj corr: %f, Pj rrmse: %f, Pj rmse: %f\n',Ca_mean_error_p_j,Ca_max_error_p_j,Ca_corr_p_j,Ca_rrmse_pj,Ca_rmse_pj);

fprintf('Ca Pij mean error: %f, Pij max error: %f, Pij corr: %f, Pij rrmse: %f, Pij rmse: %f\n',Ca_mean_error_p_ij,Ca_max_error_p_ij,Ca_corr_p_ij,Ca_rrmse_pij,Ca_rmse_pij);
fprintf('Ca wij mean error: %f, wij max error: %f, wij corr: %f, wij rrmse: %f, wij rmse: %f\n',Ca_mean_error_w_ij,Ca_max_error_w_ij,Ca_corr_w_ij,Ca_rrmse_wij,Ca_rmse_wij);
fprintf('Ca bj mean error: %f, bj max error: %f, bj corr: %f, bj rrmse: %f, bj rmse: %f\n',Ca_mean_error_b_j,Ca_max_error_b_j,Ca_corr_b_j,Ca_rrmse_bj,Ca_rmse_bj);

% Cadence comparison of separative signal
figure(2)
subplot(511),
x=1:points;
plot(x*0.001,si,'Color',blue1, 'linewidth',0.8);hold on
plot(x*0.001,sj,'Color',green1, 'linewidth',0.8);
title('Cadence comparison of separative signal')
subplot(512),
plot(x*0.001,zi,'Color',blue2, 'linewidth',1.5);hold on
plot(x*0.001,zi_modified_sep,'Color',blue1, 'linewidth',1.5);hold on
plot(x*0.001,zj,'Color',green2, 'linewidth',1.5);hold on
plot(x*0.001,zj_modified_sep,'Color',green1, 'linewidth',1.5);

subplot(513);
plot(x*0.001,pi,'Color',blue2, 'linewidth',1.5);hold on
plot(x*0.001,pi_modified_sep,'Color',blue1, 'linewidth',1.5);hold on
plot(x*0.001,pj,'Color',green2, 'linewidth',1.5);hold on
plot(x*0.001,pj_modified_sep,'Color',green1, 'linewidth',1.5);hold on
plot(x*0.001,pij,'Color',orange2, 'linewidth',1.5);hold on
plot(x*0.001,pij_modified_sep,'Color',orange1, 'linewidth',1.5);

subplot(514);
plot(x*0.001,wij,'Color',orange2, 'linewidth',1.5);hold on
plot(x*0.001,wij_modified_sep,'Color',orange1, 'linewidth',1.5);

subplot(515);
plot(x*0.001,bj,'Color',green2, 'linewidth',1.5);hold on
plot(x*0.001,logpj_modified_sep,'Color',green1, 'linewidth',1.5);

