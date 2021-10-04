%% Pt-Hf-Ti memristor parameter
w_init1 = 1; % the initial state condition [0:1] 
D1 = 10;
P=1;
j = 1;
Roff1 = 2.5e3;
Ron1 = 100;
alpha_on1 = 3;
alpha_off1 = 1;
k_on1 = -8e10;
k_off1 = 40.3;
v_on1=-0.53;
v_off1=0.5;
delta_t = 1e-3;
points=200;
% For the Li window function:
Li_X=zeros(1,200);
Li_X_dot=zeros(1,200);
Li_X(1)=w_init1*D1;
J_Li=1;
alpha_Li=0;
a_Li=1;
beta_Li=-0.3;
gama_Li=0.3;
P_Li=1;
%% Pt-Hf-Ti memristor input
v1=zeros(1,200);
v1(1:200)=-0.53119;
v1(101:4:200)=50;
%% Pt-Hf-Ti memristor
for p=2:200
    if  (v1(p) > 0) && (v1(p) > v_off1)
        Li_X_dot(p)=k_off1*(v1(p)/v_off1-1)^alpha_off1;
        Li_X(p)=Li_X(p-1)+delta_t*Li_X_dot(p).*(J_Li*(1-(alpha_Li*(Li_X(p-1)/D1)^3 + (a_Li^2)*(Li_X(p-1)/D1)^2 + (1-a_Li^2) + beta_Li*(Li_X(p-1)/D1)^2 + gama_Li*(Li_X(p-1)/D1))^P_Li));
    elseif (v1(p) <= 0) && (v1(p) < v_on1)
        Li_X_dot(p)=k_on1*((v1(p))/v_on1-1)^alpha_on1;
        Li_X(p)=Li_X(p-1)+delta_t*Li_X_dot(p).*(J_Li*(1-(alpha_Li*(Li_X(p-1)/D1)^3 + (a_Li^2)*(Li_X(p-1)/D1 -1)^2 + (1-a_Li^2) + beta_Li*(Li_X(p-1)/D1)^2 + gama_Li*(Li_X(p-1)/D1))^P_Li));
    else
        Li_X(p)=Li_X(p-1);
        Li_X_dot(p)=0;
    end  
end
Li_R=Roff1.*Li_X./D1+Ron1.*(1-Li_X./D1);
%% Ferroelectric memristor parameter
w_init2 = 0; % the initial state condition [0:1] 
D2 = 10;
Roff2 = 5e7;
Ron2 = 1.5e5;
alpha_on2 = 5;
alpha_off2 = 5;
k_on2 = -3e10;
k_off2 = 1e5;
v_on2=-5.7;
v_off2=1.4;
delta_t=1e-3;
% For the New window function:
v2=zeros(1,200);
X2=zeros(1,200);
X_dot2=zeros(1,200);
X2(1)=w_init2*D2;
%% Ferroelectric memristor input
v2(1:200)=-586.5e-2;
v2(1:4:100)=200e-2;
v2 = [-587e-2,v2(1:199)]; 
%% Ferroelectric memristor
for i=2:200
    if (v2(i) > 0) && (v2(i) > v_off2)

        X_dot2=k_off2*(v2(i)/v_off2-1)^alpha_off2*j*(1-(X2(i-1)/D2)^(2*P));
        X2(i)=X2(i-1)+delta_t*X_dot2;
    elseif (v2(i) <= 0) && (v2(i) < v_on2)

        X_dot2(i)=k_on2*(v2(i)/v_on2-1)^alpha_on2*j*(1-(X2(i-1)/D2-1)^(2*P));
        X2(i)=X2(i-1)+delta_t*X_dot2(i);
    else

        X2(i)=X2(i-1);
        X_dot2(i)=0;
    end
end
R2=Roff2.*X2./D2+Ron2.*(1-X2./D2);
R2=[1.5e5,R2(1:199)];
%% figure
blue3 = [101/255 158/255 206/255];
blue2 = [46/255 117/255 181/255];
blue3 = [86/255 224/255 224/255];
green1 = [169/255 209/255 142/255];
green2 = [83/255 129/255 53/255];
green3 = [69/255 239/255 61/255];
orange1 = [244/255 177/255 131/255];
orange2 = [197/255 90/255 17/255];
orange3 = [180/255 53/255 31/255];
gray = [165/255 165/255 165/255];

t=1:200;

figure(1);
plot(t(1:100)*delta_t,Li_R(1:100),'Color',green1, 'linewidth',2);hold on
plot(t(100:200)*delta_t,Li_R(100:200),'Color',green2, 'linewidth',2);hold on
title('The resistance of the Pt-Hf-Ti memristor')
hold off

figure(2);
plot(t(1:99)*delta_t,R2(1:99),'Color',blue2, 'linewidth',2);hold on
plot(t(99:200)*delta_t,R2(99:200),'Color',blue1, 'linewidth',2);
set(gca,'Color','white') 
title('The resistance of the Ferroelectric memristor')
hold off

figure(3);
plot(t(1:100)*delta_t,v1(1:100),'Color',green1,'LineWidth',2);hold on
plot(t(100:200)*delta_t,v1(100:200),'Color',gray,'LineWidth',2);
title('input voltage of the Pt-Hf-Ti memristor')
hold off

figure(4);
plot(t(1:99)*delta_t,v2(1:99),'Color',gray,'LineWidth',2);hold on
plot(t(99:200)*delta_t,v2(99:200),'Color',blue1,'LineWidth',2);
title('input voltage of the Ferroelectric memristor')
hold off
