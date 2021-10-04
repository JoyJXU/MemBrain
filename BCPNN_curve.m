%% Initial value of BCPNN
pi(1) = 0;
ei(1) = 0;
zi(1) = 0;
si(1) = 0;
pj(1) = 0;
ej(1) = 0;
zj(1) = 0;
sj(1) = 0;
pij(1) = 0;
eij(1) = 0;
pi1(1) = 0;
ei1(1) = 0;
zi1(1) = 0;
si1(1) = 0;
pj1(1) = 0;
ej1(1) = 0;
zj1(1) = 0;
sj1(1) = 0;
pij1(1) = 0;
eij1(1) = 0;
x(1) = 1;
%% BCPNN parameter
kp = 1/500;
ke = 1/60;
kz = 1/11;
kfti = 1/(11);
eps=0.01;
%% BCPNN input
for i= 2: 1000
    x(i) = i;
    si(i) = 0;
    sj(i) = 0;
    if ((i > 90) && (i < 260))
        if( i < 155)
        if (rand>0.2)
            if ((si(i-1)== 1)||(si(i-2)== 1)||(si(i-3)== 1))
                si(i) = 0;
            else
                si(i) = 1;
            end
        end
        end
        sj(i) = 0;
        if (i > 220)
            si(i) = 0;
            sj(i) = 0;
            if (rand>0.2)
                if ((sj(i-1)== 1)||(sj(i-2)== 1)||(sj(i-3)== 1))
                    sj(i) = 0;
                else
                    sj(i) = 1;
                end
            end
        end
    end
        if ((i > 500)&&(i < 700))
            sj(i) = 0;
            si(i) = 0;
            if (i > 540)
                if (rand>0.3)
                    if ((sj(i-1)== 1)||(sj(i-2)== 1)||(sj(i-3)== 1))
                        sj(i) = 0;
                    else
                        sj(i) = 1;
                    end
                end
            end
            if(i<650)
                if (rand>0.25)
                    if ((si(i-1)== 1)||(si(i-2)== 1)||(si(i-3)== 1))
                        si(i) = 0;
                    else
                        si(i) = 1;
                    end
                end
            end
        end
    % original BCPNN
    pi(i) = pi(i-1) * (1 - kp) + ei(i-1)*kp;
    ei(i) = ei(i-1) *( 1 - ke) + zi(i-1)*ke;
    zi(i) = zi(i-1) * (1 - kz) + si(i-1)*kfti;
    pj(i) = pj(i-1) * (1 - kp) + ej(i-1)*kp;
    ej(i) = ej(i-1) *( 1 - ke) + zj(i-1)*ke;
    zj(i) = zj(i-1) * (1 - kz) + sj(i-1)*kfti;
    pij(i) = pij(i-1) * (1 - kp) + eij(i-1)*kp;
    eij(i) = eij(i-1) * (1 - ke) + zi(i-1)*zj(i-1)*ke;
    wij(i) = log((pij(i-1)+eps^2)/((pi(i-1)+eps)*(pj(i-1)+eps)));
    % simplified BCPNN
    pi1(i) = pi1(i-1) *( 1 - kp) + zi1(i-1)*kp;
    zi1(i) = zi1(i-1) * (1 - kz) + si(i-1)*kfti;
    pj1(i) = pj1(i-1) *( 1 - kp) + zj1(i-1)*kp;
    zj1(i) = zj1(i-1) * (1 - kz) + sj(i-1)*kfti;
    pij1(i) = pij1(i-1) * (1 - kp) + zi1(i-1)*zj1(i-1)*kp;
    wij1(i) = log((pij1(i-1)+eps^2)/((pi1(i-1)+eps)*(pj1(i-1)+eps)));
end
%% BCPNN figure
%color setting
blue1 = [101/255 158/255 206/255];
blue2 = [46/255 117/255 181/255];
green1 = [169/255 209/255 142/255];
green2 = [83/255 129/255 53/255];
orange1 = [244/255 177/255 131/255];
orange2 = [197/255 90/255 17/255];
%original BCPNN
figure(1)
subplot(511),
x=1:1000;
plot(x*0.001,si,'Color',blue1, 'linewidth',0.8);hold on
plot(x*0.001,sj,'Color',green1, 'linewidth',0.8);
title('original BCPNN')

subplot(512),
plot(x*0.001,zi,'Color',blue1, 'linewidth',1.5);hold on
plot(x*0.001,zj,'Color',green1, 'linewidth',1.5);

subplot(513);
plot(x*0.001,ei,'Color',blue1, 'linewidth',1.5);hold on
plot(x*0.001,ej,'Color',green1, 'linewidth',1.5);hold on
plot(x*0.001,eij,'Color',orange1, 'linewidth',1.5);

subplot(514);
plot(x*0.001,pi,'Color',blue1, 'linewidth',1.5);hold on
plot(x*0.001,pj,'Color',green1, 'linewidth',1.5);hold on
plot(x*0.001,pij,'Color',orange1, 'linewidth',1.5);

subplot(515);
plot(x*0.001,wij,'Color',orange1, 'linewidth',1.5);

% simplify BCPNN
figure(2)
subplot(511),
x=1:1000;
plot(x*0.001,si,'Color',blue1, 'linewidth',0.8);hold on
plot(x*0.001,sj,'Color',green1, 'linewidth',0.8);
title('simplify BCPNN')

subplot(512),
plot(x*0.001,zi1,'Color',blue1, 'linewidth',1.5);hold on
plot(x*0.001,zj1,'Color',green1, 'linewidth',1.5);

subplot(514);
plot(x*0.001,pi1,'Color',blue1, 'linewidth',1.5);hold on
plot(x*0.001,pj1,'Color',green1, 'linewidth',1.5);hold on
plot(x*0.001,pij1,'Color',orange1, 'linewidth',1.5);

subplot(515);
plot(x*0.001,wij1,'Color',orange1, 'linewidth',1.5);
