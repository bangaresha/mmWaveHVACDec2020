clear;
close all;
c=3*10^8;
% radius=0.305/2;
radius=0.127/2;
freq = 59.4E9:1E7:60.4E9;
omega = 2*pi*freq;

load('besDerZerMat5k.mat');
load('besZerMat5k.mat');

for fi=1:length(freq)
    n_TE=[];
    m_TE=[];
    fc_TE=[];
    coWnTE=[];
    for m=1:5001
        for n=1:5000
            tt = besDerZerMat(m,n);
            fc_TE_temp=(c/(2*pi*radius))*besDerZerMat(m,n);
            if fc_TE_temp <= freq(fi)
                n_TE = [n_TE n];
                m_TE = [m_TE m-1];
                fc_TE = [fc_TE fc_TE_temp];
                coWnTE = [coWnTE besDerZerMat(m,n)/radius];
            end
        end
    end
end

for fi=1:length(freq)
    n_TM=[];
    m_TM=[];
    fc_TM=[];
    coWnTM=[];
    for m=1:5001
        for n=1:5000
            fc_TM_temp=(c/(2*pi*radius))*besZerMat(m,n);
            if fc_TM_temp <= freq(fi)
                n_TM = [n_TM n];
                m_TM = [m_TM m-1];
                fc_TM = [fc_TM fc_TM_temp];
                coWnTM = [coWnTM (besZerMat(m,n))/radius];
            end
        end
    end
end

Zo=50;
eta = 377;
% l = 0.0306;
l = 0.0025;
sigma = 10^6;
mu = 4*pi*10^-7;
R = ((2*pi*freq*mu)./(2*sigma)).^0.5;
k = 2*pi*freq./c;

[radresTES, radresTMS, gammaTES, gammaTMS,p,fc] = radResCyl_multitone(m_TE,n_TE,m_TM,n_TM,...
    radius,freq,fc_TE,fc_TM,c,k,R,eta,l,coWnTE,coWnTM);
att = [];
sigRssi = [];
WGlenS= [1 2 4 8];
% WGlenS = 4;
for l = 1:length(WGlenS)
    [attT, sigRssiT] = chImpRespStr(freq,WGlenS(l),Zo,radresTES,radresTMS,gammaTES,gammaTMS);
    att = [att; attT];
    sigRssi = [sigRssi; sigRssiT];
end
% sigRssi = att-14;
channel = 10.^(att./10);

figure
plot(freq/1e9,att);
title('RSSI versus frequency Bend');

figure
plot(freq/1e9,channel);
title('RSSI versus frequency Bend');
% 
% figure
% plot(freq/1E9,att);
% title('Attenuation versus frequency');

figure
subplot(4,1,1)
plot(freq/1E9,sigRssi(1,:),'k-*')
subplot(4,1,2)
plot(freq/1E9,sigRssi(2,:),'r-o')
subplot(4,1,3)
plot(freq/1E9,sigRssi(3,:),'b-d')
subplot(4,1,4)
plot(freq/1E9,sigRssi(4,:),'g-s')
% title('RSSI versus frequency Bend');

rayTRssi(1,:) = [-27.03*ones(1,20) -26.45*ones(1,20) -27.76*ones(1,20) -26.63*ones(1,20) ...
    -27.43*ones(1,20) -26.73];
rayTRssi(2,:) = [-28.8*ones(1,20) -30.64*ones(1,20) -30.3*ones(1,20) -29.43*ones(1,20) ...
    -28.6*ones(1,20) -28.73];
rayTRssi(3,:) = [-34.9*ones(1,20) -36*ones(1,20) -35.9*ones(1,20) -35.43*ones(1,20) ...
    -36.13*ones(1,20) -36.73];
rayTRssi(4,:) = [-40.43*ones(1,20) -39.53*ones(1,20) -40.96*ones(1,20) -40.93*ones(1,20) ...
    -39.16*ones(1,20) -40.73];
figure
subplot(4,1,1)
plot(freq/1E9,sigRssi(1,:),'k',freq/1E9,rayTRssi(1,:),'r')
subplot(4,1,2)
plot(freq/1E9,sigRssi(2,:),'k',freq/1E9,rayTRssi(2,:),'r')
subplot(4,1,3)
plot(freq/1E9,sigRssi(3,:),'k',freq/1E9,rayTRssi(3,:),'r')
subplot(4,1,4)
plot(freq/1E9,sigRssi(4,:),'k',freq/1E9,rayTRssi(4,:),'r')
% title('RSSI versus frequency Bend');
% 
% figure
% plot(freq,att(1,:),'k-*',freq,att(2,:),'r-o',freq,att(3,:),'b-d',freq,att(4,:),'g-s');
% title('Attenuation versus frequency');


