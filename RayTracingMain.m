%% Ray Tracing Engine Parameters    
tic
optimizationMode = 0;   % if 1, then reduces the Rx mesh size to measurement locations only
plotMode = 1;           % Ray Tracing engine plot mode
demoMode = 1;           % Ray Tracing Engine Plot Mode

losFlag = 1; 
reflectionFlag = 1;                     % whether or not calculate First reflections
secondReflectionFlag = 1;
reflectExaggerationFac = 1e0;          % Must be 1 unless to emphasize reflection for demonstration purposes
                                       % whether or not calculate LoS
disableIncidentAngle = 0;               % 1 Disables the incident angle calculation, if disableIncidentAngle= 1
solidIncidentAngle = 45;                % if disableIncidentAngle =1, then assign this which overwrites all the incident angles! This is unnecessary feature
polarizationSwap = 1;               % (See notes in HOW THIS WORSK)  % 1, Applies TE to walls and TM to ceiling. 0, applies TM to the walls and TE to the ceiling

imageRSSIScale = 5;         % increase this if number of meshes nodes are small
imageEFScale = 5;         % increase this if number of meshes nodes are small
grayScaleImage = 0;

refDistance = 0.04;            % Reference distance from Tx
FPSLRefLoss = 0;

%use 40
antennaGainRes = 40;
antennaEffiLoss = -11.5;         % dB antenna efficiency, missmatch and loss all together

ceilingEnable = 0; % Allowing to define ceiling and floor
groundLevel = 0;
ceilingLevel = 1000;  % Height of the ceiling

% For LOS dis = length of ceiling
%timeLOS = ceilingLevel/vel;

% use 100
mesh_.xNodeNum = 10;   % Keep the x and y mesh size the same, increase the size for better resolution and especially if you're increasing the frequency
mesh_.yNodeNum = 10;
mesh_.zNodeNum = 1;

%% Antenna Gain pattern calculation

[TxAntennaGainAE] = AntennaTemp(antennaGainRes,demoMode) + antennaEffiLoss;  % TxAntennaGainAE needs to be in dB
RxAntennaGainAE = TxAntennaGainAE;

% Tx.xyz = [ 200, 200, 1]; % Location of the transmitter (s)
% Tx.xyz = [0.03, 0.03, 0]; % Location of the transmitter (s)
Tx.xyz = [45, 45, 0]; % Location of the transmitter (s)
Tx.power =  67; % power of the transmitter dB(m) % start from 47 at 2.4 GHz
% start from 2 at 60 GHz and decrease by 20 for duct
% start from 2 at 60 GHz and decrease by 2 for duct

% Defining the boundary of the analysis (something like a boundary condition) 
boundary = [0,90
            0,90
            0,1000];    
% boundary = [0,125.0
%             0,125.0
%             0,1000]; 

%% CLOCK WISE WALL DEFINITION        
[wallxyz1, wallxyz2, wallxyz3, wallxyz4,wallX,wallY,wallZ] = CSV23D_V1(demoMode,groundLevel,ceilingLevel,Tx.xyz);

wall.xyz1 = wallxyz1;
wall.xyz2 = wallxyz2;
wall.xyz3 = wallxyz3;
wall.xyz4 = wallxyz4;

wall.X = wallX;
wall.Y = wallY;
wall.Z = wallZ;

if ceilingEnable == 1
    ceillFloor.xyz1 = [0,0,ceilingLevel
                       0,0,groundLevel
                        ];

    ceillFloor.xyz2 = [0,0.1250,ceilingLevel
                       0,0.1250,groundLevel
                        ];

    ceillFloor.xyz3 = [0.1250,0.1250,ceilingLevel
                        0.1250,0.1250,groundLevel
                        ];

    ceillFloor.xyz4 = [0.1250,0,ceilingLevel
                        0.1250,0,groundLevel
                        ];
else
    ceillFloor.xyz1 = [];
    ceillFloor.xyz2 = [];
    ceillFloor.xyz3 = [];
    ceillFloor.xyz4 = [];
end
wall.relativePerm = 10*ones(size(wall.xyz1,1)+size(ceillFloor.xyz1,1),1);
%% Adding Ceillilng and Floor to the structure
for i = 1:size(ceillFloor.xyz1,1) 
    wall.xyz1 = [wall.xyz1;ceillFloor.xyz1(i,:)];
    wall.xyz2 = [wall.xyz2;ceillFloor.xyz2(i,:)];
    wall.xyz3 = [wall.xyz3;ceillFloor.xyz3(i,:)];
    wall.xyz4 = [wall.xyz4;ceillFloor.xyz4(i,:)];
end

frequency = 60.2e9;  % frequency in hertz
% frequency = (2.4:0.3125:2.5).*1e9;  % frequency in hertz 0.0003125
% frequency = 2.5*1e9;  
lightVel = 3e8;
lambda = lightVel./frequency;
velocity = frequency.*lambda;
% RayTracingEng_V01
alphaLOS = 1;
alphaFirstRef = 0.75;
alphaSecondRef = 0.5;
for v = 1:length(velocity)
    freq = frequency(v);
    lambda = lightVel./freq;
    vel = freq.*lambda;
    RayTracingEng_V01
%     for i = 1:size(Rx.xyz,1)
%         losDistance(i) = sqrt((Tx.xyz(1,1) - Rx.xyz(i,1))^2 + (Tx.xyz(1,2) - Rx.xyz(i,2))^2 +...
%             (Tx.xyz(1,3) - Rx.xyz(i,3))^2);
%         losTime(i) = losDistance(i)/velocity(v);
%         losPhase(i) = exp(-1i*pi*losTime(i));
%         losAmp(i) = alphaLOS/losDistance(i);
%         impRespLOS(i) = losAmp(i) * losPhase(i);
%         firstRefDistance(i) = sqrt((TxReflection1Rx(1,1,i) - TxReflection1Rx(1,4,i))^2 + ...
%             (TxReflection1Rx(1,2,i) - TxReflection1Rx(1,5,i))^2 +...
%             (TxReflection1Rx(1,3,i) - TxReflection1Rx(1,6,i))^2) + ...
%             sqrt((TxReflection1Rx(1,4,i) - TxReflection1Rx(1,7,i))^2 + ...
%             (TxReflection1Rx(1,5,i) - TxReflection1Rx(1,8,i))^2 +...
%             (TxReflection1Rx(1,6,i) - TxReflection1Rx(1,9,i))^2);
%         firstRefTime(i) = firstRefDistance(i)/velocity(v);
%         firstRefPhase(i) = exp(-1i*pi*firstRefTime(i));
%         firstRefAmp(i) = -alphaFirstRef/firstRefDistance(i);
%         impRespFirstRef(i) = firstRefAmp(i) * firstRefPhase(i);
%         secondRefDistance(i) = sqrt((TxReflection22Rx(1,1,i) - TxReflection22Rx(1,4,i))^2 + ...
%             (TxReflection22Rx(1,2,i) - TxReflection22Rx(1,5,i))^2 +...
%             (TxReflection22Rx(1,3,i) - TxReflection22Rx(1,6,i))^2) + ...
%             sqrt((TxReflection22Rx(1,4,i) - TxReflection22Rx(1,7,i))^2 + ...
%             (TxReflection22Rx(1,5,i) - TxReflection22Rx(1,8,i))^2 +...
%             (TxReflection22Rx(1,6,i) - TxReflection22Rx(1,9,i))^2) + ...
%             sqrt((TxReflection22Rx(1,7,i) - TxReflection22Rx(1,10,i))^2 + ...
%             (TxReflection22Rx(1,8,i) - TxReflection22Rx(1,11,i))^2 +...
%             (TxReflection22Rx(1,9,i) - TxReflection22Rx(1,12,i))^2);
%         secondRefTime(i) = secondRefDistance(i)/velocity(v);
%         secondRefPhase(i) = exp(-1i*pi*secondRefTime(i));
%         secondRefAmp(i) = alphaSecondRef/secondRefDistance(i);
%         impRespSecondRef(i) = secondRefAmp(i) * secondRefPhase(i);
%     end
%     impResp(v) = sum(impRespLOS) + sum(impRespFirstRef) + sum(impRespSecondRef);
end
PropagationMaps
PropagationMapsEF
toc

% L = [1 2 4 8];
% RwP = [-34.1212 -37.7700 -42.4794 -47.2173];
% RP = [-36.8592 -40.5080 -45.2173 -50.9552];
% RMP = [-94.0744 -97.6685 -101.3379 -105.0773];
% PP1 = [-82 -83 -82  -83];
% PP2 = [-80 -80 -80 -80];
% fspl = [-54 -60 -66 -72];
% figure 
% plot(L, RwP, 'r', L, RP, 'y',L, PP1, 'g',L, PP2, 'b',L, fspl, 'k');


length = [1 2 4 8];
pressure = [17.05 23.38 31.76 42.39];
temp = [288.15 293.15 298.15 303.15];
densityT1P1 = [3.2075 4.4905 6.415];
densityT2P2 = [4.325 6.055 8.65];
densityT3P3 = [5.75 8.05 11.5];
densityT4P4 = [7.6 10.64 15.2];

RSSID1T1P1 = [129.46 134.4 141.03 147.03];
RSSID2T1P1 = [130.43 135.4 142.03 147.96];
RSSID3T1P1 = [131.5 136.43 143.03 149.03];

RSSID1T2P2 = [132.59 137.53 144.16 150.16];
RSSID2T2P2 = [133.56 138.53 145.16 151.1];
RSSID3T2P2 = [134.56 139.53 146.16 152.16];

RSSID1T3P3 = [135.2 140.6 147.2 153.2];
RSSID2T3P3 = [136.2 141.6 148.2 154.2];
RSSID3T3P3 = [137.2 142.6 149.2 155.2];

RSSID1T4P4 = [138.8 143.72 150.36 156.32];
RSSID2T4P4 = [139.8 144.72 151.36 157.32];
RSSID3T4P4 = [140.8 145.72 152.36 158.32];

figure
subplot(2,2,1)
plot(length,RSSID1T1P1,length,RSSID2T1P1,length,RSSID3T1P1)
title('Temp = 288.15 and Pressure = 17.05');
legend('density = 3.2075','density = 4.4905','density = 6.415');  

subplot(2,2,2)
plot(length,RSSID1T2P2,length,RSSID2T2P2,length,RSSID3T2P2)
title('Temp = 293.15 and Pressure = 23.38');    
legend('density = 4.325','density = 6.055','density = 8.65'); 

subplot(2,2,3)
plot(length,RSSID1T3P3,length,RSSID2T3P3,length,RSSID3T3P3)
title('Temp = 298.15 and Pressure = 31.76');
legend('density = 5.75','density = 8.05','density = 11.5'); 

subplot(2,2,4)
plot(length,RSSID1T4P4,length,RSSID2T4P4,length,RSSID3T4P4)
title('Temp = 303.15 and Pressure = 42.39');
legend('density = 7.6','density = 10.64','density = 15.2'); 

RayT1 = [-29.6 -34.4 -40.9 -46.97];
RayT2 = [-29.1 -33.9 -40.5 -46.43];
RayT3 = [-28.5 -33.43 -40.06 -46.03];
RayT4 = [-27.9 -32.8 -39.4 -45.5];
RayT0 = [-24.43	-28.6 -36.13 -39.06];
FSPL = [-48.03 -54.05 -60.07 -66.09];

figure
plot(length,RayT1,length,RayT2,length,RayT3,length,RayT4,length,RayT0,length,FSPL)
legend('T=288.15','T=293.15','T=298.15','T=303.15','Without T','FSPL'); 

RayP1 = [-38.4 -43.4 -49.93 -55.96];
RayP2 = [-37.9 -42.83 -49.46 -55.43];
RayP3 = [-37.2 -42.1 -48.84 -54.87];
RayP4 = [-36.5 -41.5 -48.1 -54.13];
RayP0 = [-24.43	-28.6 -36.13 -39.06];
FSPL = [-48.03 -54.05 -60.07 -66.09];

figure
plot(length,RayP1,length,RayP2,length,RayP3,length,RayP4,length,RayP0,length,FSPL)
legend('D=1200','D=1190','D=1165','D=1140','Without moisture vapour','FSPL'); 

figure
x = ['FSPL','Duct w/o air','Duct with dry air','Duct with moist air',...
    'Max. ray beamforming [10]','With Channel Model [11]'];
y = [-66.09, -39.06, -46.43, -55.43, -80, -70];
b = bar(y,0.4);
xtips1 = ['1','2','3','4','5','6']; %b.XEndPoints;
ytips1 = ['FSPL','Duct w/o air','Duct with dry air',...
    'Duct with moist air','Max. ray beamforming [10]','With Channel Model [11]']; %b.YEndPoints;
% labels1 = string(b(1).YData);
% text(xtips1, ytips1)


freq = 59.4E9:1E7:60.4E9;
rayTRssi(1,:) = [-34.03*ones(1,20) -33.45*ones(1,20) -34.76*ones(1,20) -33.63*ones(1,20) ...
    -34.43*ones(1,20) -33.73];
rayTRssi(2,:) = [-35.8*ones(1,20) -37.64*ones(1,20) -37.3*ones(1,20) -36.43*ones(1,20) ...
    -35.6*ones(1,20) -35.73];
rayTRssi(3,:) = [-41.9*ones(1,20) -43*ones(1,20) -42.9*ones(1,20) -42.43*ones(1,20) ...
    -43.13*ones(1,20) -43.73];
rayTRssi(4,:) = [-48.43*ones(1,20) -47.53*ones(1,20) -48.96*ones(1,20) -48.93*ones(1,20) ...
    -47.16*ones(1,20) -48.73];
expRSSI(1,:) = -35*ones(1,101);
expRSSI(2,:) = -38*ones(1,101);
expRSSI(3,:) = -44*ones(1,101);
expRSSI(4,:) = -49*ones(1,101);
figure
subplot(4,1,1)
plot(freq/1E9,expRSSI(1,:),'b',freq/1E9,rayTRssi(1,:),'r')
subplot(4,1,2)
plot(freq/1E9,expRSSI(2,:),'b',freq/1E9,rayTRssi(2,:),'r')
subplot(4,1,3)
plot(freq/1E9,expRSSI(3,:),'b',freq/1E9,rayTRssi(3,:),'r')
subplot(4,1,4)
plot(freq/1E9,expRSSI(4,:),'b',freq/1E9,rayTRssi(4,:),'r')
