function [attS, sigRssiS] = chImpRespStr(freq,WGlen,Zo,radresTES,radresTMS,gammaTES,gammaTMS)

radreacTES = imag(hilbert(radresTES));
radreacTMS = imag(hilbert(radresTMS));

for fi=1:length(freq)
    antresTE(fi) = sum(radresTES(fi,:));
    antresTM(fi) = sum(radresTMS(fi,:));
    antreacTE(fi) = sum(radreacTES(fi,:));
    antreacTM(fi) = sum(radreacTMS(fi,:));
    antimpTE(fi) = antresTE(fi) +1i*antreacTE(fi);
    antimpTM(fi) = antresTM(fi) +1i*antreacTM(fi);
    for n = 1:length(radresTES(fi,:))
        TEmodeimp(fi,n)=radresTES(fi,n) + 1i*radreacTES(fi,n);
    end
    for n = 1:length(radresTMS(fi,:))
        TMmodeimp(fi,n)=radresTMS(fi,n) + 1i*radreacTMS(fi,n);
    end
    TsTE = diag(exp(-1*gammaTES(fi,:)*WGlen));
    chTEmodeS = TEmodeimp(fi,:)*TsTE;
    TsTM = diag(exp(-1*gammaTMS(fi,:)*WGlen));
    chTMmodeS = TMmodeimp(fi,:)*TsTM;
    totPower(fi) = sum(sum(TsTE)) + sum(sum(TsTE));
    channelS(fi) = ((2*Zo)/(abs(antimpTM(fi)+Zo+antimpTE(fi))^2))*...
        (sum(chTEmodeS)+sum(chTMmodeS));
    if isnan(channelS(fi)) == 1
        channelS(fi) = 0;
    end
    
    attS(fi) = 10*log10(abs(channelS(fi)));
%     sigRssiS(fi) = attS(fi) - 14;
    
    p = 1013;
    T = 293.15;  %293 to 300
    rho = 1300; %1300
    rp = p/1013;
    rt = 288/T;
    % gammaY = ((rt*3.27E-2) + (1.67E-3 * (rho*(rt^7)/rp)) + (7.7E-4 * sqrt(freq)) +...
    %     (3.79/(((freq - 22.235)^2) + 9.81*rt*rp*rp)) + ((11.71*rt)/(((freq - 183.31)^2) + 11.85*rt*rp*rp)) + ...
    %     ((4.01*rt)/(((freq - 325.153)^2) + 10.44*rt*rp*rp))) * freq*freq*rho*rp*rt*1E-4;
    gammaD = 14.94*(rp^2)*(rt^8.5);
     gammaY = 1;
    % gammaD = 1;
    sigRssiS(fi) = attS(fi) -10*log10(gammaD)-10*log10(gammaY);
%     sigRssiS(fi) = attS(fi) - 14 -10*log10(gammaD)-10*log10(gammaY);
end
end