function [radresTE, radresTM, gammaTE, gammaTM,p,fc] = radResCyl_multitone(m_TE,n_TE,m_TM,n_TM,radius,...
    freq,fc_TE,fc_TM,c,k,R,eta,l,coWnTE,coWnTM)
gammaTE = zeros(length(freq), length(n_TE));
gammaTM = zeros(length(freq), length(n_TM));
radresTE = zeros(length(freq), length(n_TE));
radresTM = zeros(length(freq), length(n_TM));
 
for fr=1:length(freq)
    for n = 1:length(n_TE)
        gTESq = (coWnTE(n))^2;
        betaTE = (sqrt(k(fr)^2 - gTESq));%*1E-3;
        alphaTempTE = (R(fr)*k(fr))/(eta*radius*betaTE);
        alphaTE=8.85*(alphaTempTE*(((m_TE(n)^2)/...
            (((radius*(coWnTE(n)))^2)- ((m_TE(n))^2)))...
            +(((c*coWnTE(n))/(2*pi*freq(fr)))^2)));
        gammaTE(fr,n)=alphaTE+1i*betaTE;
        constRadTE = eta*(k(fr))*((m_TE(n))^2)/...
            (betaTE*pi*((sin(l*(k(fr))*pi/180))^2)*...
            ((besselj(m_TE(n),radius*coWnTE(n)))^4)*...
            ((((radius*coWnTE(n))^2)-((m_TE(n))^2))^2));
        radTEfunc=@(zeta) ((besselj(m_TE(n),radius.*zeta.*coWnTE(n)).*...
            sin((k(fr).*radius.*((l/radius)-1+zeta))*pi/180))./zeta);
        radTEnum = (integral(radTEfunc,1-(l/radius),1)).^2;
        radresTE(fr,n)= constRadTE*radTEnum;
        if imag(radresTE(fr,n)) ~= 0
            radresTE(fr,n) =0;
            gammaTE(fr,n) = 0;
        end
    end
    for n = 1:length(n_TM)
        gTMSq = (coWnTM(n))^2;
        betaTM = (sqrt(k(fr)^2 - gTMSq));%*1E-3;
        alphaconstTM = (R(fr)*k(fr))/(eta*radius*betaTM);
        alphaTM=8.85*alphaconstTM;
        gammaTM(fr,n)=alphaTM+1i*betaTM;
        constRadTM = eta*betaTM/(pi*((sin((k(fr))*l*pi/180))^2)*...
            (k(fr))*((0.5*((besselj(m_TM(n)-1,radius*coWnTM(n)))-...
            (besselj(m_TM(n)+1,radius*coWnTM(n)))))^2));
        
        radTMfunc=@(zeta) 0.5*((besselj((m_TM(n))-1,zeta.*...
            (coWnTM(n))*radius))-(besselj(m_TM(n)+1,zeta.*...
            (coWnTM(n))*radius))).*...
            sin((k(fr).*radius.*((l./radius)-1+zeta))*pi/180);
        radTMnum =(integral(radTMfunc,1-(l./radius),1)).^2;
        radresTM(fr,n)= constRadTM*radTMnum;
        if imag(radresTM(fr,n)) ~= 0
            radresTM(fr,n) =0;
            gammaTM(fr,n) = 0;
        end
    end
end
pTE_f = [];
fcTE_f = [];
pTM_f = [];
fcTM_f = [];
sumRadResTM = sum(radresTM(60,:));
sumRadResTE = sum(radresTE(60,:));
for i = 1:length(radresTE(60,:))
    pTE_f = [pTE_f; 100*radresTE(60,i)/sumRadResTE];
    fcTE_f = [fcTE_f; fc_TE(i)];
end
for i = 1:length(radresTM(60,:))
    pTM_f = [pTM_f; 100*radresTM(60,i)/sumRadResTM];
    fcTM_f = [fcTM_f; fc_TM(i)];
end
p = [pTE_f; pTM_f];
fc = [fcTE_f; fcTM_f];
sc = stem(fc,p);
title('Mode Power(%Total Power) Versus Mode Cutoff Frequency');
% modePower = [sum(radRes) counting];
% modePower = [p; counting];
end

% for fr=1:length(freq)
%     sumRadResTE(fr) = sum(radresTE(fr,:));
%     for i = 1:length(radresTE(fr,:))
%         pTE(fr,i) = 100*radresTE(fr,i)/sumRadResTE(fr);
% %         if pTE(fr,i) >= 10
%             pTE_f = [pTE_f; pTE(fr,i)];
%             fcTE_f = [fcTE_f; fc_TE(i)];
% %         end
%     end
% end
% 
% for fr=1:length(freq)
%     sumRadResTM(fr) = sum(radresTM(fr,:));
%     for i = 1:length(radresTM(fr,:))
%         pTM(fr,i) = 100*radresTM(fr,i)/sumRadResTM(fr);
% %         if pTM(fr,i) >= 10
%             pTM_f = [pTM_f; pTM(fr,i)];
%             fcTM_f = [fcTM_f; fc_TM(i)];
% %         end
%     end
% end

