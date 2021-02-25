function [radresTE, radresTM, gammaTE, gammaTM] = radResBend(m_TE,n_TE,m_TM,n_TM,radius,freq,c,k,R,eta,l,coWnTE,coWnTM,Rb)
for fr=1:length(freq)
    for n = 1:length(n_TE)
        gTE = coWnTE(n);
        gbTE = gTE*sqrt(1 + (3/((2*Rb*gTE)^2)));
        gTESq = gbTE^2;
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
        gTM = coWnTM(n);
        gbTM = gTM*sqrt(1 + (3/((2*Rb*gTM)^2)));
        gTMSq = gbTM^2;
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
end