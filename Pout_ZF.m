clear all
Mt = 2; % Tx antenna (Mt > 1)
Mr = 2; % Rx antenna (Mr >= Mt)
Me = 2; % eavesdropping antenna (Me >= Mt)
SNRr_dB = [0:3:30]; % average SNR at Rx (in dB)
SNRr = 10.^(SNRr_dB/10);
SNRe_dB = 10*ones(size(SNRr_dB)); % average SNR at eavesdropper (in dB)
SNRe = 10.^(SNRe_dB/10);
R = 1; % target secrecy rate of 1 bits/s/Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
po1 = zeros(size(SNRr)); % numerical integration for P_out,ZF
po2 = zeros(size(SNRr)); % exact closed-form expression (19) for P_out,ZF
for s = 1:length(SNRr)
    temp1 = 0;
    temp2 = 0;
    for k = 0:Mr-Mt
        temp1 = temp1 + integral(@(x) x.^(Me-Mt).*exp(-x/SNRe(s))./(factorial(Me-Mt)*SNRe(s)^(Me-Mt+1)).*...
                exp(-((2^R)*x+(2^R)-1)/SNRr(s)).*...
                1/(factorial(k)*SNRr(s)^k).*((2^R)*x+(2^R)-1).^k,0,inf);
        for l = 0:k
            temp2 = temp2 + 1/(factorial(k)*SNRr(s)^k)*nchoosek(k,l)*2^(l*R)*(2^R-1)^(k-l)*factorial(l+Me-Mt)/(1/SNRe(s)+2^R/SNRr(s))^l;
        end
    end
    po1(s) = 1 - temp1;
    po2(s) = 1 - exp(-(2^R-1)/SNRr(s))/(factorial(Me-Mt)*(1+2^R*SNRe(s)/SNRr(s))^(Me-Mt+1))*temp2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ps1 = zeros(size(SNRr)); % asymptotic numerical integration for P_out,ZF
ps2 = zeros(size(SNRr)); % asymptotic closed-form expression (28) for P_out,ZF
for s = 1:length(SNRr)
    temp = 0;
    ps1(s) = ps1(s) + integral(@(x) exp(-x/SNRe(s)).*x.^(Me-Mt)./(factorial(Me-Mt)*SNRe(s)^(Me-Mt+1)).*...
             ((2^R)*x+(2^R)-1).^(Mr-Mt+1)/(factorial(Mr-Mt+1)*SNRr(s)^(Mr-Mt+1)),0,inf);                 
    for m = 0:Mr-Mt+1
        temp = temp + nchoosek(Mr-Mt+1,m)*2^(m*R)*(2^R-1)^(Mr-Mt+1-m)*SNRe(s)^m*factorial(Me-Mt+m);    
    end
    ps2(s) = temp/(factorial(Me-Mt)*factorial(Mr-Mt+1)*SNRr(s)^(Mr-Mt+1));
end
semilogy(SNRr_dB,po1,SNRr_dB,po2,'o',SNRr_dB,ps1,SNRr_dB,ps2,'o','linewidth',1) % plot P_out vs SNR