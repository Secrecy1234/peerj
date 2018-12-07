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
po1 = zeros(size(SNRr)); % numerical integration for P_out,MMSE
po2 = zeros(size(SNRr)); % exact closed-form expression (20) for P_out,MMSE
for s = 1:length(SNRr)
    for m1 = 0:Me-1
        tempe = 0;
        for n1 = max([0 m1-Mt+1]):m1
            tempe = tempe + nchoosek(Mt-1,m1-n1)/(factorial(n1)*SNRe(s)^n1);
        end
        ae(m1+1) = tempe;
    end
    for m2 = 0:Mr-1
        tempr = 0;
        for n2 = max([0 m2-Mt+1]):m2
            tempr = tempr + nchoosek(Mt-1,m2-n2)/(factorial(n2)*SNRr(s)^n2);
        end
        ar(m2+1) = tempr;
    end
    temp = 0;
    for m1 = 0:Me-1
        for m2 = 0:Mr-1
            temp = temp + integral(@(x) exp(-x/SNRe(s))./((1+x).^Mt)*ae(m1+1).*x.^(m1-1).*(x.^2/SNRe(s)+(Mt+1/SNRe(s)-m1-1)*x-m1).*...
                   exp(-(2^R*x+2^R-1)/SNRr(s))./((2^R*x+2^R).^(Mt-1))*ar(m2+1).*(2^R*x+2^R-1).^m2,0,inf);   
        end
    end            
    po1(s) = 1 - temp;
    
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    for m = 0:Me-1
        for n = 0:Mr-1
            for k = 0:n
                for l1 = 0:m+1
                    temp1 = temp1 + ae(m+1)*ar(n+1)*nchoosek(n,k)*(-1)^k*2^((n-k)*R)*...
                            1/SNRe(s)*nchoosek(m+1,l1)*(-1)^l1*igamma(m+n-k-l1-2*Mt+3,2^R/SNRr(s)+1/SNRe(s))/(2^R/SNRr(s)+1/SNRe(s))^(m+n-k-l1-2*Mt+3);
                end
                for l2 = 0:m
                    temp2 = temp2 + ae(m+1)*ar(n+1)*nchoosek(n,k)*(-1)^k*2^((n-k)*R)*...
                            (Mt+1/SNRe(s)-m-1)*nchoosek(m,l2)*(-1)^l2*igamma(m+n-k-l2-2*Mt+2,2^R/SNRr(s)+1/SNRe(s))/(2^R/SNRr(s)+1/SNRe(s))^(m+n-k-l2-2*Mt+2);                       
                end
                for l3 = 0:m-1
                    temp3 = temp3 + ae(m+1)*ar(n+1)*nchoosek(n,k)*(-1)^k*2^((n-k)*R)*...
                            m*nchoosek(m-1,l3)*(-1)^l3*igamma(m+n-k-l3-2*Mt+1,2^R/SNRr(s)+1/SNRe(s))/(2^R/SNRr(s)+1/SNRe(s))^(m+n-k-l3-2*Mt+1);
                end
            end
        end
    end
    po2(s) = 1 - exp(1/SNRr(s)+1/SNRe(s))/2^((Mt-1)*R)*(temp1 + temp2 - temp3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ps1 = zeros(size(SNRr)); % asymptotic numerical integration for P_out,MMSE
ps2 = zeros(size(SNRr)); % asymptotic closed-form expression (32) for P_out,MMSE
for s = 1:length(SNRr)
    for m = 0:Me-1
        ps1(s) =  ps1(s) + integral(@(x) exp(-x/SNRe(s))./((x+1).^Mt)*ae(m+1).*(x.^(m+1)/SNRe(s)+(Mt+1/SNRe(s)-m-1).*(x.^m)-m*x.^(m-1)).*...
                  1/(factorial(Mr-Mt+1)*SNRr(s)^(Mr-Mt+1)).*((2^R*x+2^R-1).^Mr)./(2^R*x+2^R).^(Mt-1),0,inf);
    end
    
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    for m = 0:Me-1
        for n = 0:Mr
            for k1 = 0:m+1
                temp1 = temp1 + ae(m+1)*nchoosek(Mr,n)*(-1)^n*SNRe(s)^(m-n+Mr-2*Mt+1)/2^(n*R)*...
                        SNRe(s)*nchoosek(m+1,k1)*(-1/SNRe(s))^k1*igamma(m-n-k1+Mr-2*Mt+3,1/SNRe(s));
            end
            for k2 = 0:m
                temp2 = temp2 + ae(m+1)*nchoosek(Mr,n)*(-1)^n*SNRe(s)^(m-n+Mr-2*Mt+1)/2^(n*R)*...
                        SNRe(s)*(Mt+1/SNRe(s)-m-1)*nchoosek(m,k2)*(-1/SNRe(s))^k2*igamma(m-n-k2+Mr-2*Mt+2,1/SNRe(s));
            end
            for k3 = 0:m-1
                temp3 = temp3 + ae(m+1)*nchoosek(Mr,n)*(-1)^n*SNRe(s)^(m-n+Mr-2*Mt+1)/2^(n*R)*...
                        m*nchoosek(m-1,k3)*(-1/SNRe(s))^k3*igamma(m-n-k3+Mr-2*Mt+1,1/SNRe(s));
            end
        end
    end
    ps2(s) = exp(1/SNRe(s))*2^((Mr-Mt+1)*R)/factorial(Mr-Mt+1)*(temp1 + temp2 - temp3)/SNRr(s)^(Mr-Mt+1);
end
semilogy(SNRr_dB,po1,SNRr_dB,po2,'o',SNRr_dB,ps1,SNRr_dB,ps2,'o','linewidth',1) % plot P_out vs SNR