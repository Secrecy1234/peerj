clear all
Mt = 2; % Tx antenna
Mr = 2; % Rx antenna
K = max(Mt,Mr); % relevant parameter defined before (5)
L = min(Mt,Mr); % relevant parameter defined before (5)
Me = 2; % eavesdropping antenna
SNRr_dB = [0:3:30]; % average SNR at Rx (in dB)
SNRr = 10.^(SNRr_dB/10);
SNRe_dB = 10*ones(size(SNRr_dB)); % average SNR at eavesdropper (in dB)
SNRe = 10.^(SNRe_dB/10);
R = 1; % target secrecy rate of 1 bits/s/Hz

% Genarating the coefficients 'a' defined in (6), according to the
% following reference:
% Dighe, P. A., Mallik, R. K., and Jamuar, S. S. (2003).
% Analysis of transmit-receive diversity in Rayleigh fading.
% IEEE Transactions on Communications, 51(4):694–703.
if (K == 2) & (L == 1)
    a = [0 1];
elseif (K == 2) & (L == 2)
    a = [2 -2  2; -1 0 0];
elseif (K == 3) & (L == 1)
    a = [0 0 1];
elseif (K == 3) & (L == 2)
    a = [0 3 -4  3; 0 -3/4 -1/4 0];
elseif (K == 3) & (L == 3)
    a = [3 -6 12 -12 6; -3 3/2 -3/4 -3/8 -3/8; 1 0 0 0 0];
elseif (K == 4) & (L == 2)
    a = [0 0 4 -6  4; 0 0 -1/2 -3/8 -1/8];
elseif (K == 4) & (L == 3)
    a = [0 6 -16 27 -24 10 0; 0 -3 1 3/8 -3/4 -5/32 -15/32; 0 2/3 8/27 1/27 0 0 0]; 
elseif (K == 4) & (L == 4)
    a = [4 -12 36 -68 84 -60 20 0 0; -6 6 -6 1 -1 5/2 -5/2 35/32 -35/32; 4 -4/3 4/9 28/81 92/243 100/729 20/729 0 0; -1 0 0 0 0 0 0 0 0];
elseif (K == 6) & (L == 2)
    a = [0 0 0 0 6 -10 6 0 0; 0 0 0 0 -3/16 -5/16 -9/32 -21/128 -7/128];
elseif (K == 6) & (L == 3)
    a = [0 0 0 15 -48 75 -60 21 0 0 0;
        0 0 0 -15/8 -3/4 15/32 15/32 -21/128 -21/32 63/512 -315/512;
        0 0 0 5/27 8/27 65/243 40/243 49/729 112/6561 14/6561 0];
elseif (K == 6) & (L == 6)
    a = [ 6  -30   150   -540     1440    -2832       4080      -4200         2940        -1260           252               0              0                0               0                 0                 0                  0                0; ...
        -15   30   -60   135/2   -225/4    303/4    -1527/8    12621/32    -10941/16     33579/32     -44289/32       776853/512    -1399167/1024     999999/1024   -1072071/2048      1576575/8192      -315315/8192              0                0; ...
         20  -20    20   -40/9     40/9   -808/81     808/81   -1064/243     1064/243    -6328/6561     6328/6561     -45584/19683     45584/19683   -112112/59049    112112/59049     -560560/531441     560560/531441    -9529520/43046721 9529520/43046721; ...
        -15  15/2 -15/4  -45/32  -225/128   39/128    159/1024  -231/8192   -3129/16384   -441/2048   -15813/65536   -387387/2097152 -933471/8388608  -27027/524288 -1198197/67108864 -2207205/536870912 -315315/536870912         0                0; ...
          6  -6/5   6/25  36/125  216/625 3288/15625 6792/78125 9408/390625  8652/1953125 4788/9765625   252/9765625        0              0                0               0                 0                 0                  0                0; ...
         -1 zeros(1,18)];
elseif (K == 8) & (L == 2)
    a = [0 0 0 0 0 0 8 -14 8 0 0 0 0; 0 0 0 0 0 0 -1/16 -21/128 -15/64 -15/64 -45/256 -99/1024 -33/1024]; 
elseif (K == 8) & (L == 4)    
    a = [0 0 0 0 56 -280 696 -1036 952 -504 120 0 0 0 0 0 0 0 0; ...
        0 0 0 0 -21/4 0 27/8 63/64 -147/64 -189/64 -45/32 2871/256 -13959/512 21021/512 -39039/1024 45045/2048 -15015/2048 0 0; ...
        0 0 0 0 56/81 280/243 232/243 868/2187 -392/6561 -392/2187 40/19683 8272/59049 10912/59049 334048/1594323 992992/4782969  760760/4782969 3875872/43046721 544544/14348907 544544/43046721; ...
        0 0 0 0 -7/128 -35/256 -201/1024 -3395/16384 -5705/32768 -1953/16384 -8745/131072 -31779/1048576 -46035/4194304 -6435/2097152 -21021/33554432 -45045/536870912 -3003/536870912 0 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
po1 = zeros(size(SNRr)); % numerical integration for P_out,TR
po2 = zeros(size(SNRr)); % exact closed-form expression (17) for P_out,TR
for s = 1:length(SNRr)
    temp = 0;
    for m = 1:L
        for n = abs(Mt-Mr):(Mt+Mr-2*m)*m

            po1(s) = po1(s) + integral(@(x) x.^(Me-1).*exp(-x/SNRe(s))/(SNRe(s)^Me*factorial(Me-1)).*...
                     a(m,n+1).*gammainc(m*(2^R*x+2^R-1)/SNRr(s),n+1),0,inf);
                    
            for k = 0:n
                for l = 0:k          
                    temp = temp + a(m,n+1)*exp(-m*(2^R-1)/SNRr(s))*(m/SNRr(s))^k*...
                           factorial(l+Me-1)*2^(l*R)*(2^R-1)^(k-l)/(factorial(l)*factorial(k-l))*...
                           (m*2^R/SNRr(s)+1/SNRe(s))^(-l-Me);
                end
            end
        end
    end
    po2(s) = 1 - 1/(SNRe(s)^Me*factorial(Me-1))*temp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ps1 = zeros(size(SNRr)); % asymptotic numerical integration for P_out,TR
ps2 = zeros(size(SNRr)); % asymptotic closed-form expression (24) for P_out,TR
for s = 1:length(SNRr)
    temp = 0;
    ps1(s) = ps1(s) + integral(@(x) x.^(Me-1).*exp(-x/SNRe(s))/(SNRe(s)^Me*factorial(Me-1)).*...
             prod(factorial(L-[1:L]))/prod(factorial(Mt+Mr-[1:L])).*((2^R*x+2^R-1)/SNRr(s)).^(Mt*Mr),0,inf);

    for n = 0:Mt*Mr         
        temp = temp + nchoosek(Mt*Mr,n)*factorial(n+Me-1)*2^(n*R)*(2^R-1)^(Mt*Mr-n)*...
                    SNRe(s)^n;
    end
    ps2(s) = prod(factorial(L-[1:L]))/(prod(factorial(Mt+Mr-[1:L]))*factorial(Me-1))*SNRr(s)^(-Mt*Mr)*temp;
end
semilogy(SNRr_dB,po1,SNRr_dB,po2,'o',SNRr_dB,ps1,SNRr_dB,ps2,'o','linewidth',1) % plot P_out vs SNR