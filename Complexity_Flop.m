clear all
Mt = [1:10]; % Tx antenna
N1 = 1; % iteration number for the power iteration method is set to 1
N2 = 10; % iteration number for the power iteration method is set to 10

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mr = Mt; % Rx antenna (case 1: Mt = Mr = Me)
Me = Mr; % eavesdropping antenna
TR1 = 2*Mt.*Mr.^2+2*Mt.*Mr+2*Mt.*Me+2*Mt+(2*N1-1)*Mr.^2+2*N1*Mr+2*Me; % TR complexity for N1
TR2 = 2*Mt.*Mr.^2+2*Mt.*Mr+2*Mt.*Me+2*Mt+(2*N2-1)*Mr.^2+2*N2*Mr+2*Me; % TR complexity for N2
ZF = 2*Mt.^2+4*Mt.*Mr+4*Mt.*Me-Mr-Me+2; % ZF complexity
MMSE = 2*Mt.^2+4*Mt.*Mr+4*Mt.*Me-Mr-Me+4; % MMSE complexity
plot(Mt,TR1,Mt,TR2,Mt,ZF,Mt,MMSE,'o','LineWidth',1); hold on % plot Complexity vs Mt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mr = 2*Mt; % Rx antenna (case 2: Mr = Me = 2Mt)
Me = Mr; % eavesdropping antenna
TR1 = 2*Mt.*Mr.^2+2*Mt.*Mr+2*Mt.*Me+2*Mt+(2*N1-1)*Mr.^2+2*N1*Mr+2*Me; % TR complexity for N1
TR2 = 2*Mt.*Mr.^2+2*Mt.*Mr+2*Mt.*Me+2*Mt+(2*N2-1)*Mr.^2+2*N2*Mr+2*Me; % TR complexity for N2
ZF = 2*Mt.^2+4*Mt.*Mr+4*Mt.*Me-Mr-Me+2; % ZF complexity
MMSE = 2*Mt.^2+4*Mt.*Mr+4*Mt.*Me-Mr-Me+4; % MMSE complexity
plot(Mt,TR1,'--',Mt,TR2,'--',Mt,ZF,'--',Mt,MMSE,'o','LineWidth',1); % plot Complexity vs Mt