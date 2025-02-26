%% <<<<<<<<<< radio communication system simulation script >>>>>>>>>>>> %%

clc; clear variables; % clear all variables
close all;
addpath( 'subfunctions'); % add directory "subfunctions" to path

% global simulation parameters
ebN0dB = 0:30; % SNR (per bit) in dB

% parameters for bits generation
nMinErr=100;
nBitsPerLoop = 10e3;
nMaxBits= 10*nBitsPerLoop;

% Modulation parameters
modulationFormat = 'PSK'; % Choose between PSK or QAM
bits_per_symbol = 2; % 2 = QPSK, 4 = 16QAM ...
modulationOrder = 2^bits_per_symbol; % modulationOrder: M = 2^m ;M = 4 for QPSK; m = number of bits to be modulated;
constellation = generateConstellation(modulationFormat, modulationOrder); % constellation of the modulation format

% print if print_flag = 1
print_flag = 0;

%%% Input for channel, antenna and signal combining config %%%
%K = input("Please type in your K value: "); % K = P_LOS / P_NLOS --> for Rayleigh K = 0
K = [0, 2, 5, 10];
K_numeric = K; % this is just for the title at the end of this script
nrAntennas = input("Please type in the number of antennas: ");
transmitDiversity = ["MRC"; "EGC"; "SDC"; "sum"];
transmitDiversitySchemeInput = input("Please select your transmit diversity scheme according to their corresponding numbers: MRC=1, EGC=2, SDC=3 or SUM=4: ");
transmitDiversityScheme = transmitDiversity(transmitDiversitySchemeInput);

% VARIABLES FOR ERROR COUNTING
nErr = zeros(numel(ebN0dB),1);
ber = zeros(numel(ebN0dB),length(K));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% here goes the simulation loop...
for snr_loop = 1:numel(ebN0dB)
    for j = 1:length(K)
    snr_bit = ebN0dB(snr_loop);
    loopCnt = 0;
    nBits = nBitsPerLoop;
    nTotalErrors = 0;

        while nTotalErrors < nMinErr &&  nBits < nMaxBits 
        numberOfSymbols = nBits / bits_per_symbol;    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% transmitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = generateBits(nBits);
            txSym = mapper(data, constellation);
    
            % (potential visualization of results)
            if print_flag == 1
                scatterplot(constellation)
                if ((strcmp(modulationFormat,'PSK') == 1) && modulationOrder == 4) 
                    title(['Q',modulationFormat, ' Constellation Diagram of txSym'])
                else
                    title([int2str(modulationOrder),'-',modulationFormat, ' Constellation Diagram txSym'])
                end
                grid on;
            end
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        chSym = zeros(nrAntennas, numberOfSymbols);
        channelCoefficients = radioFadingChannel(numberOfSymbols,K(j),nrAntennas);
        chSym_radio = txSym .* channelCoefficients; 
    
    
        %%% SNR AWGN %%%%
            
            % old version of added awgn -> setSNR new function
            %chSym(:, :) = add_awgn(chSym_rayleigh, constellation, snr_bit);
            chSym(:, :) = setSNR(chSym_radio, snr_bit+10*log10(bits_per_symbol));
    
            % (potential visualization of results) -> only first antenna
            if print_flag == 1
                scatterplot(chSym(1,:));
                title("Transmitted symbols with channel coefficients and awgn");
            end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% receiver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
            % Antenna combining according to 3rd argument -> MRC, EGC, SDC, sum
            rxSym = antennaCombining(chSym, channelCoefficients, transmitDiversityScheme);
            
            % (potential visualization of results) -> only first antenna
            if print_flag == 1
                scatterplot(rxSym(1,:));
                title("Received symbols without channel coefficients");
            end
                 
            rxBits = zeros(1, nBits);
            [rxSym(:,:), rxBits(:,:)] = receiver(rxSym(:,:), constellation);
    
            % (potential visualization of results)
            if print_flag == 1
                scatterplot(rxSym(1,:))
                if ((strcmp(modulationFormat,'PSK') == 1) && modulationOrder == 4) 
                    title(['Q',modulationFormat, ' Constellation Diagram of rxSym'])
                else
                    title([int2str(modulationOrder),'-',modulationFormat, ' Constellation Diagram of rxSym'])
                end
                grid on;
            end
    
    
    
        %%%%%%%%%%%%%%%%%%%%%%%% determination of number of Errors and BER %%%%%%%%%%%%%%%%%%%%%%%%%

            [nErr(snr_loop,:), ber(snr_loop,K(j)+1)] = countErrors(rxBits(:,:), data);

            nTotalErrors = nTotalErrors + nErr(snr_loop,:);
            %calculate how many times you've run through the loop
            loopCnt = loopCnt+1;
            %calculate nBits for next loop
            nBits = nBitsPerLoop * (loopCnt+1);
        end
    end
end


%%%%%%%%%%% visualization of end results (e.g. BER vs. SNR) %%%%%%%%%%%
% calculate analytic awgn
ebN0dB_lin = 10.^(ebN0dB/10);
analytic_awgn = 0.5 * erfc(sqrt(ebN0dB_lin));

% Calculate analytic rayleigh
analytic_rayleigh = 0.5 .* (1 - sqrt(ebN0dB_lin./(1+ebN0dB_lin)));

% Calculate analytic uncle bens rice channels with K=0; 2; 5; 10
if nrAntennas == 1
    K=2;
    fun = @(theta) (((1+K)*sin(theta).^2) ./ ((1+K)*sin(theta).^2 + ebN0dB_lin)) .* exp(-1 .* ((K * ebN0dB_lin) ./ (((1+K) * sin(theta).^2) + ebN0dB_lin)));
    analytic_rice2 = (1/pi) * integral(fun, 0, pi/2, 'ArrayValued', true);
    K=5;
    fun = @(theta) (((1+K)*sin(theta).^2) ./ ((1+K)*sin(theta).^2 + ebN0dB_lin)) .* exp(-1 .* ((K * ebN0dB_lin) ./ (((1+K) * sin(theta).^2) + ebN0dB_lin)));
    analytic_rice5 = (1/pi) * integral(fun, 0, pi/2, 'ArrayValued', true);
    K=10;
    fun = @(theta) (((1+K)*sin(theta).^2) ./ ((1+K)*sin(theta).^2 + ebN0dB_lin)) .* exp(-1 .* ((K * ebN0dB_lin) ./ (((1+K) * sin(theta).^2) + ebN0dB_lin)));
    analytic_rice10 = (1/pi) * integral(fun, 0, pi/2, 'ArrayValued', true);
else %if nr Antennas > 1
    analytic_rayleigh = berfading(ebN0dB./nrAntennas, 'qam', modulationOrder, nrAntennas, 0); % K=0
    analytic_rice2 = berfading(ebN0dB, 'qam', modulationOrder, nrAntennas, 2);
    analytic_rice5 = berfading(ebN0dB, 'qam', modulationOrder, nrAntennas, 5);
    analytic_rice10 = berfading(ebN0dB, 'qam', modulationOrder, nrAntennas, 10);
end

% Calculate BER for SDC with QPSK in a Rayleigh
rayleighSDCfun = @(gamma) 0.5 * erfc(sqrt(2*ebN0dB_lin.*gamma)./sqrt(2)) * (1-exp(-gamma)^(nrAntennas-1))*exp(-gamma);
analyticRayleighSdcBer = nrAntennas * integral(rayleighSDCfun, 0, Inf, 'ArrayValued',true);

% Calculate BER for MRC with QPSK in a Rayleigh
rayleighMRCfun = @(gamma) 0.5 * erfc(sqrt(2*ebN0dB_lin.*gamma)./sqrt(2)) * gamma^(nrAntennas-1)*exp(-gamma);
analyticRayleighMrcBer = 1/gamma(nrAntennas) * integral(rayleighSDCfun, 0, Inf, 'ArrayValued',true);

% Plot numeric and analytic
figure
semilogy(ebN0dB, ber(:,1), '.','MarkerSize', 20, 'DisplayName', 'Numeric K=0', 'color', 'blue');
hold on;
semilogy(ebN0dB, ber(:,3), '.','MarkerSize', 20, 'DisplayName', 'Numeric K=2', 'color', 'red');
semilogy(ebN0dB, ber(:,6), '.','MarkerSize', 20, 'DisplayName', 'Numeric K=5', 'color', 'green');
semilogy(ebN0dB, ber(:,11), '.','MarkerSize', 20, 'DisplayName', 'Numeric K=10', 'color', 'magenta');
semilogy(ebN0dB, analytic_awgn, 'DisplayName', 'AWGN', 'LineWidth', 2);
semilogy(ebN0dB, analytic_rayleigh, 'DisplayName', 'Analytic Rayleigh', 'LineWidth', 2, 'color', 'blue');
semilogy(ebN0dB, analytic_rice2, 'DisplayName', 'Analytic Rice K = 2', 'LineWidth', 2, 'color', 'red');
semilogy(ebN0dB, analytic_rice5, 'DisplayName', 'Analytic Rice K = 5', 'LineWidth', 2, 'color', 'green');
semilogy(ebN0dB, analytic_rice10, 'DisplayName', 'Analytic Rice K = 10', 'LineWidth', 2, 'color', 'magenta');
semilogy(ebN0dB, analyticRayleighSdcBer, 'DisplayName', 'Analytic Rayleigh SDC', 'LineWidth', 2);
semilogy(ebN0dB, analyticRayleighMrcBer, 'DisplayName', 'Analytic Rayleigh MRC', 'LineWidth', 2);
legend('Location', 'best');
xlabel("SNR per Bit in dB");
ylabel("BER");
xlim([0 30])
ylim([10e-7 0.5])
titleString1 = "Comparision of analytic data with"; 
titleString2 = " " + modulationOrder + '-' + modulationFormat + " with K = " + K_numeric + " for " + nrAntennas + " antennas with " + transmitDiversityScheme;
title(titleString1, titleString2);
%title([int2str(modulationOrder),'-',modulationFormat,' with K = ',int2str(K),' for ',int2str(nrAntennas),' antennas with ',transmitDiversityScheme], 'Interpreter', 'none')
grid on;