%% <<<<<<<<<< radio communication system simulation script >>>>>>>>>>>> %%

clc; clear variables; % clear all variables
close all;
addpath( 'subfunctions'); % add directory "subfunctions" to path

% global simulation parameters
ebN0dB = 0:30; % SNR (per bit) in dB
ebN0dB_lin = 10.^(ebN0dB/10); % Linear SNR

% parameters for bits generation
nMinErr=100;
nBitsPerLoop = 10e3;
nMaxBits= 10*nBitsPerLoop;

% Modulation parameters
modulationFormat = 'PSK'; % Choose between PSK or QAM
bits_per_symbol = 2; % 2 = QPSK, 4 = 16QAM ...
modulationOrder = 2^bits_per_symbol; % modulationOrder: M = 2^m ;M = 4 for QPSK; m = number of bits to be modulated;
constellation = generateConstellation(modulationFormat, modulationOrder); % constellation of the modulation format

%%% Input for channel, antenna and signal combining config %%%
%K = input("Please type in your K value: "); % K = P_LOS / P_NLOS --> for Rayleigh K = 0
K = 0;
K_numeric = K; % this is just for the title at the end of this script
%nrAntennas = input("Please type in the number of antennas: ");
nrAntennas = [1, 2, 5, 10];
transmitDiversity = ["MRC"; "EGC"; "SDC"; "sum"];
%transmitDiversitySchemeInput = input("Please select the Combination Method: MRC=1, EGC=2, SDC=3 or SUM=4: ");
transmitDiversityScheme = transmitDiversity(1); % always 1 = MRC for this function

% VARIABLES FOR ERROR COUNTING
nErr = zeros(numel(ebN0dB),1);
BER_Analytic = zeros(numel(ebN0dB),length(nrAntennas));
BER_Simulation = zeros(numel(ebN0dB),length(nrAntennas));
analyticRayleighMrcBer = zeros(length(nrAntennas), numel(ebN0dB));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% here goes the simulation loop...
for snr_loop = 1:numel(ebN0dB)
    for k_counter = 1:length(nrAntennas)
    snr_bit = ebN0dB(snr_loop);
    loopCnt = 0;
    nBits = nBitsPerLoop;
    nTotalErrors = 0;

        while nTotalErrors < nMinErr &&  nBits < nMaxBits 
            numberOfSymbols = nBits / bits_per_symbol;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% transmitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = generateBits(nBits);
            txSym = mapper(data, constellation);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
            chSym = zeros(nrAntennas(k_counter), numberOfSymbols);
            channelCoefficients = radioFadingChannel(numberOfSymbols,K,nrAntennas(k_counter));
            chSym_radio = txSym .* channelCoefficients; 
    
        %%% SNR AWGN %%%%
            chSym(:, :) = setSNR(chSym_radio, snr_bit+10*log10(bits_per_symbol));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% receiver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
            % Antenna combining according to 3rd argument -> MRC, EGC, SDC, sum
            rxSym = antennaCombining(chSym, channelCoefficients, transmitDiversityScheme);

            rxBits = zeros(1, nBits);
            [rxSym(:,:), rxBits(:,:)] = receiver(rxSym(:,:), constellation);

    %%%%%%%%%%%%%% determination of number of Errors and BER %%%%%%%%%%%%%%
            [nErr(snr_loop,:), BER_Simulation(snr_loop,nrAntennas(k_counter)+1)] = countErrors(rxBits(:,:), data);
            nTotalErrors = nTotalErrors + nErr(snr_loop,:);

            %calculate how many times you've run through the loop
            loopCnt = loopCnt+1;
            %calculate nBits for next loop
            nBits = nBitsPerLoop * (loopCnt+1);

            if (nrAntennas(k_counter) == 1)
                % check if QPSK is used
                if ((strcmp(modulationFormat, 'PSK') == 1) && (modulationOrder == 4))
                    % Calculate analytic uncle bens rices
                    analytic_rice_fun = @(theta) (((1+K)*sin(theta).^2) ./ ((1+K)*sin(theta).^2 + ebN0dB_lin)) .* exp(-1 .* ((K * ebN0dB_lin) ./ (((1+K) * sin(theta).^2) + ebN0dB_lin)));
                    BER_Analytic(:,K+1) = transpose((1/pi) * integral(analytic_rice_fun, 0, pi/2, 'ArrayValued', true));
                    
                    % Calculate BER for SDC with QPSK in a Rayleigh
                    rayleighSDCfun = @(gamma) 0.5 * erfc(sqrt(2*ebN0dB_lin.*gamma)./sqrt(2)) * (1-exp(-gamma)^(nrAntennas(k_counter)-1))*exp(-gamma);
                    analyticRayleighSdcBer = nrAntennas(k_counter) * integral(rayleighSDCfun, 0, Inf, 'ArrayValued',true);
                    
                    % Calculate BER for MRC with QPSK in a Rayleigh
                    rayleighMRCfun = @(gamma) 0.5 * erfc(sqrt(2*ebN0dB_lin.*gamma)./sqrt(2)) * gamma^(nrAntennas(k_counter)-1)*exp(-gamma);
                    analyticRayleighMrcBer(k_counter, :) = 1/gamma(nrAntennas(k_counter)) * integral(rayleighSDCfun, 0, Inf, 'ArrayValued',true);
                else
                    print('Analytical functions for ', modulationOrder,modulationFormat, 'are not available');
                    continue;
                end
            else
                %print('Analtical functions are only for SISO system.');
                continue;
            end
        end
    end
end

%%%%%%%%%%%%% visualization of end results (e.g. BER vs. SNR) %%%%%%%%%%%%%
% calculate analytic awgn
analytic_awgn = 0.5 * erfc(sqrt(ebN0dB_lin));

% Calculate analytic rayleigh
analytic_rayleigh = 0.5 .* (1 - sqrt(ebN0dB_lin./(1+ebN0dB_lin)));


% Plot numeric and analytic
figure
semilogy(ebN0dB, BER_Simulation(:,2), '.','MarkerSize', 20, 'DisplayName', 'Nr = 1', 'color', 'blue');
hold on;
semilogy(ebN0dB, BER_Simulation(:,3), '.','MarkerSize', 20, 'DisplayName', 'Nr = 2', 'color', 'red');
semilogy(ebN0dB, BER_Simulation(:,6), '.','MarkerSize', 20, 'DisplayName', 'Nr = 5', 'color', 'green');
semilogy(ebN0dB, BER_Simulation(:,11), '.','MarkerSize', 20, 'DisplayName', 'Nr = 10', 'color', 'magenta');
semilogy(ebN0dB, analytic_awgn, 'DisplayName', 'AWGN', 'LineWidth', 2);
semilogy(ebN0dB, analytic_rayleigh, 'DisplayName', 'Analytic Rayleigh', 'LineWidth', 2, 'color', 'blue');

% semilogy(ebN0dB, BER_Analytic(:,3), 'DisplayName', 'Analytic Rice K = 2', 'LineWidth', 2, 'color', 'red');
% semilogy(ebN0dB, BER_Analytic(:,6), 'DisplayName', 'Analytic Rice K = 5', 'LineWidth', 2, 'color', 'green');
% semilogy(ebN0dB, BER_Analytic(:,11), 'DisplayName', 'Analytic Rice K = 10', 'LineWidth', 2, 'color', 'magenta');
legend('Location', 'best');
xlabel("SNR per Bit in dB");
ylabel("BER");
xlim([0 30])
ylim([10e-7 0.5])
titleString1 = "Comparision of analytic data with"; 
titleString2 = " " + modulationOrder + '-' + modulationFormat + " with K = " + K_numeric + " for " + nrAntennas(k_counter) + " antennas with " + transmitDiversityScheme;
title(titleString1, titleString2);
%title([int2str(modulationOrder),'-',modulationFormat,' with K = ',int2str(K),' for ',int2str(nrAntennas),' antennas with ',transmitDiversityScheme], 'Interpreter', 'none')
grid on;