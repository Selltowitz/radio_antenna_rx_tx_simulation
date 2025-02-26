function y = setSNR(x, snrdB)

% mean power of the signal x
nrAntennas = length(x(:,1));
P_Signal = zeros(nrAntennas, 1);
    for j = 1:nrAntennas
        P_Signal(j,:) = mean(abs(x(j,:)) .^2); 
        %P_Signal = mean(abs(x) .^2);

        % Linear Noise power PN ( SNRdB = 10*log10(Ps/PN) )
        P_Noise = P_Signal * 10.^(-snrdB / 10);

        x_size=length(x); %number of elements of the signal x

        % imaginary and real gauss distribution of the Noise
        Noise = randn(1, x_size) + 1j.*randn(1,x_size);

        % mean power of the Gauss distributed noise
        P_Rauschen = mean (abs(Noise) .^2);

        % scaling Noise with ?(Noise Power/mean power of gauss distibuted noise)
        Alpha = sqrt(P_Noise / P_Rauschen);     % scaling factor

        Noise = Alpha.*Noise;

        y = x + Noise;
    end
end