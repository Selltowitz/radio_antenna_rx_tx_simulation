function y = radioFadingChannel(nSamp,k, Nr)
y = zeros(Nr, nSamp);
for j = 1:Nr
    % for a Rayleigh channel k=0, otherwise its a rice channel
    
    %Channel coefficients(NLOS)
    h_NLOS = randn(1,nSamp)+1i*randn(1,nSamp);  
    
    %Normalizing of LOS amplitude to P_LOS = K * P_NLOS 
    P_NLOS = mean(abs(h_NLOS).^2);   % given (linear) mean power NLOS component
    P_LOS = P_NLOS*k ;       
    Amp_h_LOS = sqrt(P_LOS);             %Amplitude von LOS Channel coefficients
    Phase = rand(size(h_NLOS)).*2*pi;    % phase distribution
    
    h_LOS = Amp_h_LOS.* exp(1j .* Phase);  %Channel coefficients (LOS)
                                        %Phase turns because of Doppler frequency; Amplitude constant 
    h_gesamt = h_LOS+h_NLOS;     % Add NLOS and LOS components
    Pgesamt = mean(abs(h_gesamt).^2);        % Normalized, complex channel coefficients
    y(j,:) = h_gesamt.*sqrt(1/(Nr*Pgesamt));
end
end