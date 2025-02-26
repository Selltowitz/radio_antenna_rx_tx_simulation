function y = add_awgn(tx_symbols, constellation, snr_bit)
    rows = numel(tx_symbols(:,1));
    M = numel(constellation);
    for counter = 1:rows
        snr_linear = 10.^(snr_bit/10); 
        snr_symbol_linear = snr_linear * log2(M);
        snr_symbol_db = 10 * log10(snr_symbol_linear); % convert back to dB
    
        y(counter,:) = awgn(tx_symbols(counter,:), snr_symbol_db, 'measured');
    end
end