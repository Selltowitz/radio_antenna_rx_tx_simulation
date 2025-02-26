function [nErr, ber] = countErrors(rx_bits,tx_bits)
    nErr = sum(abs(tx_bits - rx_bits));
    ber = nErr / numel(rx_bits);
end