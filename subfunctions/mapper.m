function y = mapper(bits, constellation)
    % Check if the modulation order is a power of 2
    if ~isPowerOfTwo(numel(constellation))
        error('Modulation order must be a power of 2');
    end

    % wie viele Bits zusammengefasst werden pro Symbol
    divider_bits = log2(numel(constellation));

    remaining_bits = mod(numel(bits), divider_bits);

    if remaining_bits == 0
        bits_to_add = 0;
    else
        bits_to_add = divider_bits - remaining_bits;
    end

    % if mod(length(bits), log2(numel(constellation))) ~= 0
    %     bits = [bits, 0]; % Append a zero to make it even
    % end

    % fehlende Bits hinzufügen
    bits = [bits, zeros(1, bits_to_add)];
    
    % Berechne Zeilen
    zeilen = numel(bits)/divider_bits;
    
    % Reshape bits modulation specific wise
    bits_rearranged = reshape(bits, log2(numel(constellation)), zeilen).';
    
    % convert bits into decimals -> LSB/MSB fragen!!!!!
    bits_decimal = bi2de(bits_rearranged, 'left-msb');
    
    % init result vector y
    y = zeros(zeilen, 1);
    
    % map from decimal values to complex values
    for counter = 1:zeilen
        y(counter) = constellation(bits_decimal(counter)+1); 
    end
    y = transpose(y);
end

function result = isPowerOfTwo(n)
    result = (n > 0) && bitand(n, n - 1) == 0;
end