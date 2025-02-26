 %%% Mapper bits to complex values for QPSK, 8PSK and 16QAM %%%%%
 function y = transmitter(bits, constellation)

% Überprüfe wie viele Bits zusammengefasst werden pro Symbol
if numel(constellation) == 4 %QPSK
    divider_bits = 2;
elseif numel(constellation) == 8 % 8PSK
    divider_bits = 3;
else
    divider_bits = 4; % 16-QAM
end

remaining_bits = mod(numel(bits), divider_bits);
if remaining_bits == 0
    bits_to_add = 0;
else 
    bits_to_add = divider_bits - remaining_bits;
end
% fehlende Bits hinzufügen
bits = [bits, zeros(1, bits_to_add)];
%bits = transpose(bits);
% Berechne Zeilen
zeilen = numel(bits)/divider_bits;

% Reshape bits modulation specific wise
bits_rearranged = reshape(bits, [], zeilen)';

% convert bits into decimals -> LSB/MSB fragen!!!!!
bits_decimal = bi2de(bits_rearranged, 'left-msb');

% init result vector y
y = zeros(zeilen, 1);

% map from decimal values to complex values
for counter = 1:zeilen
   y(counter) = constellation(bits_decimal(counter)+1); 
end

% uncomment following line to print the symbols
%scatterplot(y) ;

% Switch case for title
switch numel(constellation)
    case 4
        title("QPSK");
    case 8
        title("8PSK");
    case 16
        title("16-QAM");
end
y = transpose(y);
grid on;
end