% generate constellation points, by default gray-coded
% modulationOrder: M = 2^m ;M = 4 for QPSK; m = number of bits to be modulated;
% modulationFormat: Choose between PSK or QAM
function constellationPoints = generateConstellation(modulationFormat, modulationOrder)
    % Check if the modulation order is a power of 2
    if ~isPowerOfTwo(modulationOrder)
        error('Modulation order must be a power of 2');
    end

    % Map Gray-coded binary to constellation points
    switch modulationFormat
        case 'PSK'
            constellationPoints = pskmod(0:modulationOrder-1 , modulationOrder , pi/modulationOrder, 'gray');
        case 'QAM'
            constellationPoints = qammod(0:modulationOrder-1 , modulationOrder , 'gray');
        otherwise
            error('Modulation format is not correct');
    end
end

function result = isPowerOfTwo(n)
    result = (n > 0) && bitand(n, n - 1) == 0;
end