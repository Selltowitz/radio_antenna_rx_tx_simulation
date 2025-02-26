function [y_symbol, y_bits] = receiver(x,const)

Transp_const = const.';        % transpose the vector
Dist = sqrt((real(x)-real(Transp_const)).^2 + (imag(x)-imag(Transp_const)).^2);  % Calculating the Distance between Signal and Constellationpoints
%in this case we obtain a matrix of calculated distance between every value of the signal x and all the 
%elements of constellation points vector.in every column is calculated the
%distance between an element of the signal and all the constellation points
                                                               
[~,y] = min(Dist);     % searching index of smallest Distance for every value of the recieves signal
                  %that means the position (which line) of the smallest value of every column

y_symbol=const(y);   % every element of the signal x will be replaced with the constellation points within have the smallest ditance h
% convert indices into bits
y_bits = logical(dec2bin(y-1)-'0');
% convert from character vector to 1xN binary row vector
y_bits = y_bits';
y_bits = reshape(y_bits,1,[]);
