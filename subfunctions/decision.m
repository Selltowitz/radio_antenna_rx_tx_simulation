function y = decision(x, constellation)
    % Normalize the received signal based on the average power of the constellation
    normilizedSignal = x / mean(abs(constellation).^2);
    % Initialize the output vector
    y = zeros(1, numel(x));
    for counter = 1:numel(x)
        % Calculate the distance between the received symbol and each constellation point
        distances = abs(normilizedSignal(counter) - constellation);

        % Find the index of the closest constellation point
        [~, index] = min(distances);

        % Assign the decided symbol to the output vector
        y(counter) = constellation(index);
    end
    %y = y.';
end