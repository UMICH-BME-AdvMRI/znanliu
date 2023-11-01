function result = gradMatrix(S, gradNum)
% S: previous static matrix
% gradNum: number of gradient applied

    [m, n] = size(S);
    if abs(S(1, end)) < 10^-16
        n = n + gradNum - 1;
    else
        n = n + gradNum;
        S = [S, [0;0;0]];
    end

    newS = zeros(m, n);
    newS(1, 2:end) = S(1, 1:end-1);
    newS(2, 1:end-1) = S(2, 2:end);
    newS(3, :) = S(3, :);
    newS(1, 1) = conj(newS(2, 1));

    result = newS;

end