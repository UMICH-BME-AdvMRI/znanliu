function result = RFMatrixX(alpha, phase)

    result = zeros(3, 3);
    
    result(1, 1) = cos(alpha/2)^2;
    result(1, 2) = -exp(2*-1i*phase)*sin(alpha/2)^2;
    result(1, 3) = exp(-1i*phase)*sin(alpha);

    result(2, 1) = conj(result(1, 2));
    result(2, 2) = result(1, 1);
    result(2, 3) = conj(result(1, 3));

    result(3, 1) = -0.5 * conj(result(1, 3));
    result(3, 2) = -0.5 * result(1, 3);
    result(3, 3) = cos(alpha);

end