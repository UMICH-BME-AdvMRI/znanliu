
T_1 = 200:100:1500;
T_2 = 50:30:300;

df = 0;
num = 64;
[t1, t2] = meshgrid(T_1, T_2);
combT = [t1(:), t2(:)];
FA = [180; 120; 60];
echo = [6, 16, 32, 64];

refocT = 5;
delayT = refocT/2;

S = zeros(length(FA), length(combT), 64);

M0 = [0; 0; 1];

excitation = 90;
excitation = excitation/180 * pi;
M1 = RFMatrixX(excitation, 0)*M0;

for k = 1:length(FA)
    alpha = FA(k);
    firstFA = alpha+(90 - alpha/2);
    alpha = alpha/180*pi;
    firstFA = firstFA/180*pi;

    for j = 1:length(combT)
        T1 = combT(j, 1);
        T2 = combT(j, 2);
        [A, B] = freeprecess(delayT, T1, T2, df);
        M2 = A*M1+B;
        M3 = gradMatrix(M2, 1);
        M4 = RFMatrixY(firstFA, 0)*M3;
        for i = 1:num-1
            M5 = gradMatrix(M4, 1);
            M6 = A*M5+B;
            S(k, j, i) = abs(M6(1, 1));
            M7 = A*M6+B;
            M8 = gradMatrix(M7, 1);
            M4 = RFMatrixY(alpha, 0)*M8;
        end
        M5 = A*M4+B;
        M6 = gradMatrix(M5, 1);
        S(k, j, end) = abs(M6(1, 1));
    end
end

for k = 1:length(FA)
    for j = 1:length(echo)
        tempS = S(k, :, echo(j));
        tempV = reshape(tempS, length(T_2), length(T_1));
        figure;
        pcolor(t1, t2, tempV);
        colorbar
        caxis([0, 0.9])
        shading flat
        xlabel('T1 value')
        ylabel('T2 value')
        title(sprintf('\\alpha=%d, echo=%d', FA(k), echo(j)))
    end
end


