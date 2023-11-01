clc
clear all
close all

df = 0;
num = 64;

combT = [200 50; 500 100; 800 150; 1200 250;1500 300];
FA = [180; 120; 60];

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

S180 = squeeze(S(1, :, :));
figure(1)
grid on;
plot(S180', 'Linewidth', 2)
legend('T1=200, T2=50', 'T1=500, T2=100', 'T1=800, T2=150', 'T1=1200, T2=250', 'T1=1500, T2=300', 'Fontsize', 16)
axis([0 65 0 1])
xlabel('# echo', 'Fontsize', 16)
ylabel('Signal Intensity', 'Fontsize', 16)
title('\alpha=180', 'Fontsize', 16)

S120 = squeeze(S(2, :, :));
figure(2)
grid on;
plot(S120', 'Linewidth', 2)
legend('T1=200, T2=50', 'T1=500, T2=100', 'T1=800, T2=150', 'T1=1200, T2=250', 'T1=1500, T2=300', 'Fontsize', 16)
axis([0 65 0 1])
xlabel('# echo', 'Fontsize', 16)
ylabel('Signal Intensity', 'Fontsize', 16)
title('\alpha=120', 'Fontsize', 16)

S60 = squeeze(S(3, :, :));
figure(3)
grid on;
plot(S60', 'Linewidth', 2)
legend('T1=200, T2=50', 'T1=500, T2=100', 'T1=800, T2=150', 'T1=1200, T2=250', 'T1=1500, T2=300', 'Fontsize', 16)
axis([0 65 0 1])
xlabel('# echo', 'Fontsize', 16)
ylabel('Signal Intensity', 'Fontsize', 16)
title('\alpha=60', 'Fontsize', 16)


% S6 = S(:, :, 6);
% S16 = S(:, :, 16);
% S32 = S(:, :, 32);
% S48 = S(:, :, 48);
% 
% T1value = combT(:, 1);
% T2value = combT(:, 2);
% 
% figure(4)
% grid on
% scatter3(T1value, T2value, S6, 'LineWidth', 2);
% xlabel('T1', 'Fontsize', 16);
% ylabel('T2', 'FontSize', 16);
% zlabel('Signal Intensity', 'FontSize', 16);
% title('#6 echo', 'FontSize', 16)
% legend('\alpha=180', '\alpha=120', '\alpha=60', 'Fontsize', 16)
% view(25, 30);
% 
% figure(5)
% grid on
% scatter3(T1value, T2value, S16, 'LineWidth', 2);
% xlabel('T1', 'Fontsize', 16);
% ylabel('T2', 'FontSize', 16);
% zlabel('Signal Intensity', 'FontSize', 16);
% title('#16 echo', 'FontSize', 16)
% legend('\alpha=180', '\alpha=120', '\alpha=60', 'Fontsize', 16)
% view(25, 30);
% 
% figure(6)
% grid on
% scatter3(T1value, T2value, S32, 'LineWidth', 2);
% xlabel('T1', 'Fontsize', 16);
% ylabel('T2', 'FontSize', 16);
% zlabel('Signal Intensity', 'FontSize', 16);
% title('#32 echo', 'FontSize', 16)
% legend('\alpha=180', '\alpha=120', '\alpha=60', 'Fontsize', 16)
% view(25, 30);
% 
% figure(7)
% grid on
% scatter3(T1value, T2value, S48, 'LineWidth', 2);
% xlabel('T1', 'Fontsize', 16);
% ylabel('T2', 'FontSize', 16);
% zlabel('Signal Intensity', 'FontSize', 16);
% title('#48 echo', 'FontSize', 16)
% legend('\alpha=180', '\alpha=120', '\alpha=60', 'Fontsize', 16)
% view(25, 30);

