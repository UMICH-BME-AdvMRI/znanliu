clear all
close all
clc
%% Settings
T1 = 1000; % ms
T2 = 2; % ms


gammaBar = 42.58*10^6; % Hz/T
TBW = 8;
t = [-1:0.001:1]* (10^-3); % s
position = [-1:0.01:1]*10^-2; % m
dt = t(2) - t(1);

rf_FA = [10, 30, 90];
rf_FA = rf_FA / 180 * pi;
rf_tau = 2*10^-3; % s

BW = TBW / rf_tau;
rf_sinc = sinc(BW*t);
rf_area = sum(rf_sinc)*dt;
rf_amp = rf_FA / (gammaBar * 2 * pi * rf_area)*10^2;
rf_func = zeros(length(rf_FA), length(t));
rf_ffunc = zeros(length(rf_FA), 2*length(t));

for i = 1:length(rf_FA)
    rf_func(i, :) = rf_amp(i) .* rf_sinc;
    rf_ffunc(i, :) = [rf_func(i, :) zeros(1, length(t))];
end

gradZ_amp = 18.8*10^-3; % T/m
gradZ = gradZ_amp * ones(1, length(t));
gradZ = [gradZ -0.5*gradZ];

df = 0;

Msig = zeros(length(rf_FA), length(position));
MM = zeros(3, length(position));


%% 
for i = 1:length(rf_FA)
    for j = 1:length(position)
        M = [0 0 1]';
        [A, B] = freeprecess(1000*dt/2, T1, T2, df);
        for k = 1:length(rf_ffunc)
            M = A*M + B;
            phi = 2*pi*gammaBar*position(j)*gradZ(k)*dt/2;
            M = zrot(phi)*M;
            M = throt(abs(rf_ffunc(i, k)), angle(rf_ffunc(i, k))) * M;
            M = A*M+B;

            M = zrot(phi)*M;
        end
        Msig(i, j) = M(1) + 1i*M(2);

        if i == 1
            MM(:, j) = M;
        end
    end 
end

rf_FFT_10 = fftshift(fft(rf_func(1, :)));
rf_FFT_90 = fftshift(fft(rf_func(3, :)));

% sample = linspace(1, length(rf_FFT_10), length(position));
% rf_FFT_10_Sample = interp1(rf_FFT_10, sample);
% rf_FFT_90_Sample = interp1(rf_FFT_90, sample);


%% Figure

% figure(1)
% hold on; grid on;
% plot(position*10^4, MM(1, :), 'Linewidth', 2)
% plot(position*10^4, MM(2, :), 'Linewidth', 2)
% plot(position*10^4, MM(3, :), 'Linewidth', 2)
% legend('Mx', 'My', 'Mz', 'Fontsize', 16)
% xlabel('position (mm)')
% ylabel('Magnetization')
% title('Magnetization Vector, FA=10degree')

figure(1)
hold on; grid on;
plot(position*10^4, abs(Msig(1, :)), 'Linewidth', 2)
plot(position*10^4, abs(Msig(2, :)), 'Linewidth', 2)
plot(position*10^4, abs(Msig(3, :)), 'Linewidth', 2)
legend('10 degree', '30 degree', '90 degree', 'Fontsize', 16)
xlabel('position (mm)')
ylabel('Magnetization')
title('Transverse Signal Mxy')

figure(2)
hold on; grid on;
plot(position*10^4, abs(Msig(1, :)), 'Linewidth', 2)
plot(position*10^4, abs(Msig(3, :)), 'Linewidth', 2)
plot(position*10^4, abs(rf_FFT_10(900:1100)), 'Linewidth', 2)
plot(position*10^4, abs(rf_FFT_90(900:1100)), 'Linewidth', 2)
legend('10 degree Bloch', '90 degree Bloch', '10 degree FFT', '90 degree FFT', 'Fontsize', 16)
xlabel('position (mm)')
ylabel('Magnetization')
title('Transverse Signal Mxy')


