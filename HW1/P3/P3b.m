clear all
close all
clc
%% Settings
T1 = 1000; % ms
T2 = 100; % ms


gammaBar = 42.58*10^6; % Hz/T
TBW = 8;
t = [-1:0.001:1]* (10^-3); % s
position = [-1:0.01:1]*10^-2; % m
dt = t(2) - t(1);

rf_FA = 90;
rf_FA = rf_FA / 180 * pi;
rf_tau = 2*10^-3; % s

BW = TBW / rf_tau;
rf_sinc = sinc(BW*t);
rf_area = sum(rf_sinc)*dt;
rf_amp = rf_FA / (gammaBar * 2 * pi * rf_area)*10^2;
rf_func = rf_amp * rf_sinc;
rf_func = [rf_func zeros(1, length(t))];

gradZ_amp = 18.8*10^-3; % T/m
gradZ = gradZ_amp * ones(1, length(t));
gradZ = [gradZ -0.5*gradZ];

df = [0 200];

Msig = zeros(1, length(position));
MM = zeros(3, length(position));


%% 
for i = 1:length(df)
    for j = 1:length(position)
        M = [0 0 1]';
        [A, B] = freeprecess(1000*dt/2, T1, T2, df(i));
        for k = 1:length(rf_func)
            M = A*M + B;
            phi = 2*pi*gammaBar*position(j)*gradZ(k)*dt/2;
            M = zrot(phi)*M;
            M = throt(abs(rf_func(k)), angle(rf_func(k))) * M;
            M = A*M+B;

            M = zrot(phi)*M;
        end
        Msig(i, j) = M(1) + 1i*M(2);

        if i == 1
            MM(:, j) = M;
        end
    end 
end

%% Figure
figure(1)
grid on;
tt = linspace(0, 4*10^-3, length(rf_func));
subplot(2, 1, 1)
plot(tt, rf_func);
xlabel('Time (s)')
ylabel('B1 Strength (T)')
title('RF Pulse')
subplot(2, 1, 2)
plot(tt, gradZ * 10^3)
xlabel('Time (s)')
ylabel('Gz (mT/m)')
title('Slice Select Gradient')

figure(2)
hold on; grid on;
plot(position*10^4, MM(1, :), 'Linewidth', 2)
plot(position*10^4, MM(2, :), 'Linewidth', 2)
plot(position*10^4, MM(3, :), 'Linewidth', 2)
legend('Mx', 'My', 'Mz', 'Fontsize', 16)
xlabel('position (mm)')
ylabel('Magnetization')
title('Magnetization Vector')

figure(3)
hold on; grid on;
plot(position*10^4, abs(Msig(1, :)), 'Linewidth', 2)
plot(position*10^4, abs(Msig(2, :)), 'Linewidth', 2)
legend('0Hz', '200Hz', 'Fontsize', 16)
xlabel('position (mm)')
ylabel('Magnetization')
title('Transverse Signal Mxy')


