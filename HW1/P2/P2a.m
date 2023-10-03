

%% Settings
T1 = 1000;
T2 = 100;
df = [-200:1:200];

TR = [5, 10, 20];
TE = TR/2;

flip = 60;
flip = flip / 180 * pi;

Msig = zeros(length(TR), length(df));

%% Main
for i = 1:length(TR)
    for j = 1:length(df)
        Rflip = yrot(flip);
        [Atr, Btr] = freeprecess(TR(i)-TE(i), T1, T2, df(j));
        [Ate, Bte] = freeprecess(TE(i), T1, T2, df(j));
    
        Mss = inv(eye(3) - Ate*Rflip*Atr) * (Ate*Rflip*Btr + Bte);
    
        Msig(i, j) = Mss(1) + 1i*Mss(2);
    end
end


%% Figure
figure(1)
plot(df, abs(Msig), 'Linewidth', 2)
xlabel('Frequency (Hz)', 'Fontsize', 12)
ylabel('Magnitude', 'Fontsize', 12);
title('bSSFP Frequency Responce', 'Fontsize', 16)
legend('TR=5ms', 'TR=10ms', 'TR=20ms', 'Fontsize', 16)
grid on;

figure(2)
plot(df, angle(Msig), 'Linewidth', 2)
xlabel('Frequency (Hz)', 'Fontsize', 12)
ylabel('Phase (Radians)', 'Fontsize', 12)
title('bSSFP Phase Responce', 'Fontsize', 16)
axis([min(df) max(df) -pi pi])
legend('TR=5ms', 'TR=10ms', 'TR=20ms', 'Fontsize', 16)
grid on;



