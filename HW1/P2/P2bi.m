
%% Settings
T1 = 1000;
T2 = 100;
df = [-200:1:200];

TR = 10;
TE = TR/2;

flip = 10;
flip = flip / 180 * pi;

Msig = zeros(1, length(df));
SG = [0 0 0; 0 0 0; 0 0 1];

%%
for i = 1:length(df)
    Rflip = yrot(flip);
    [Ate, Bte] = freeprecess(TE, T1, T2, df(i));
    [Atr, Btr] = freeprecess(TR-TE, T1, T2, df(i));
    Atr = SG*Atr;
    Mss = inv(eye(3) - Ate*Rflip*Atr) * (Ate*Rflip*Btr + Bte);
    Msig(1, i) = Mss(1) + 1i*Mss(2);
end

%% Figure

figure(1)
plot(df, abs(Msig))
xlabel('Frequency (Hz)')
ylabel('Magnitude');
title('TR=5ms, TE=2.5ms')
grid on;

figure(2)
plot(df, angle(Msig))
xlabel('Frequency (Hz)')
ylabel('Phase (Radians)')
title('TR=5ms, TE=2.5ms')
axis([min(df) max(df) -pi pi])
grid on;


% clear;
% T1 = 1000;
% T2 = 100;
% alpha = 10 * pi/180;
% TE = 5;
% TR = 10;
% df = [-200:200];
% 
% %% One solution
% signal = zeros(length(df), 1);
% 
% for jj = 1:length(df)
%     M = M_ss_flash(alpha, T1, T2, TE, TR, df(jj));
%     signal(jj) = M(1) + 1j*M(2);
% end
% 
% figure
% for ii=1:3
%     plot(df, abs(signal));
%     hold on
% end
% xlabel('Frequency [Hz]')
% ylabel('Steady state signal')
% ylim([0 0.07])
% 
% function M = M_ss_flash(alpha, T1, T2, TE, TR, df)
% R = yrot(alpha);
% P = [0 0 0; 0 0 0; 0 0 1];
% [A_TE, B_TE] = freeprecess(TR-TE, T1, T2, df);
% [A_TR, B_TR] = freeprecess(TE, T1, T2, df);
% M = inv(eye(3)-A_TE*R*A_TR*P) * (A_TE*R*B_TR+B_TE);
% end
