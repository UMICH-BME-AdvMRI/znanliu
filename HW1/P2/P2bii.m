
%% Settings
T1 = 1000;
T2 = 100;

df = [-200:1:200];
Grad = [2*pi, 4*pi, 8*pi, 16*pi];

TR = 10;
TE = TR/2;

flip = 10;
flip = flip / 180 * pi;

Msig = zeros(length(Grad), length(df));


for i = 1:length(Grad)
    phi = linspace(0, Grad(i)-Grad(i)/length(df), length(df));

    for j = 1:length(df)
        Rflip = yrot(flip);
        [Ate, Bte] = freeprecess(TE, T1, T2, df(j));
        [Atr, Btr] = freeprecess(TR-TE, T1, T2, df(j));
        Atr = zrot(phi(j)) * Atr;

        Mss = inv(eye(3) - Ate*Rflip*Atr) * (Ate*Rflip*Btr + Bte);
        Msig(i, j) = Mss(1) + 1i*Mss(2);
    end
end

%% Figures

figure(1)
hold on; grid on;
for i = 1:length(Grad)
    signal = abs(sum(Msig(i, :)));
    scatter(Grad(i), signal, 'filled', 'sizedata', 100)
end
xticks([Grad]);

xlabel('Dephasing', 'FontSize', 16)
ylabel('Steady State Magnitude', 'FontSize', 16);
title('Gradient Spoiler', 'FontSize', 16)


