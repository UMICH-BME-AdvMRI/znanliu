

%% Settings
T1 = 1000;
T2 = 100;

df = 0;

TR = 10;
TE = TR / 2;

flip = [0:0.01:1] * pi;
inc = 117/180 * pi;
Nex = 100;

Msig = zeros(1, length(flip));

%%
for k = 1:length(flip)
    Nf = 100;
    rf_phi = [1:Nf] / Nf *2 * pi;
    M = zeros(3, Nf, Nex+1);

    [Ate, Bte] = freeprecess(TE, T1, T2, df);
    [Atr, Btr] = freeprecess(TR-TE, T1, T2, df);

    M = [zeros(2, Nf); ones(1, Nf)];
    on = ones(1, Nf);

    Rfph = 0;
    Rfinc = inc;

    for n = 1:Nex
        A = Ate*throt(flip(k), Rfph);
        B = Bte;
        M = A*M+B*on;

        Msig(k) = mean(M(1, :) + 1i*M(2, :)) *exp(-1i*Rfph);
        M = Atr*M+Btr*on;

        for i = 1:Nf
            M(:, i) = zrot(rf_phi(i))*M(:, i);
        end

        Rfph = Rfph + Rfinc;
%         Rfinc = Rfinc + Rfinc;
    end
end

%% Figures
plot(flip, abs(Msig))
xlabel('FA (Rad)')
ylabel('Magnitude')



