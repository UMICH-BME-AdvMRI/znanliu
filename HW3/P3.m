%%%%%%%%%%%%%%
% GRAPPA     %
% 11/19/2023 %
%%%%%%%%%%%%%%

clc
clear all
close all

% load data
dataP3 = load("Data_Assignment3_Problem2.mat");
data = dataP3.kspaceData;
coil = dataP3.coilmaps;

%% GRAPPA K-space

ACL = 24; % autocolibration lines
kX = 3;
kY = 2;

[L, ~, N] = size(data);
coil_image = zeros(L, L, N); % image with coil sensitive
newKSpace = zeros(L, L, N); % data convoluted with fft2(coil)

for i = 1:N
    % i-th data
    temp_data = data(:, :, i);
    coil_image(:, :, i) = ifftshift(ifft2(fftshift(temp_data))) .* coil(:, :, i);
    newKSpace(:, :, i) = fftshift(fft2(ifftshift(coil_image(:, :, i))));
end

%% GRAPPA weight
start_ACL = (L - ACL)/2 + 1;
end_ACL = start_ACL + ACL - 1;

% collect 3x3 patch
patch_ACL = zeros((ACL-2)*(L-2), kX*(kY+1), N);
% select 3x2 source
s_src = zeros(size(patch_ACL, 1), kX*kY, N);
% select 1 target
s_tar = zeros(size(patch_ACL, 1), 1);

for k = 1:N
    num = 1;
    for i = start_ACL:(end_ACL-2)
        for j = 1:(L-2)
            patch_ACL(num, :, k) = reshape(newKSpace(i:i+2, j:j+2, k), [1, 9])';
            num = num + 1;
        end
    end
    s_src(:, :, k) = patch_ACL(:, [1:3, 7:9], k);
    s_tar(:, k) = patch_ACL(:, 5, k);
end

S = [s_src(:, :, 1) s_src(:, :, 2) s_src(:, :, 3) s_src(:, :, 4) s_src(:, :, 5) s_src(:, :, 6) s_src(:, :, 7) s_src(:, :, 8)];
T = [s_tar(:, 1) s_tar(:, 2) s_tar(:, 3) s_tar(:, 4) s_tar(:, 5) s_tar(:, 6) s_tar(:, 7) s_tar(:, 8)];
W = pinv(S) * T;

%% GRAPPA recon
patch_temp = zeros((L/2-1)*(L-2), kX*(kY+1), N);
t_src = zeros(size(patch_temp, 1), kX*kY, N);
t_tar = zeros(L/2-1, L-2, N);

for k = 1:N
    num = 1;
    for i = 1:2:(L-2)
        for j = 1:(L-2)
            patch_temp(num, :, k) = reshape(newKSpace(i:i+2, j:j+2, k), [1, 9])';
            num = num + 1;
        end
    end
    t_src(:, :, k) = patch_temp(:, [1:3, 7:9], k);
end

new_S = [t_src(:, :, 1) t_src(:, :, 2) t_src(:, :, 3) t_src(:, :, 4) t_src(:, :, 5) t_src(:, :, 6) t_src(:, :, 7) t_src(:, :, 8)];
new_T = new_S * W;


%% GRAPPA Image
GRAPPA_matrix = zeros(L, L, N);
image = zeros(L, L);

for i = 1:N
    % reshape to 99*198
    t_tar(:, :, i) = reshape(new_T(:, i), [L-2, L/2-1])';
    GRAPPA_matrix(:, :, i) = newKSpace(:, :, i);
    GRAPPA_matrix(2:2:end-1, 2:end-1, i) = t_tar(:, :, i);
    image = image + ifftshift(ifft2(fftshift(GRAPPA_matrix(:, :, i))));
end

figure(1);
imagesc(abs(image));
title("GRAPPA Recon")
colorbar;

% plot_name = "p3_GRAPPA";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");




