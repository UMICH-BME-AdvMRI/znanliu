clc
clear all
close all

% load data
dataP2 = load("Data_Assignment3_Problem2.mat");
data = dataP2.kspaceData;
coil = dataP2.coilmaps;

% Fully sampled Image
[L, ~, N] = size(data);
coil_image = zeros(L, L, N);
image = zeros(L, L);

for i = 1:N
    temp_data = data(:, :, i);
    coil_image(:, :, i) = ifftshift(ifft2(temp_data));
    image = image + coil_image(:, :, i).*coil(:, :, i);
end

figure(1)
imagesc(abs(image))
colorbar;

% Aliased R=2
undersampled_data = zeros(L, L, N);
undersampled_coil_image = zeros(L, L, N);
undersampled_image = zeros(L, L);
for i = 1:(L/2)
    undersampled_data(:, i*2-1, :) = data(:, i*2-1, :);
end

for i = 1:N
    temp_data = undersampled_data(:, :, i);
    undersampled_coil_image(:, :, i) = ifftshift(ifft2(temp_data));
    undersampled_image = undersampled_image + undersampled_coil_image(:, :, i) .* coil(:, :, i);
end

figure(2)
imagesc(abs(undersampled_image))
colorbar;

% SENSE R=2 Reconstruction
SENSE_R2_image = zeros(L, L);

for i = 1:L
    for j = 1:L/2
        temp_image = squeeze(undersampled_coil_image(i, j, :));
        temp_coil = [squeeze(coil(i, j, :)), squeeze(coil(i, j+100, :))];
        temp_signal = pinv(temp_coil) * temp_image;
        SENSE_R2_image(i, j) = temp_signal(1);
        SENSE_R2_image(i, j+100) = temp_signal(2);
    end
end

figure(3)
imagesc(abs(SENSE_R2_image))
colorbar;

% SENSE R=4 Reconstruction
R = 4;
undersampled_data_R4 = zeros(L, L, N);
undersampled_coil_image_R4 = zeros(L, L, N);
undersampled_image_R4 = zeros(L, L);
for i = 1:(L/R)
    undersampled_data_R4(:, i*R-1, :) = data(:, i*R-1, :);
end

for i = 1:N
    temp_data = undersampled_data_R4(:, :, i);
    undersampled_coil_image_R4(:, :, i) = ifftshift(ifft2(temp_data));
    undersampled_image_R4 = undersampled_image_R4 + undersampled_coil_image_R4(:, :, i) .* coil(:, :, i);
end

figure(4)
imagesc(abs(undersampled_image_R4))
colorbar;

SENSE_R4_image = zeros(L, L);

for i = 1:L
    for j = 1:L/R
        temp_image_1 = squeeze(undersampled_coil_image(i, j, :));
        temp_image_2 = squeeze(undersampled_coil_image(i, j+L/R, :));
        temp_coil = [squeeze(coil(i, j+(L*0)/R, :)), squeeze(coil(i, j+(L*1)/R, :)), squeeze(coil(i, j+(L*2)/R, :)), squeeze(coil(i, j+(L*3)/R, :)),];
        temp_signal_1 = pinv(temp_coil) * temp_image_1;
        temp_signal_2 = pinv(temp_coil) * temp_image_2;
        SENSE_R4_image(i, j+(L*0)/R) = temp_signal_1(1);
        SENSE_R4_image(i, j+(L*1)/R) = temp_signal_2(2);
        SENSE_R4_image(i, j+(L*2)/R) = temp_signal_1(3);
        SENSE_R4_image(i, j+(L*3)/R) = temp_signal_2(4);
    end
end

figure(5)
imagesc(abs(SENSE_R4_image))
colorbar;














