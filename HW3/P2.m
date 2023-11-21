%%%%%%%%%%%%%%
% SENSE      %
% 11/19/2023 %
%%%%%%%%%%%%%%

clc
clear all
close all

% load data
dataP2 = load("Data_Assignment3_Problem2.mat");
data = dataP2.kspaceData; % K-space data
coil = dataP2.coilmaps; % coil map data

%% P2a
%% Fully sampled Image
[L, ~, N] = size(data);
coil_image = zeros(L, L, N);
image = zeros(L, L);

for i = 1:N
    temp_data = data(:, :, i); % i-th K-space data
    coil_image(:, :, i) = conj(ifftshift(ifft2(temp_data))); % i-th image
    image = image + coil_image(:, :, i).*coil(:, :, i); % coil-combined image
end

figure(1)
imagesc(abs(image))
colorbar;
title("fully-sampled image")

% plot_name = "p2a_fully_sampled_image";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");

%% P2b
%% Aliased R=2
undersampled_data = zeros(L, L, N);
undersampled_coil_image = zeros(L, L, N);
undersampled_image = zeros(L, L);

% phase encoding direction is oriented vertically
for i = 1:(L/2)
    % undersampled K-space data
%     undersampled_data(:, i*2-1, :) = data(:, i*2-1, :);
    undersampled_data(i*2-1, :, :) = data(i*2-1, :, :);
end

for i = 1:N
    % i-th undersampled K-space data
    temp_data = undersampled_data(:, :, i);
    % i-th undersampled image data
    undersampled_coil_image(:, :, i) = conj(ifftshift(ifft2(temp_data)));
    % coil-combined undersampled image data
    undersampled_image = undersampled_image + undersampled_coil_image(:, :, i) .* coil(:, :, i);
end

figure(2)
imagesc(abs(undersampled_image))
colorbar;
title("undersampled image R=2")

% plot_name = "p2b_undersampled_image";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");

%% P2c
%% SENSE R=2 Reconstruction
R=2;
SENSE_R2_image = zeros(L, L);

% i = frequency encoding
% for i = 1:L
%     % j = phase encoding
%     for j = 1:L/R
%         % signal matrix at point (i, j), R^{8x1}
%         temp_image = conj(squeeze(undersampled_coil_image(i, j, :)));
%         % coil-sensitive matrix at point (i, j) and (i, j+L/2), R^{8x2}
%         temp_coil = [squeeze(coil(i, j, :)), squeeze(coil(i, j+100, :))];
%         % pseudo-inverse calculate the true signal at point (i, j) and (i,
%         % j+L/R)
%         temp_signal = pinv(temp_coil) * temp_image;
%         SENSE_R2_image(i, j) = temp_signal(1);
%         SENSE_R2_image(i, j+100) = temp_signal(2);
%     end
% end

for j = 1:L
    % j = phase encoding
    for i = 1:L/R
        % signal matrix at point (i, j), R^{8x1}
        temp_image = conj(squeeze(undersampled_coil_image(i, j, :)));
        % coil-sensitive matrix at point (i, j) and (i, j+L/2), R^{8x2}
        temp_coil = [squeeze(coil(i, j, :)), squeeze(coil(i+100, j, :))];
        % pseudo-inverse calculate the true signal at point (i, j) and (i,
        % j+L/R)
        temp_signal = pinv(temp_coil) * temp_image;
        SENSE_R2_image(i, j) = temp_signal(1);
        SENSE_R2_image(i+100, j) = temp_signal(2);
    end
end

figure(3)
imagesc(abs(SENSE_R2_image))
colorbar;
title("SENSE-R2-image")

% plot_name = "p2c_sense_reconstruct";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");

%% P2d
%% SENSE R=4 Reconstruction
R = 4;
undersampled_data_R4 = zeros(L, L, N);
undersampled_coil_image_R4 = zeros(L, L, N);
undersampled_image_R4 = zeros(L, L);
for i = 1:(L/R)
%     undersampled_data_R4(:, i*R-1, :) = data(:, i*R-1, :);
    undersampled_data_R4(i*R-1, :, :) = data(i*R-1, :, :);
end

for i = 1:N
    temp_data = undersampled_data_R4(:, :, i);
    undersampled_coil_image_R4(:, :, i) = conj(ifftshift(ifft2(temp_data)));
    undersampled_image_R4 = undersampled_image_R4 + undersampled_coil_image_R4(:, :, i) .* coil(:, :, i);
end

figure(4)
imagesc(abs(undersampled_image_R4))
colorbar;
title("undersampled image R=4")

% plot_name = "p2d_undersampled_image_R4";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");

figure(5)
imagesc(abs(image) - 2*abs(SENSE_R2_image))
colorbar;
title("diff SENSE-R2-image")

SENSE_R4_image = zeros(L, L);

% for i = 1:L
%     for j = 1:L/R
%         temp_image_1 = squeeze(undersampled_coil_image(i, j, :));
%         temp_image_2 = squeeze(undersampled_coil_image(i, j+L/R, :));
%         temp_coil = [squeeze(coil(i, j+(L*0)/R, :)), squeeze(coil(i, j+(L*1)/R, :)), squeeze(coil(i, j+(L*2)/R, :)), squeeze(coil(i, j+(L*3)/R, :)),];
%         temp_signal_1 = pinv(temp_coil) * temp_image_1;
%         temp_signal_2 = pinv(temp_coil) * temp_image_2;
%         SENSE_R4_image(i, j+(L*0)/R) = temp_signal_1(1);
%         SENSE_R4_image(i, j+(L*1)/R) = temp_signal_2(2);
%         SENSE_R4_image(i, j+(L*2)/R) = temp_signal_1(3);
%         SENSE_R4_image(i, j+(L*3)/R) = temp_signal_2(4);
%     end
% end

for j = 1:L
    for i = 1:L/R
        temp_image = conj(squeeze(undersampled_coil_image_R4(i, j, :)));
        temp_coil = [squeeze(coil(i+(L*0)/R, j, :)), squeeze(coil(i+(L*1)/R, j, :)), squeeze(coil(i+(L*2)/R, j, :)), squeeze(coil(i+(L*3)/R, j, :))];
        temp_signal_1 = pinv(temp_coil) * temp_image;
        SENSE_R4_image(i+(L*0)/R, j) = temp_signal_1(1);
        SENSE_R4_image(i+(L*1)/R, j) = temp_signal_1(2);
        SENSE_R4_image(i+(L*2)/R, j) = temp_signal_1(3);
        SENSE_R4_image(i+(L*3)/R, j) = temp_signal_1(4);
    end
end

figure(6)
imagesc(abs(SENSE_R4_image))
colorbar;
title("SENSE-R4-image")

% plot_name = "p2d_sense_reconstruct_R4";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");

figure(7)
imagesc(abs(image) - 4*abs(SENSE_R4_image))
colorbar;
title("diff SENSE-R4-image")











