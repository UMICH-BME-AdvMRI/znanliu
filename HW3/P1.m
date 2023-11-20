%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial Fourier Imaging %
% 11/19/2023              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear all
close all

% load K-space data
data = load("Data_Assignment3_Problem1.mat");
data = data.kspaceData_SingleCoil;
N = size(data, 1);

%% Pr1a
%% reconstruct the fully-sampled K-space data
image = ifftshift(ifft2(data)); % convert from Fourier domain to image domain
figure(1)
imagesc(abs(image));
colorbar;
title("fully-sampled image", 'FontSize', 16)

% plot_name = "p1a_fully_sampled_image";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");

figure(2)
imagesc(angle(image));
colorbar;
title("fully-sampled image (phase)", 'FontSize', 16)

% plot_name = "p1a_fully_sampled_image_phase";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");

%% recostruct the zero-filled K-space data
partial_fourier_factor = 5/8;
zero_filled_data = zeros(N, N);
partial_N = N * partial_fourier_factor;
zero_filled_data(1:partial_N, :) = data(1:partial_N, :);

zero_filled_image = ifftshift(ifft2(zero_filled_data));
figure(3)
imagesc(abs(zero_filled_image));
colorbar;
title("zero-filled image", 'FontSize', 16)

% plot_name = "p1a_zero_filled_image";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");

figure(4)
imagesc(angle(zero_filled_image));
colorbar;
title("zero-filled image (phase)", 'FontSize', 16)

% plot_name = "p1a_zero_filled_image_phase";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");

%% the difference of fully-sampled and zero-filled
diff_image = abs(image) - abs(zero_filled_image);
diff_phase = angle(image) - angle(zero_filled_image);
figure(5)
imagesc(abs(diff_image));
colorbar;
title("difference", 'FontSize', 16)

% plot_name = "p1a_difference";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");

figure(6)
imagesc(diff_phase)
colorbar;
title("difference (phase)", 'FontSize', 16)

% plot_name = "p1a_difference_phase";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");



%% Pr1b
%% POCS
% set conditions to break out the loop
error = 1e-6;
Niter = 100;
pre = 0;

% 2D hanning window
hanning_2D = hanning(N)*hanning(N)';
hanning_data = hanning_2D .* data;
hanning_data(1:3/8*N, :) = 0;
hanning_data(5/8*N:end, :) = 0;
% initial phase guess
initial_phase = angle(ifftshift(ifft2(hanning_data)));

POCS_data = zeros(N, N);
POCS_data(1:partial_N, :) = data(1:partial_N, :);
POCS_image = ifftshift(ifft2(POCS_data));

% initial image
initial_image = abs(POCS_image) .* exp(1j * initial_phase);
temp_POCS_image = initial_image;

for i = 1:Niter
    temp_POCS_KSpace = fftshift(fft2(temp_POCS_image));
    % replace by the true partial fourier part
    temp_POCS_KSpace(1:partial_N, :) = data(1:partial_N, :);
    
    % compute the 2-norm of difference
    post = norm(temp_POCS_KSpace - data);
    % disp(abs(pre-post))
    if abs(pre-post) < error
        break;
    end
    pre = norm(temp_POCS_KSpace - data);
    temp_POCS_image = abs(ifftshift(ifft2(temp_POCS_KSpace)));
end
post_POCS_image = ifftshift(ifft2(temp_POCS_KSpace));

figure(7)
imagesc(abs(post_POCS_image));
colorbar;
title("POCS image")

% plot_name = "p1b_POCS_image";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");

figure(8)
imagesc(angle(post_POCS_image))
colorbar;
title("POCS image phase")

% plot_name = "p1b_POCS_image_phase";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");

%% the difference of fully-sampled and POCS
diff_image = abs(image) - abs(post_POCS_image);
diff_phase = angle(image) - angle(post_POCS_image);

figure(9)
imagesc(abs(diff_image));
colorbar;
title("difference", 'FontSize', 16)

% plot_name = "p1b_difference_POCS";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");

figure(10)
imagesc(diff_phase)
colorbar;
title("difference (phase)", 'FontSize', 16)

% plot_name = "p1b_difference_POCS_phase";
% folder_name = "plots/";
% format_name = ".jpeg";
% file_name = folder_name + plot_name + format_name;
% exportgraphics(gcf, file_name, "ContentType","image");
