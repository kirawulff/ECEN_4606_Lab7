
%% Problem 1
%Problem 1
% a.
test = -10:9;
test2 = fftshift(test)
test3 = ifftshift(test2)

% b
figure
X = fft2(fftshift(SquareWave));

% c
Y = abs(log(ifftshift(X).^2));
aticklen = 200e-3*633e-9/(4e-6);
% atickstart = (length(X)-1)*aticklen/2;
% atick = linspace(-atickstart,atickstart, 5);
Z = imresize(Y, aticklen);
imagesc(Z)

% d
r_lens = .254/2;
r_lens_px = r_lens/aticklen;
mask = zeros(size(Z));
center = length(mask)/2;

% circle time
for index1 = 0:length(mask)
  for index2 = 0:length(mask)
      dist = (center - index1)^2 + (center - index2)^2;
      %fprintf("dist %d", dist);
      if (dist < r_lens_px^2)
          %fprintf("hi")
          mask(index2, index1) = 1;
      end
  end
end
imagesc(mask);
final = mask.*Z;
imagesc(final);

%% Part 2

% Part 2
k = 1;
theta = .005;
SquareWave2 = zeros(size(SquareWave));
for index1 = 1:length(SquareWave)
    % assuming a vertical tilt
    for index2 = 1:length(SquareWave)
        SquareWave2(index2, index1) = SquareWave(index2, index1) * cos(1i*k*theta*index1);
    end
end

% b
figure
X = fft2(fftshift(SquareWave2));

% c
Y = abs(log(ifftshift(X).^2));
aticklen = 200e-3*633e-9/(4e-6);
% atickstart = (length(X)-1)*aticklen/2;
% atick = linspace(-atickstart,atickstart, 5);
Z = imresize(Y, aticklen);
imagesc(Z)

% d
r_lens = .254/2;
r_lens_px = r_lens/aticklen;
mask = zeros(length(Z));
center = length(mask)/2;

% circle time
for index1 = 0:length(mask)
  for index2 = 0:length(mask)
      dist = (center - index1)^2 + (center - index2)^2;
      %fprintf("dist %d", dist);
      if (dist < r_lens_px^2)
          %fprintf("hi")
          mask(index2, index1) = 1;
      end
  end
end
d = 200e-3 * tan(.005);
d_px = d/aticklen
imagesc(mask);
final = mask.*Z;
imagesc(final);
        

%% Part 3
Hologram = imread('CT_Hologram.jpg');
Reference = imread('CT_ReferenceBeam.jpg');
Object = imread('CT_ObjectBeam.jpg');

% Without subtraction
Z = rgb2gray(Hologram);
figure
X = fft2(mat2gray(fftshift(Z)));
Y_Hologram = abs(log10(ifftshift(X).^2));
imagesc(Y_Hologram);

% With subtraction
Z = rgb2gray(Hologram)-(rgb2gray(Reference)+rgb2gray(Object));
figure
X = fft2(mat2gray(fftshift(Z)));
Y_Hologram = abs(log10(ifftshift(X).^2));
imagesc(Y_Hologram);

%% Problem 4

Hologram = imread('USAF_Hologram.jpg');
Reference = imread('USAF_Reference.jpg');
Object = imread('USAF_ObjectBeam.jpg');

% With subtraction
Z = rgb2gray(Hologram)-(rgb2gray(Reference)+rgb2gray(Object));
figure
X = fft2(mat2gray(fftshift(Z)));
Y_Hologram = abs(log10(ifftshift(X).^2));
imagesc(Y_Hologram);

x = linspace(-1,1,length(Hologram));
y = linspace(-1,1,length(Hologram));
[X_grid, Y_grid] = meshgrid(x,y);

% To find best looking phase
% for c = -2:.2:2
% grid_phase = exp(c * (X_grid.^2 + Y_grid.^2));
% 
% G_Hologram = double(rgb2gray(Hologram)) .* grid_phase;
% Z = uint8(G_Hologram) - (rgb2gray(Reference)+rgb2gray(Object));
% figure
% 
% X = fft2(mat2gray(fftshift(Z)));
% Y_Hologram = abs(log10(ifftshift(X).^2));
% imagesc(Y_Hologram);
% title(['c = ', num2str(c)]);
% end

c = 1.4;
grid_phase = exp(c * (X_grid.^2 + Y_grid.^2));

G_Hologram = double(rgb2gray(Hologram)) .* grid_phase;
Z = uint8(G_Hologram) - (rgb2gray(Reference)+rgb2gray(Object));
figure

X = fft2(mat2gray(fftshift(Z)));
Y_Hologram = abs(log10(ifftshift(X).^2));
imagesc(Y_Hologram);
title(['c = ', num2str(c)]);
