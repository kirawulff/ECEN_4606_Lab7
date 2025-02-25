
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
Y = abs(log(ifftshift(X)).^2);
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
Y = abs(log(ifftshift(X)).^2);
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
        