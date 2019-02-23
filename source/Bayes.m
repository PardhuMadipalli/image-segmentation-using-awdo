function MED = Bayes(image)
I=im2double(image);
% imshow(I);
% title('origial image');
%BW1 = edge(I,'canny');
%% Adding noise
% I_noise = imnoise(I,'speckle',0.1);
% figure; imshow(I_noise);
% title('Image with speckle noise');
% imwrite(I_noise, 'noised.jpg', 'JPEG');
%title('speckle noise');
%% BAYES
MED= bayesEstimateDenoise(I,'sigmaSpatial', 2, 'windowSize', 2,'sigmaFactor', 2, 'sigmaMethod', 'firstrows');
MED=im2uint8(MED);
