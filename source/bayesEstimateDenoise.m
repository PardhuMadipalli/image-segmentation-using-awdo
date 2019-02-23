function res = bayesEstimateDenoise(img, varargin)
% RES = BAYESESTIMATEDENOISE(IMG, VARARGIN);
%
% This function implements the denoising algorithm described in 
% "General Bayesian estimation for speckle noise reduction in optical
% coherence tomography retinal imagery" by Alexander Wong, Akshaya Mishra,
% Kostadinka Bizheva and David A.Clausi, Optics Express, April 2010, 
% Vol. 18, No. 8, 8338--8352
%
% Parameters:
% img: Input image (preferably log-scaled OCT image)
% varargin: Additional parameters. If not set, default
% variables are used.
%   'sigmaSpatial': Standard deviation of the gaussian spatial 
%       sampling distribution
%   'windowSize': (Width-1)/2 of the neighbourhood comparison window
%   'numSamples': Number of samples for the probability estimation
%   'sigmaFactor': Factor + Sigma = treshold for sample inclusion 
%   'sigmaMethod': 'local' or 'firstrows'
% res: Denoised image
%
% Calling example:
%
% res = bayesEstimateDenoise(img, 'sigmaSpatial', 2, 'windowSize', 2, ...
%   'sigmaFactor', 2)
%
% This version of the code was revised and commented in December 2011
%
% Authors: Markus Mayer, Martin Kraus, Pattern Recognition Lab, University
% of Erlangen-Nuremberg, markus.mayer@informatik.uni-erlangen.de
%
% You may use this code as you want. I would be grateful if you would go to
% my homepage look for articles that you find worth citing in your next
% publication:
% http://www5.informatik.uni-erlangen.de/en/our-team/mayer-markus
% Thanks, Markus

% Standard parameter set
params.sigmaSpatial = 1.5;          %1.5%
params.windowSize = 3;              %3%
params.numSamples = 300;            %70%
params.sigmaFactor = 2;             %2%
params.sigmaMethod = 'local';

% Get the optional parameters from the varargin cell
if (~isempty(varargin) && iscell(varargin{1}))
    varargin = varargin{1};
end
for k = 1:2:length(varargin)
    if (strcmp(varargin{k}, 'sigmaSpatial'))
        params.sigmaSpatial = varargin{k+1};
    elseif (strcmp(varargin{k}, 'windowSize'))
        params.windowSize = varargin{k+1};
    elseif (strcmp(varargin{k}, 'numSamples'))
        params.numSamples = varargin{k+1};
    elseif (strcmp(varargin{k}, 'sigmaFactor'))
        params.sigmaFactor = varargin{k+1};
    elseif (strcmp(varargin{k}, 'sigmaMethod'))
        params.sigmaMethod = varargin{k+1};
    end
end

% for i = 1:size(img,1)
%     for j = 1:size(img,2)
%         if img(i,j)==0
%             img(i,j)=img(i,j);
%         end
%         img(i,j)=-255*log(img(i,j));
%         
%     end
% end

img=log(img+1);
max(max(img));
new=zeros(size(img));


% 
% for i = 1:size(img,1)
%     for j = 1:size(img,2)
%         img(i,j)=exp(img(i,j));
%     end
% end



% Mirror the image for mean/std computation
imgMirr = [img(:,params.windowSize + 1:-1:2) img img(:,end-1:-1:end-params.windowSize)];
imgMirr = [imgMirr(params.windowSize + 1:-1:2, :); imgMirr; imgMirr(end-1:-1:end-params.windowSize, :)];
% imgMirr = double(imgMirr);
%figure;imshow(imgMirr);
%title('mirrored image'); 

res = img;

% Estimate a global standard deviation value of the noise from the first
% rows of a image (useful for retinal OCT-Data)
if strcmp(params.sigmaMethod, 'firstrows')
    regionwidth = 10;
    temp = [reshape(img(1:regionwidth,:),1, size(img,2) * regionwidth)...
    reshape(img(end-regionwidth+1:end,:),1, size(img,2) * regionwidth)];
    %size(img, 2)
    

    sigma = std(double(temp(:)));
end

imgMean = zeros(size(img));
imgStd = zeros(size(img));

% Precompute mean and stdImages
for i = 1:size(img,1)
    for j = 1:size(img,2)
        region = imgMirr(i:i + 2 * params.windowSize, j:j + 2* params.windowSize);
        if(i==1 && j==1)
        size(region);
        mean(region);
        end
        imgMean(i,j) = mean(region(:));
        imgStd(i,j) = std(region(:));
    end
end

% Go over the image and perform the bayesian estimation of the pixel values
% without noise. See the original paper for details.
for i = 1:size(img,1)
    for j = 1:size(img,2)
       if  strcmp(params.sigmaMethod, 'local')
            sigma = imgStd(i,j);
        end
        compareMean = imgMean(i,j);
        
        weightSum = 0;
        intensitySum = 0;
        
        if sigma == 0
            weightSum = 1;
            intensitySum = img(i,j);
        else
            for k = 1:params.numSamples
                criteria = 2 * params.sigmaFactor * sigma;
                
                x = 1;
                y = 1;
                
                while criteria > params.sigmaFactor * sigma
                    x = round(randn(1) * params.sigmaSpatial) + i;
                    y = round(randn(1) * params.sigmaSpatial) + j;
                    
                    % Mirror the sampling points if necessary (border lines
                    % doubled). This has nothing to do with the actual paper
                    % implementation (as it is not mentioned), but seemed
                    % reasonable to us.
                    if x < 1
                        if mod(floor(-x / size(img,1)),2) == 0
                            x = abs(x + floor(-x / size(img,1)) * size(img,1)) + 1;
                        else
                            x = x + (floor(-x / size(img,1)) + 1) * size(img,1);
                        end
                    elseif x > size(img,1)
                        if mod(floor(x /size(img,1)),2) == 0
                            x = mod(x,size(img,1));
                        else
                            x = size(img,1) - mod(x,size(img,1)) + 1;
                        end
                    end
                    
                    if y < 1
                        if mod(floor(-y / size(img,2)),2) == 0
                            y = abs(y + floor(-y / size(img,2)) * size(img,2)) + 1;
                        else
                            y = y + (floor(-y / size(img,2)) + 1) * size(img,2);
                        end
                    elseif y > size(img,2)
                        if mod(floor(y /size(img,2)),2) == 0
                            y = mod(y,size(img,2));
                        else
                            y = size(img,2) - mod(y,size(img,2)) + 1;
                        end
                    end
                    
                    % Compute values out of the sampling sites position
                    siteMean = imgMean(x,y);
                    criteria = abs(compareMean - siteMean);
                end
                
                % Add weight to histogramm (histogram is implicit, as we
                % directly compute the weighted mean of the histogram)
                weight = exp(- (criteria / (2 * sigma * sigma)));
                weightSum = weightSum + weight;
                intensitySum = intensitySum + img(x,y) * weight;
            end
        end
        
        % Estimate noise-free value
        if weightSum == 0
            res(i,j) = 0;
        else
            res(i,j) = intensitySum / weightSum;
        end
    end
end
    
    for i = 1:size(res,1)
        for j = 1:size(res,2)
            res(i,j)=exp(res(i,j))-1;
            %new(i,j)=res(i,j)*255;
        end
    end
max(max(res));
end
