function uv = estimateHSflow(frame1,frame2,lambda)
% estimate the HS flow

if ~exist('lambda','var')
    labmda = 80;
end
[H, W, chs] = size(frame1);

% build the image pyramid
pyramid_spacing = 2.0;
pyramid_levels  = 1 + floor( log(min(W,H)/16) / log(pyramid_spacing) );
smooth_sigma = sqrt(2);   
f            = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma);

pyramid1 = cell(pyramid_levels,1);
pyramid2 = cell(pyramid_levels,1);
pyramid1{1} = frame1;
pyramid2{1} = frame2;

for m=2:pyramid_levels
   % TODO #1: build Gaussian pyramid for coarse-to-fine optical flow estimation
   pyramid1{m-1} = imfilter(pyramid1{m-1},f);
   pyramid2{m-1} = imfilter(pyramid2{m-1},f);
   pyramid1{m} = imresize(pyramid1{m-1},0.8);
   pyramid2{m} = imresize(pyramid2{m-1},0.8);
   % END TODO
end

% coarst-to-fine compute the flow
uv = zeros(H,W,2);
for levels = pyramid_levels:-1:1
    fprintf('level: %d \n',levels);
    [H1, W1, chs] = size(pyramid1{levels});
    uv = resample_flow(uv,[H1,W1]);
    uv = estimateHSflowlayer(pyramid1{levels},pyramid2{levels},uv,lambda);
end
end