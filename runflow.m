% run flow

% read the images
addpath('.\tools')

%%
% RubberWhale
frame1 = imread('data\other-data\Urban3\frame10.png');
frame2 = imread('data\other-data\Urban3\frame11.png');
gtflow = readFlowFile('data\other-gt-flow\Urban3\flow10.flo');

% % Grove2
% frame1 = imread('data\other-data\Grove2\frame10.png');
% frame2 = imread('data\other-data\Grove2\frame11.png');
% gtflow = readFlowFile('data\other-gt-flow\Grove2\flow10.flo');
% 
% % Grove3
% frame1 = imread('data\other-data\Grove3\frame10.png');
% frame2 = imread('data\other-data\Grove3\frame11.png');
% gtflow = readFlowFile('data\other-gt-flow\Grove3\flow10.flo');
% 
% % Hydrangea
% frame1 = imread('data\other-data\Hydrangea\frame10.png');
% frame2 = imread('data\other-data\Hydrangea\frame11.png');
% gtflow = readFlowFile('data\other-gt-flow\Hydrangea\flow10.flo');


% more data are included in the data folder

frame1 = double(rgb2gray(frame1));
frame2 = double(rgb2gray(frame2));
gtflow(gtflow>1e9) = NaN;
%%

% parameters
lambda = 10;  % you can try differnt value for different cases

% estimate optical flow using 
uv = estimateHSflow(frame1,frame2,lambda);

% warp frame2 to frame1 accroding to the estimated optical flow
[H W chs] = size(frame1);
[x,y] = meshgrid(1:W,1:H);
x1 = x+uv(:,:,1);
y1 = y+uv(:,:,2);
warpimg2 = [];
for ch = 1:size(frame2,3)
    warpimg2(:,:,ch) = interp2_bicubic(frame2(:,:,ch),x1,y1);
end
warpimg2(isnan(warpimg2)) = 0;

% compute the error with groundtruch
[aae stdae aepe] = flowAngErr(gtflow(:,:,1), gtflow(:,:,2), uv(:,:,1), uv(:,:,2), 0);

% print errors
fprintf('\nAAE %3.3f average EPE %3.3f \n', aae, aepe);

% show the results
figure,
subplot(2,2,1), imshow(uint8(frame1)),title('Frame1')
subplot(2,2,2), imshow(uint8(frame2)),title('Frame2')
subplot(2,2,3), imshow(flowToColor(uv)),title('Flow')
subplot(2,2,4), imshow(uint8(warpimg2)),title('Warped Image')









