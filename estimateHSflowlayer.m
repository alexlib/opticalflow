function uv = estimateHSflowlayer(frame1,frame2,uv,lambda, maxwarping)
% compute Horn Shark optical flow

if ~exist('lambda','var')
    lambda = 80;
end

if ~exist('maxwarping','var')
    maxwarping = 10;
end

H = size(frame1,1);
W = size(frame1,2);
npixels = H*W;
[x, y] = meshgrid(1:W,1:H);

%TODO#3: build differential matrix and Laplacian matrix according to image size
e = ones(npixels, 1);
dy = spdiags([-e e],0:1,npixels,npixels);
dx = spdiags([-e e],[0, H],npixels,npixels);
L = dx.'*dx+dy.'*dy;
% END TODO

% Kernel to get gradient
h = [1 -8 0 8 -1]/12;

for i=1:maxwarping
    % TODO#2: warp image using the flow vector
    x1 = x+uv(:,:,1);
    y1 = y+uv(:,:,2);
    warpimg2 = [];
    warpimg2 = interp2_bicubic(frame2,x1,y1);
    % END TODO
    warpimg2(isnan(warpimg2)) = 0;
    
    % TODO#4: compute image gradient Ix, Iy, and Iz 
    Iz = warpimg2 - frame1;
    Ix = imfilter(warpimg2,h);
    Iy = imfilter(warpimg2,h.');
    % END TODO
    
    % TODO#5: build linear system to solve HS flow
    Ix = reshape(Ix,npixels,1);
    Iy = reshape(Iy,npixels,1);
    Iz = reshape(Iz,npixels,1);
    U = reshape(uv(:,:,1),npixels,1);
    V = reshape(uv(:,:,2),npixels,1);
    
    Ix = spdiags(Ix,0,npixels,npixels);
    Iy = spdiags(Iy,0,npixels,npixels);
    
    A = [Ix*Ix+lambda*L Ix*Iy; Ix*Iy Iy*Iy+lambda*L];
    b = -[Ix*Iz+lambda*L*U; Iy*Iz+lambda*L*V];
    % solve linear system
    deltauv     = reshape(A\b, size(uv));
%     [xx ~] = pcg(A,b,[],100);
%     deltauv     = reshape(xx, size(uv));
    
    deltauv(deltauv>1) = 1;
    deltauv(deltauv<-1) = -1;
    
    uv = uv+deltauv;
    
    % TODO#6: use median filter to smooth the flow map
    uv(:,:,1) = medfilt2(uv(:,:,1),[7 7]);
    uv(:,:,2) = medfilt2(uv(:,:,2),[7 7]);
    
    fprintf('Warping step: %d, Incremental norm: %3.5f \n', i, norm(deltauv(:)));
end

end
