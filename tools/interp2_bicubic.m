function [ZI] = interp2_bicubic(Z, XI, YI, Dxfilter)
  if nargin < 4
      Dxfilter = [1 -8 0 8 -1]/12;
  end;
  Dyfilter = Dxfilter';
  Dxyfilter = conv2(Dxfilter, Dyfilter, 'full');
  
  input_size = size(XI);
  
  % Reshape input coordinates into a vector
  XI = reshape(XI, 1, prod(input_size));
  YI = reshape(YI, 1, prod(input_size));
  
  % Bound coordinates to valid region
  sx = size(Z, 2);
  sy = size(Z, 1);
  
  fXI = floor(XI);
  cXI = fXI + 1;
  fYI = floor(YI);
  cYI = fYI + 1;
 
  indx = (fXI<1) | (cXI>sx) | (fYI<1) | (cYI>sy);
  
  fXI = max(1, min(sx, fXI));
  cXI = max(1, min(sx, cXI));
  fYI = max(1, min(sy, fYI));
  cYI = max(1, min(sy, cYI));

  
  % Image at 4 neighbors
  Z00 = Z(fYI + sy * (fXI - 1));
  Z01 = Z(cYI + sy * (fXI - 1));
  Z10 = Z(fYI + sy * (cXI - 1));
  Z11 = Z(cYI + sy * (cXI - 1));
  
  % x-derivative at 4 neighbors
%   DX = imfilter(Z, [-0.5, 0, 0.5], 'symmetric', 'corr');
  DX = imfilter(Z, Dxfilter, 'symmetric', 'corr'); % 
  DX00 = DX(fYI + sy * (fXI - 1));
  DX01 = DX(cYI + sy * (fXI - 1));
  DX10 = DX(fYI + sy * (cXI - 1));
  DX11 = DX(cYI + sy * (cXI - 1));

  % y-derivative at 4 neighbors
%   DY = imfilter(Z, [-0.5, 0, 0.5]', 'symmetric', 'corr'); 
  DY = imfilter(Z, Dyfilter, 'symmetric', 'corr'); % 
  DY00 = DY(fYI + sy * (fXI - 1));
  DY01 = DY(cYI + sy * (fXI - 1));
  DY10 = DY(fYI + sy * (cXI - 1));
  DY11 = DY(cYI + sy * (cXI - 1));

  % xy-derivative at 4 neighbors
%   DXY = imfilter(Z, 0.25 * [1, 0, -1; 0 0 0; -1 0 1], 'symmetric', 'corr');
  DXY = imfilter(Z, Dxyfilter, 'symmetric', 'corr'); % modified by dqsun
  DXY00 = DXY(fYI + sy * (fXI - 1));
  DXY01 = DXY(cYI + sy * (fXI - 1));
  DXY10 = DXY(fYI + sy * (cXI - 1));
  DXY11 = DXY(cYI + sy * (cXI - 1));

  
  W = [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;
        0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0;
       -3,  0,  0,  3,  0,  0,  0,  0, -2,  0,  0, -1,  0,  0,  0,  0;  
        2,  0,  0, -2,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,  0;
        0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0;
        0,  0,  0,  0, -3,  0,  0,  3,  0,  0,  0,  0, -2,  0,  0, -1;
        0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0,  0,  1,  0,  0,  1;
       -3,  3,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;
        0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -2, -1,  0,  0;
        9, -9,  9, -9,  6,  3, -3, -6,  6, -6, -3,  3,  4,  2,  1,  2;
       -6,  6, -6,  6, -4, -2,  2,  4, -3,  3,  3, -3, -2, -1, -1, -2;
        2, -2,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;
        0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  1,  1,  0,  0;
       -6,  6, -6,  6, -3, -3,  3,  3, -4,  4,  2, -2, -2, -2, -1, -1;
        4, -4,  4, -4,  2,  2, -2, -2,  2, -2, -2,  2,  1,  1,  1,  1];

  
  V = [Z00(:)'; Z10(:)'; Z11(:)'; Z01(:)'; ...
       DX00(:)'; DX10(:)'; DX11(:)'; DX01(:)'; ...
       DY00(:)'; DY10(:)'; DY11(:)'; DY01(:)'; ...
       DXY00(:)'; DXY10(:)'; DXY11(:)'; DXY01(:)'];
  
  C = W * V;
  
  alpha_x = reshape(XI - fXI, input_size);
  alpha_y = reshape(YI - fYI, input_size);
  
  % Clip out-of-boundary pixels to boundary
  alpha_x(indx) = 0;
  alpha_y(indx) = 0;

  
  
  % Interpolation

  ZI = zeros(input_size);
  
  idx = 1;
  for i = 0:3
    for j = 0:3
      ZI = ZI + reshape(C(idx, :), input_size) .* alpha_x.^i .* ...
           alpha_y.^j;
      idx = idx + 1;
    end
  end  

ZI(indx) = NaN;
end
