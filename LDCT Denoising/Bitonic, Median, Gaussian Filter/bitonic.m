function [X, Y, Z] = bitonic(A, f, o)
%BITONIC preserve bitonicity in the input data
%   X = BITONIC(A, f) filters A with an integer extent given by f. In 1D
%   the window length is 2f+1, in 2D the window is a disc of radius f. A
%   can either be a vector, a 2D array or an RGB image als a 3D array. In
%   2D, the filter length defines the radius of the disc used in
%   morphological operations.
%
%   X = BITONIC(A, f, o) uses o as the base centile for the
%   morphological operations. The default (which is nearly always the best
%   value) is 10. Lower values can preserve isolated fine details at the cost of
%   some slight nonlinear distortion. A value of 20 or higher is effective
%   at reducing salt & pepper noise, but will also eliminate more isolated details.
%  
%   [X, Y, Z] = BITONIC(A, f) will also output a Gaussian filtered image
%   in Y, and a median filtered image in Z, using the same size filters and
%   morphological structural elements.
%
%   The input A can be of either double or uint8 class. If the rankfilt2
%   function has been compiled, processing is considerably faster when
%   using uint8, which is the default for image data.
%
%   Example:
%     A = imread('varying_noise.png','PNG');
%     B = bitonic(A, 5);
%     imagesc(B);
%     colormap(gray);

% check for appropriate inputs
narginchk(2, 3);
if (nargin<3) 
  o = 0.1;
else
  o = o/100;
end
% r = 0.65;
r = 3;

% also check for appropriate outputs
nargoutchk(1, 3);

% set up a gaussian linear filter of appropriate length
gf = floor(r*2*f);
gn = 2*gf+1;
gx = [0:(gn-1)]-gf;
gauss = exp(-gx.^2/(r*f)^2);
rgauss = gauss/sum(gauss);
if (gf == 1)
  % slight adjustment for very short filters
  gauss = 1-(1-gauss)/4;
end
gauss = gauss/sum(gauss);

% set up shape for rank filter and centiles
if isvector(A)
  n = 2*f+1;
  if iscolumn(A)
    shape = ones(n,1);
  else
    shape = ones(1,n);
  end
  c1 = round(o*(n-1)) + 1;
  c2 = n+1-c1;
else
  x = [0:(2*f)]-f;
  n = 2*f+1;
  shape = double(((x'*ones(1,n)).^2+(ones(n,1)*x).^2)<(f+0.1)^2);
  n = sum(sum(shape));
  f = (n-1)/2;
  c1 = round(o*(n-1)) + 1;
  c2 = n+1-c1;
end

% get extension vectors for filtering beyond image edges
m = [ones(1,(length(gauss)-1)/2) 1:size(A,1) size(A,1)*ones(1,(length(gauss)-1)/2)];
n = [ones(1,(length(gauss)-1)/2) 1:size(A,2) size(A,2)*ones(1,(length(gauss)-1)/2)];

% check for existence of optimised mex file
% or ordfilt2 from image processing toolbox
if (exist('rankfilt2','file')==3)
  mex_file = true;
else
  disp(sprintf('This function is more efficient if you compile the rankfilt2 function.\nParticularly for data of type uint8.\nType "mex rankfilt2.cpp" to do so.')); 
  mex_file = false;
end

% run over all colours if RGB data
for k=1:size(A,3)
  
  % perform robust opening and closing
  if (mex_file)
    [Aopen, Aclose] = rankfilt2(A(:,:,k), c1, shape, c2);
    Aopen = rankfilt2(Aopen, c2, shape);
    Aclose = rankfilt2(Aclose, c1, shape);
  else
    Aopen = ordfilt2(ordfilt2(A(:,:,k), c1, shape, 'symmetric'), c2, shape, 'symmetric');
    Aclose = ordfilt2(ordfilt2(A(:,:,k), c2, shape, 'symmetric'), c1, shape, 'symmetric');
  end
  
  % Gaussian filter on difference of these from original
  if iscolumn(A)
    Eopen = abs(conv(double(A(m,:,k)-Aopen(m,:))', gauss, 'valid'))';
    Eclose = abs(conv(double(Aclose(m,:)-A(m,:,k))', gauss, 'valid'))';
  elseif isrow(A)
    Eopen = abs(conv(double(A(:,n,k)-Aopen(:,n)), gauss, 'valid'));
    Eclose = abs(conv(double(Aclose(:,n)-A(:,n,k)), gauss, 'valid'));
  else
    Eopen = abs(conv2(gauss, gauss, double(A(m,n,k)-Aopen(m,n)), 'valid'));
    Eclose = abs(conv2(gauss, gauss, double(Aclose(m,n)-A(m,n,k)), 'valid'));
  end
  
  % find sum of weighting values (output of gaussian filters), remove zeros
  Esum = Eopen+Eclose;
  Eopen(Esum==0) = 0.5;
  Eclose(Esum==0) = 0.5;
  Esum(Esum==0) = 1.0;
  
  % output is weighted sum of ordinal operations
  X(:,:,k) = (Eopen.*double(Aclose) + Eclose.*double(Aopen))./Esum;
  
  % might need to output a gaussian or median with same filter lengths
  if (nargout>1)
    if iscolumn(A)
      Y(:,:,k) = conv(double(A(m,:,k))', rgauss, 'valid')';
    elseif isrow(A)
      Y(:,:,k) = conv(double(A(:,n,k)), rgauss, 'valid');
    else
      Y(:,:,k) = conv2(rgauss, rgauss, double(A(m,n,k)), 'valid');
    end
  end
  if (nargout>2)
    if (mex_file)
      Z(:,:,k) = rankfilt2(A(:,:,k), f+1, shape);
    else
      Z(:,:,k) = ordfilt2(A(:,:,k), f+1, shape, 'symmetric');
    end
  end
  
end

return;