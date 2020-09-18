function [X] = gauss2(A, f)
%GAUSS2 smooth image using a Gaussian filter.
%   X = GAUSS2(A, f) filters A with an integer extent given by f.
%   The window length is 2.6f+1, and the standard deviation of the
%   Gaussian is set to 0.65f. This makes the smoothing roughly equivalent
%   to a mean filter of length 2f+1. The filter is slightly modified from
%   a Gaussian for f = 1 or 2, to improve the performance.
%
%   If A is of integer type, the returned value is uint8.
%
%   See also ANISOTROPIC2, STRUCTENSOR2

%   Author: Graham M. Treece, University of Cambridge, UK, Sept 2018.


% check for appropriate inputs
narginchk(2, 2);

% also check for appropriate outputs
nargoutchk(1, 1);

% other input checks
if ~(ndims(A) == 2) && ~((ndims(A) == 3) && (size(A,3) == 3))
  error('A must be a grey or colour image.');
end

% set up a gaussian linear filter of appropriate length
gf = floor(2*0.65*f);
gn = 2*gf+1;
gx = [0:(gn-1)]-gf;
gauss = exp(-gx.^2/(0.65*f)^2);

% Adjust slightly for very short filters
if (gf == 0)
  X = A;
  return
elseif (gf == 1)
  gauss(1) = 0.75 + 0.25*gauss(1);
  gauss(3) = 0.75 + 0.25*gauss(3);
elseif (gf == 2)
  gauss(1) = 0.1 + 0.9*gauss(1);
  gauss(2) = 0.25 + 0.75*gauss(2);
  gauss(4) = 0.25 + 0.75*gauss(4);
  gauss(5) = 0.1 + 0.9*gauss(5);
end
gauss = gauss/sum(gauss);

% get extension vectors for filtering beyond image edges
m = [ones(1,(length(gauss)-1)/2) 1:size(A,1) size(A,1)*ones(1,(length(gauss)-1)/2)];
n = [ones(1,(length(gauss)-1)/2) 1:size(A,2) size(A,2)*ones(1,(length(gauss)-1)/2)];

% perform convolution
if (size(A,3)==1)
  X = conv2(gauss, gauss, double(A(m,n)), 'valid');
else
  for k=1:size(A,3)
    X(:,:,k) = conv2(gauss, gauss, double(A(m,n,k)), 'valid');
  end
end

% ensure the right data type
if strcmp(class(A), 'uint8')
  X = uint8(X);
elseif strcmp(class(A), 'int32')
  X = int32(X);
end