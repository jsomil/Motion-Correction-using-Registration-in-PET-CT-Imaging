function [X, Y] = grad2(A)
%GRAD2 find image gradients.
%   [X, Y] = GRAD2(A) calculates the gradients in A, in the X and Y
%   directions, using a Sobel filter. They are always returned as doubles.
%
%   See also STRUCTENSOR2

%   Author: Graham M. Treece, University of Cambridge, UK, Sept 2018.


% check for appropriate inputs
narginchk(1, 1);

% also check for appropriate outputs
nargoutchk(2, 2);

% other input checks
if ~(ndims(A) == 2) && ~((ndims(A) == 3) && (size(A,3) == 3))
  error('A must be a grey or colour image.');
end

% set up filters in each dimension
gf = [-1 0 1]/2;
sf = [1 2 1]/4;

% get extension vectors for filtering beyond image edges
m = [1 1:size(A,1) size(A,1)];
n = [1 1:size(A,2) size(A,2)];

% perform convolution
if (size(A,3)==1)
  X = conv2(sf, gf, double(A(m,n)), 'valid');
  Y = conv2(gf, sf, double(A(m,n)), 'valid');
else
  for k=1:size(A,3)
    X(:,:,k) = conv2(sf, gf, double(A(m,n,k)), 'valid');
    Y(:,:,k) = conv2(gf, sf, double(A(m,n,k)), 'valid');
  end
end
