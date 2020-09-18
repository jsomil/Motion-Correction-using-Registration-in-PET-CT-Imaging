function [X, Y] = structensor2(A, f, type)
%STRUCTENSOR2 return orientation and anisotropy from the structure tensor.
%   [X, Y] = STRUCTENSOR2(A, f) calculates the local anisotropy in X, and
%   the local feature direction in Y for the image A, using the feature
%   scale given by the integer f.
%
%   [X, Y] = STRUCTENSOR2(A, f, 'classic') uses the classic definition of 
%   anisotropy as (1 - l2/l1), where l1 and l2 are the eigenvalues of the
%   structure tensor. Otherwise the default behaviour is to smooth l2
%   first, which improves the response at corners.
%
%   See also GRAD2, GAUSS2, ANISOTROPIC2, SVBITONIC2

%   Author: Graham M. Treece, University of Cambridge, UK, Sept 2018.


% check for appropriate inputs
narginchk(2, 3);

% also check for appropriate outputs
nargoutchk(2, 2);

% other input checks
if ~(ndims(A) == 2) && ~((ndims(A) == 3) && (size(A,3) == 3))
  error('A must be a grey or colour image.');
end
if nargin<3
  type = 'adapted';
else
  if ~strcmp(type, 'classic')
    type = 'adapted';
  end
end

% calculate gradients in X and Y dimensions, and subsequent components of
% the structure tensor
[X, Y] = grad2(A); % reverse X and Y to ensure orientation is correct later
Txx = X.^2;
Tyy = Y.^2;
Txy = 2*X.*Y;

% can now filter this tensor
Txx = gauss2(Txx, f);
Tyy = gauss2(Tyy, f);
Txy = gauss2(Txy, f);

% calculate anisotropy and orientation
T1 = Txx + Tyy;
T2 = Txx - Tyy;
Y = (atan2(Txy, T2) + pi)/2;
T2 = (T2.^2 + Txy.^2).^0.5;
X = zeros(size(Y));
n = (T1 + T2)>0;
if strcmp(type, 'classic')
  X(n) = 1 - (T1(n) - T2(n))./(T1(n) + T2(n));
else
  X(n) = 1 - gauss2((T1(n) - T2(n)),2*f)./(T1(n) + T2(n));
  X(X<0) = 0;
end

