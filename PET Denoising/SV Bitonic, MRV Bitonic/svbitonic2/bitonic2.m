function [X] = bitonic2(A, f, t, centile)
%BITONIC2 bitonic filter for images.
%   X = BITONIC2(A, f) filters A with a circular masks of diameter f.
%
%   The input A can be of either double, int32 or uint8 class. Processing
%   is considerably faster when using int32 or uint8, which is the default
%   for image data.
%
%   X = BITONIC2(A, f, t) also uses t to determine a maximum data range
%   for filtering This should be set to about 3.2 x the standard deviation
%   of the noise in the image. Default is t=0 which disables it.
%
%   X = BITONIC2(A, f, t, centile) uses centile as the minimum centile
%   for the morphological operations. The default (which is usually
%   the best choice) is 10.
%
%   Example:
%     A = imread('varying_noise.png','PNG');
%     B = bitonic2(A, 6);
%     imagesc(B);
%     colormap(gray);
%
%   See also RANKOPEN2, GAUSS2, SVBITONIC2

%   Author: Graham M. Treece, University of Cambridge, UK, Sept 2018.


% check for appropriate inputs
narginchk(2, 4);
f = round(f);
if (nargin<3) 
  t = 0;
end
if (nargin<4) 
  centile = 10;
end

% also check for appropriate outputs
nargoutchk(1, 1);

% check for compiled mex file
if (exist('rankopen2','file')~=3)
  error('Type "mex rankopen2.cpp" to compile this function.');
end

% Get mask and convert centile to rank value
x = [0:(2*f)]-f;
n = 2*f+1;
M = double(((x'*ones(1,n)).^2 + (ones(n,1)*x).^2)<(f+0.1)^2);
p = sum(M(:));
c1 = round((centile/100)*(p-1)) + 1;
c2 = p + 1 - c1;

% run over all colours if RGB data
for k=1:size(A,3)
  
  % perform robust opening and closing
  [Aopen, Aclose] = rankopen2(A(:,:,k), M, c1, t);
  
  % Gaussian filter on difference of these from original
  Eopen = gauss2(double(A(:,:,k)-Aopen), f);
  Eclose = gauss2(double(Aclose-A(:,:,k)), f);
  
  % Remove offsets from Aopen and Aclose
  Aclose = double(Aclose) - Eclose;
  Aopen = double(Aopen) + Eopen;
  
  % Increase weighting factor
  Eopen = abs(Eopen.^3);
  Eclose = abs(Eclose.^3);
  
  % find sum of weighting values (output of gaussian filters), remove zeros
  Esum = Eopen+Eclose;
  n = (Esum==0);
  Eopen(n) = 0.5;
  Eclose(n) = 0.5;
  Esum(n) = 1.0;
  
  % output is weighted sum of ordinal operations
  X(:,:,k) = (Eopen.*Aclose + Eclose.*Aopen)./Esum;
  
end

% ensure the right data type
if strcmp(class(A), 'uint8')
  X = uint8(X);
elseif strcmp(class(A), 'int32')
  X = int32(X);
end