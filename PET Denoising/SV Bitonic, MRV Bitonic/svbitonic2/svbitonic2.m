function [X] = svbitonic2(A, f, t, centile)
%SVBITONIC2 structurally varying bitonic filter for images.
%   X = SVBITONIC2(A, f) filters A with a set of masks within a region of
%   side length 2f+1.
%
%   The input A can be of either double, int32 or uint8 class. Processing
%   is considerably faster when using int32 or uint8, which is the default
%   for image data.
%
%   X = SVBITONIC2(A, f, t) also uses t to determine a maximum data range
%   for filtering This should be set to about 3.2 x the standard deviation
%   of the noise in the image. Default is t=0 which disables it.
%
%   X = SVBITONIC2(A, f, t, centile) uses centile as the minimum centile
%   for the morphological operations. The default (which is nearly always
%   the best choice) is 4.
%
%   Example:
%     A = imread('varying_noise.png','PNG');
%     B = svbitonic2(A, 6);
%     imagesc(B);
%     colormap(gray);
%
%   See also SVMASKS, SVRANKOPEN2, ANISOTROPIC2, MVBITONIC2, BITONIC2

%   Author: Graham M. Treece, University of Cambridge, UK, Sept 2018.


% check for appropriate inputs
narginchk(2, 4);
f = round(f);
if (nargin<3)
  t = 0;
end
if (nargin<4) 
  centile = 4;
end

% also check for appropriate outputs
nargoutchk(1, 1);

% check for compiled mex files
if (exist('svrankopen2','file')~=3)
  error('Type "mex svrankopen2.cpp" to compile this function.');
end
if (exist('anisotropic2_mex','file')~=3)
  error('Type "mex anisotropic2_mex.cpp" to compile this function.');
end

% Get set of masks and rank values
[M, Shapes, Angles, R] = svmasks(f, centile);

% Get initial anisotropy and mask directions
if size(A,3) == 1
  [Anis, Dir] = structensor2(A, f);
  [Aopen, Mopen, Aclose, Mclose] = svrankopen2(A, M, R, t, false, true);
else
  grey = mean(A,3);
  [Anis, Dir] = structensor2(grey, f);
  if strcmp(class(A), 'uint8')
    grey = uint8(grey);
  elseif strcmp(class(A), 'int32')
    grey = int32(grey);
  end
  [Aopen, Mopen, Aclose, Mclose] = svrankopen2(grey, M, R, t, false, true);
end

% some additional information relating to the mask set
Mshape = [];
for s=1:Shapes
  Mshape = [Mshape; s*ones(Angles(s), 1)];
end
Mstart = [0; cumsum(Angles(1:(end-1)))];

% Set noise floor on anisotropy
Amin = 5.5 / (f + 5) - 0.05;

% Option for mask type; 1-from anisotropy, 2-mixture, 3-from morphology
mask_type = 2;

% Use anisotropy to refine mask sets
if (mask_type ~= 3)
  s = round((1 - (Anis / 0.8)) * (Shapes - 1)) + 1;
  s(s<1) = 1;
  s(s>Shapes) = Shapes;
  a = round((Angles(s)-1).*(Dir / pi)) + 1;
  a(a>Angles(s)) = 1;
  n = (s < Mshape(Mopen)) & (Anis > Amin);
  if (mask_type == 1)
    Mopen = Mstart(s) + a;
  elseif (mask_type == 2)
    Mopen(n) = Mstart(s(n)) + a(n);
  end
  n = (s < Mshape(Mclose)) & (Anis > Amin);
  if (mask_type == 1)
    Mclose = Mstart(s) + a;
  elseif (mask_type == 2)
    Mclose(n) = Mstart(s(n)) + a(n);
  end
end

% Now work on each channel separately
X = zeros(size(A));
for k=1:size(A,3)

  % Structurally varying opening and closing with these mask sets
  [Aopen, MRopen, Aclose, MRclose] = svrankopen2(A(:,:,k), M, R, t, true, true, Mopen, Mclose);
  
  % Use new reverse masks to improve anisotropy
  if (mask_type == 2)
    s = Mshape(MRopen);
    Anis1 = 0.8 * (1 - (s - 1) / (Shapes - 1));
    Dir1 = (double(MRopen) - Mstart(s)) * pi ./ Angles(s);
    s = Mshape(MRclose);
    Anis2 = 0.8 * (1 - (s - 1) / (Shapes - 1));
    Dir2 = (double(MRclose) - Mstart(s)) * pi ./ Angles(s);
    n = (Anis < Amin) & (Anis1 < Amin) & (Anis2 < Amin);
    n1 = Anis1(n)>Anis2(n);
    Anis(n) = Anis1(n).*n1 + Anis2(n).*(~n1);
    Dir(n) = Dir1(n).*n1 + Dir2(n).*(~n1);
  end
  
  % smooth errors
  Eopen = anisotropic2(double(A(:,:,k) - Aopen), f, 0.6, 1.5*t, Anis, Dir);
  Eclose = anisotropic2(double(Aclose - A(:,:,k)), f, 0.6, 1.5*t, Anis, Dir);
  
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