function [X] = mvbitonic2(A, f, t, levels, centile, extra)
%MVBITONIC2 multi-resolution structurally varying bitonic filter.
%   X = MVBITONIC2(A, f) filters A with a set of masks within a region of
%   side length 2f+1.
%
%   The input A can be of either double, int32 or uint8 class. Processing
%   is considerably faster when using int32 or uint8, which is the default
%   for image data.
%
%   X = MVBITONIC2(A, f, t) also uses t to determine a maximum data range
%   for filtering This should be set to about 3.2 x the standard deviation
%   of the noise in the image. Default is t=0 which disables it.
%
%   X = MVBITONIC2(A, f, t, levels) sets the maximum number of levels in
%   the multi-resolution processing. The default is levels = 3.
%
%   X = MVBITONIC2(A, f, t, levels, centile) uses centile as the minimum
%   centile for the morphological operations. The default (which is nearly
%   always the best choice) is 4.
%
%   X = MVBITONIC2(A, f, t, levels, centile, extra) sets whether to apply
%   the svbitonic an additional time at the end of each level. The default
%   behaviour is to do so (extra = true). Setting extra = false halves
%   processing time, but with slightly poorer performance.
%
%   Example:
%     A = imread('varying_noise.png','PNG');
%     B = mvbitonic2(A, 6, 180);
%     imagesc(B);
%     colormap(gray);
%
%   See also SVMASKS, SVBITONIC2, PROLONGATE2, RESTRICT2

%   Author: Graham M. Treece, University of Cambridge, UK, Sept 2018.


% check for appropriate inputs
narginchk(2, 6);
f = round(f);
if (nargin<3)
  t = max(A(:))-min(A(:));
end
if (nargin<4)
  levels = 3;
end
if (nargin<5) 
  centile = 4;
end
if (nargin<6)
  extra = true;
end

% also check for appropriate outputs
nargoutchk(1, 1);

% check for compiled mex file
if (exist('svrankopen2','file')~=3)
  error('Type "mex svrankopen2.cpp" to compile this function.');
end
if (exist('anisotropic2_mex','file')~=3)
  error('Type "mex anisotropic2_mex.cpp" to compile this function.');
end

% Work out thresholds for all levels, and limit max levels
thresh = double(t)./((2.4*f).^[0:(levels-1)]);
if ~strcmp(class(A), 'double')
  thresh = thresh(thresh>=1);
  levels = min([levels length(thresh)]);
end

% initial input at first level is A
% need to use int32 rather than uint
if strcmp(class(A), 'uint8')
  B{1} = int32(A);
else
  B{1} = A;
end

% Go down through levels, processing and restricting
for l=1:levels
  
  % Apply bitonic at this level
  B2{l} = svbitonic2(B{l}, f, thresh(l), centile);
  
  % possibly restrict and correct ready for deeper levels
  if (l < levels)
    B{l+1} = restrict2(B2{l});
    B2{l} = B2{l} - prolongate2(B{l+1}, size(B2{l}));
  end
  
end

% Go up through levels, prolongating and adding result
for l=(levels-1):-1:1
  
  % add in prolongated lower level result
  B2{l} = B2{l} + prolongate2(B2{l+1}, size(B2{l}));

  % possibly additional bitonic with lower threshold
  if extra
    B2{l} = svbitonic2(B2{l}, f, thresh(l+1), centile);
  end

end

% result is processed output at top level
if strcmp(class(A), 'uint8')
  X = uint8(B2{1});
else
  X = B2{1};
end
  