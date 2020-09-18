function [X, Shapes, Angles, R] = svmasks(f, centile)
%SVMASKS creates a set of masks for structurally varying morphology
%   X = SVMASKS(f) creates a set of elliptical masks of varying aspect
%   ratio, appropriate for a window width of 2f+1. X is a 3D logical array
%   containing these masks.
%
%   [X, Shapes, Angles] = SVMASKS(f, centile) also returns the number of
%   different shapes and number of angles in each of these shapes.
%
%   [X, Shapes, Angles, R] = SVMASKS(f, centile) also calculates
%   an appropriate set of ranking values for each mask, given an initial
%   minimum centile. R are the ranking values necessary for performing
%   structurally varying opening and closing.
%
%   See also SVRANKOPEN2, SVBITONIC2

%   Author: Graham M. Treece, University of Cambridge, UK, Sept 2018.


% check for appropriate inputs
narginchk(1, 2);

% also check for appropriate outputs
nargoutchk(1, 5);

% other input checks
if (nargout > 3) && (nargin < 2)
  error('You must define input centile if output ranks are required.');
end
f = round(f);

% set up a few basic parameters
width = 2*f + 1;
stotal = round(f/4)+1;
if (stotal < 2) stotal = 2; end
rx = f + 0.25;
x = ones(width,1)*[-f:f];
y = [-f:f]'*ones(1,width);

% work out number of masks, and ratios
masks = 0;
ratio = ones(stotal, 1);
atotal = ones(stotal, 1);
for s=1:stotal
  if (s==stotal)
    ratio(s) = 1.0;
    atotal(s) = 1;
  else
    ratio(s) = 1.7 / (f + 2.0) + ((f + 0.3) / (f + 2.0)) * (s - 1) / (stotal - 1);
    atotal(s) = 4 * round((1-ratio(s))*(f + 2)*0.4);
  end
  masks = masks + atotal(s);
end

% create the masks
p = ones(masks, 1);
sh = ones(masks, 1);
X = false(width, width, masks);
m = 1;
for s=1:stotal
  ry = rx * ratio(s);
  
  for a = 1:atotal(s)
    sh(m) = s;

    ang = pi*(a-1)/atotal(s);
    xa = x*cos(ang) + y*sin(ang);
    ya = -x*sin(ang) + y*cos(ang);
    xa((xa<=0.25) & (xa>=-0.25)) = 0;
    xa(xa>0.25) = xa(xa>0.25)-0.25;
    xa(xa<-0.25) = xa(xa<-0.25)+0.25;
    ya((ya<=0.25) & (ya>=-0.25)) = 0;
    ya(ya>0.25) = ya(ya>0.25)-0.25;
    ya(ya<-0.25) = ya(ya<-0.25)+0.25;
    
    X(:,:,m) = X(:,:,m) | ((xa/rx).^2 + (ya/ry).^2 <= 1.0);
    
    p(m) = sum(sum(X(:,:,m)));

    m = m + 1;
  end
end

% return additional information
Shapes = stotal;
Angles = atotal;

% Calculate the appropriate ranking values
% Note that Rmin and Rmax start at 1
if (nargout>3)
  p1 = sum(p(sh==1))/atotal(1);
  po = p(masks);
  g = centile/100;
  pn = 2.0*(1 - (p1./p).^0.5)*((p1./po).^0.25 / (1 - (p1./po).^0.5));
  goff = pn * 0.01 * (3 + 4.7*(1 - exp(-8 * g)));
  pn = 25.0./po.^0.5;
  goff = goff + pn.* (0.01 * (3.0 + 3.3*(1 - exp(-10 * g)))*((sh-1) / (stotal-1)));
  R = round((g + goff).*(p-1)) + 1;
  
  % Limit very small or large values
  if (centile > 0.0)
    R(R<2) = 2;
  end
  n = R>floor(p/2);
  R(n) = floor(p(n)/2);
end
