function [X, Mout, X2, M2out] = svrankopen2(A, M, R, t, dual, masks, Min, M2in)
%SVRANKOPEN2 robust structurally varying opening and closing.
%   X = SVRANKOPEN2(A, M, R) performs a robust structurally varying opening
%   on A, using the set of masks in M, at the corresponding ranks R. Masks
%   are chosen from the provided set according to which best fits the data
%   at each point. Suitable sets of M and R can be calculated using the
%   function SVMASKS.
%
%   M is a mask, which must be either a logical or double l x l x m array,
%   where l is odd, and m is the total numer of masks. Each true or
%   non-zero element in each l x l slice defines the extent of the mask.
%   This can have any shape but must be bitonic in the horizontal and
%   vertical directions, i.e. any line passing through it in these
%   directions only contains one non-zero segment.
%
%   R defines the rank for each of these masks, and must be a vector
%   of length m.
%
%   The algorithm uses a fast histogram-based technique if A is uint8 or
%   int32, or a slightly less fast sorting-based technique for doubles.
%
%   X = SVRANKOPEN2(A, M, R, t) also uses t to determine a maximum
%   relative data range for a pixel to be included in the local masks.
%   t = 0 (the default) disables this.
%
%   X = SVRANKOPEN2(A, M, R, t, dual) uses 'dual' to determine whether
%   to perform a full opening (two rank filters, dual = true) or just a
%   robust erosion (one rank filter, dual = false).
%
%   [X, Mout] = SVRANKOPEN2(A, M, R, T, dual, maskout) uses 'maskout' to
%   determine whether to output (maskout = true) an array Mout containing
%   the mask index actually used for each pixel in the input data A. The
%   default is no output (maskout = false).
%
%   X = SVRANKOPEN2(A, M, R, T, dual, maskout, Min) uses Min to set which
%   mask should be used for each pixel in A. Min should therefore be an
%   array of the same size as A, containing numbers between 1 and m.
%
%   [X, X2] = RANKOPEN2(A, M, R) outputs the closing operation in X2
%   as well as the opening in X1. This is more efficient than using the
%   function twice.
%
%   [X, Mout, X2, M2out] = SVRANKOPEN2(A, M, R, T, dual, maskout) will
%   output mask arrays for both operations, if maskout = true.
%
%   See also SVMASKS, SVBITONIC2, RANKOPEN2

%   Author: Graham M. Treece, University of Cambridge, UK, Sept 2018.

