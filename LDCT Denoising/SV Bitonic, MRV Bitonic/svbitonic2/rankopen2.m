function [X, X2] = rankopen2(A, M, R, t, dual)
%RANKOPEN2 robust opening and closing.
%   X = RANKOPEN2(A, M, R) performs a robust opening on A, using the 
%   mask in M, at the specified rank R. If R is 1 this is a normal
%   opening, and if R is the same as the number of non-zero elements
%   in the mask, this is a normal closing.
%
%   M is a mask, which must be either a logical or double l x l array,
%   where l is odd, and each true or non-zero element defines the
%   extent of the mask. This can have any shape but must be bitonic in
%   the horizontal and vertical directions, i.e. any line passing through
%   it in these directions only contains one non-zero segment.
%
%   The algorithm uses a fast histogram-based technique if A is uint8 or
%   int32, or a slightly less fast sorting-based technique for doubles.
%
%   X = RANKOPEN2(A, M, R, t) also uses t to determine a maximum
%   relative data range for a pixel to be included in the local masks.
%   t = 0 (the default) disables this.
%
%   X = RANKOPEN2(A, M, R, t, dual) uses 'dual' to determine whether
%   to perform a full opening (two rank filters, dual = true) or just a
%   robust erosion (one rank filter, dual = false).
%
%   [X, X2] = RANKOPEN2(A, M, R) outputs the closing operation in X2
%   as well as the opening in X. This is more efficient than using the
%   function twice.
%
%   See also BITONIC2, SVRANKOPEN2

%   Author: Graham M. Treece, University of Cambridge, UK, Sept 2018.
