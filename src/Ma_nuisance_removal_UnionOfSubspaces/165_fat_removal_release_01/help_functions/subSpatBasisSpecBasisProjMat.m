
function y = subSpatBasisSpecBasisProjMat( x,Ax,At,lambda,W)

if nargin < 4
    lambda = 0;
end

if nargin < 5
    W= [];
end

x1 = subSpatBasisSpecBasisMat( x, Ax,At,'notransp',W);
x2 = subSpatBasisSpecBasisMat( x1, Ax,At, 'transpose',W);

y = x2 + lambda*x;
