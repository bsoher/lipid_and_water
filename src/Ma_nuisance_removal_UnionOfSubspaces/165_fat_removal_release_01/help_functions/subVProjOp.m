function y = subVProjOp( x, designPara)
% Projection operator for updating V with given U

% prepare parameters
if isfield(designPara,'lambda')
	lambda    = designPara.lambda;
else
    lambda    = 0;
end
if isfield(designPara,'BigLambda')
	BigLambda = diag(designPara.BigLambda);
else
	BigLambda = 1;
end

% AcAx
Ax            = subVAOp( x, 'notransp', designPara);
AcAx          = subVAOp(Ax, 'transp', designPara);

% reg1
reg1          = lambda*((BigLambda')*BigLambda*reshape(x,size(BigLambda,2),[]));

% cost function
y             = AcAx + reg1(:);


end