function Y = normcols(X)
%NORMCOLS  Scale each column of X to unit Euclidean norm.
%  Y = NORMCOLS(X) returns X with every column divided by its 2-norm.
%
%  This reproduces normc from the Deep Learning Toolbox bit for bit, and is
%  provided so that the package runs on base MATLAB alone.

    Y = X ./ sqrt(sum(X.^2, 1));
end
