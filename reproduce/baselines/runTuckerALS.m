function T = runTuckerALS(filename, szMode, rank, epsilon, maxIter)
addpath([pwd '/tensor_toolbox_2.6']);
addpath([pwd '/tensor_toolbox_2.6/met']);
%cd tensor_toolbox_2.5; addpath(pwd);
%cd met; addpath(pwd);
%cd ..;
%cd ..;

%filename = '../synData/syn5_100_50000';
%% filename = 'stack_overflow.tensor.top1000';
%szMode = 5;

%%
filename
X = dlmread(filename,' ');
subs = X(:, 1:szMode); %<-- Subscripts of the nonzeros.
vals = X(:, szMode+1); %<-- The values of the nonzeros.
subs = subs + 1;
X = sptensor(subs,vals);

%%
opts = struct;
%initialize factor matrices using HOSVD
%opts.init = 'eigs'
opts.maxiters = maxIter;

tic
%T = tucker_me(X,rank,1,opts);
T = tucker_als(X, rank, 'tol', epsilon, 'maxiters', maxIter);
toc

%%
% model = [T.U{1}; T.U{2}; T.U{3}];
% model2 = dlmread(modelname, ',');
%fprintf('Norm: %f\n', norm(model - model2));
