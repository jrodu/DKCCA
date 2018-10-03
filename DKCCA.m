function [ dcell, ca, reg, dirsa, kernelMats, corcells ] = DKCCA(A, B, fit_per, window_length, np, numcomps, reggrid, gsd)

    % A:                 an array of size q_A x T x N where q_A is the number of units, T
    %           the number of time points, and N the number of trials
    % B:                 an array of size q_B x T x N
    % fit_per:           number of time steps fit for each location of the sliding
    %           window.  The sliding window slides in steps of length fit_per
    % window_length:     the total length of the sliding window.  If fit_per
    %           is odd, window_length needs to be odd.  If fit_per is even,
    %           window_length needs to be even.  An error will result if this is not
    %           the case.
    % np:                number of permutations to run, 0 returns just the raw correlogram. 
    % numcomps:          the number of canonical components to extract.  The default
    %           is 1.
    % reggrid:           the regularization grid to search over for kernel CCA.  Defaults
    %           to .7 to 1 in steps of .01
    % gsd:               number of latent factors for gram schmidt
    %           decomposition for kernel CCA. Defaults to 5
    
    
    
    %Note:    currently this algorithm is impmlemented with the linear
    %           kernel.  


%calculate kernel matrices
kernelMats=cell(2, size(A, 2));
A=A-repmat(mean(A, 3), 1, 1, size(A, 3));
B=B-repmat(mean(B, 3), 1, 1, size(B, 3));

for ind=1:size(A, 2)
    kernelMats{1,ind}=linearkernel(squeeze(A(:, ind,:)));
    kernelMats{2,ind}=linearkernel(squeeze(B(:, ind,:)));

end

[ca, reg, dirsa, corcells]=getkernelcormat2(kernelMats, fit_per, window_length, numcomps, np, reggrid, gsd);

dcell=get_statistics(ca);

end


