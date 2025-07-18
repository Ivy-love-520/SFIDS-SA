
function d = sampleEntropy(seq, wlen, r, shift)
%
% Sample Entropy (matlab-version)
%
%   SampEn = sampleEntropy(INPUT, M, R, TAU)
%
%   
%   Arguments:
%       INPUT       Nx1         Input sequence.
%       M           Int         Window-length (or "dimension").一般取2
%                              
%       R           Double      Tolerance for "similarity". 一般取0.2*std
%       TAU         Int         Spacing of valid samples (for subsampling).
%                               A value of 1 corresponds to no subsampling,
%                               2 takes every other value, etc.




if shift > 1,
    seq = downsample(seq,shift);
end

% allocate space for extracted windows;
D = zeros(length(seq)-(wlen+1), wlen+1);

% extract windows with length wlen+1
for pos=1:length(seq)-wlen-1,
    D(pos,:) = seq(pos:pos+wlen);
end

% initialise
A = 0;
B = 0;

% calculate number of windows with pairwise distance of less than r, for
% two cases:
%   1) B = with windows = 1..wlen 
%   2) A = with windows = 1..wlen+1
for i=1:size(D,1),
    % Chebyshev distance is max(abs(d_ik-d_jk))
    % D(i,i) is 0, but we should not count that.
    % Also D(i,j) is symmetrical (d(i,j)=d(j,i)), therefore we just need to
    % look at D(i+1:end). Effectively we only calculate "half" of the
    % distance matrix. Due to symmetry we can ignore the rest.
    DD = bsxfun(@minus, D(i+1:end,:), D(i,:)); % subtract current window from all future windows.
    DD = abs(DD); % DD now cheb. distance
    
    v1 = max(DD(:,1:end-1),[],2); % maximum along 2nd dim (case 1)
    v2 = max(v1, DD(:,end));      % add last column (case 2)
    
    B = B + sum(v1 < r);
    A = A + sum(v2 < r);
end

% A contains half the matches,
% B contains half the matches. For estimating A/B this doesn't matter
% really.
d = -log(A/B);

end



