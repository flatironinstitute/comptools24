function i = diagind(A)
% function i = diagind(A)
%
% Return all diagonal indices i of a square matrix A.
%
% Example: A = rand(3); A(diagind(A)) = 1;   % sets diagonal to 1

% Barnett 2/6/08

N = size(A,1);
if size(A,2)~=N
  disp('input must be square!');
end
i = sub2ind(size(A), 1:N, 1:N);
