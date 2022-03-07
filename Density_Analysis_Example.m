%Author by Oralia Nolasco
%base on: 
% n x n matrix
% s -> standard deviation
% dx -> histogram bin

function [f, x, N] = Density_Analysis(n, s, dx)
  v = [];

  for i = 1:100 % 100 k-th moments

      H = s*(randn(n)+1i*randn(n)); %random Gaussian Matrix
      for j = 1:n
          H(i, j) = sqrt(2)*H(j, j);
      end
  W = (H+H')/2; %Symetric Hermitian matrix
  v = [v; eig(W)];
  lmax = max(eig(W));
  end
  v = v/sqrt(n); % normalized eigenvalues
  [f, x] = hist(v, -2*s:dx:2*s);
   N = norm(W)/(s*sqrt(n));
end