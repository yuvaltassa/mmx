function mmx_nothreads()
%mmx_nothreads - Multithreaded matrix operations on N-D matrices but without threading
%    MMX treats an N-D matrix of double precision values as a set of pages 
%    of 2D matrices, and performs various matrix operations on those pages.  
%    MMX uses multithreading over the higher dimensions to achieve good
%    performance. Full singleton expansion is available. 
% 
%    C = mmx_nothreads(A, B) is equivalent to the matlab loop
%    for i=1:N,
%        C(:,:,i) = A(:,:,i) * B(:,:,i);
%    end
%
%
%    Singleton expansion is enabled on all dimensions so for example if
%    A = randn(5,4,3,10,1);
%    B = randn(4,6,3,1 ,6);
%    C = zeros(5,6,3,10,6);
%    then C = mmx_nothreads(A,B) is equivalent to 
%    for i = 1:3
%       for j = 1:10
%          for k = 1:6
%             C(:,:,i,j,k) = A(:,:,i,j,1) * B(:,:,i,1,k);
%          end
%       end
%    end
error(sprintf('MEX file not found.\nTry ''build_mmx''.\nType ''help mmx'' for details.'));