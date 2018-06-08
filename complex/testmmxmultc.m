classdef testmmxmultc < matlab.unittest.TestCase
    properties (TestParameter)
        % number types
        A_type = {'REAL','COMPLEX'};
        B_type = {'REAL','COMPLEX'};
        % transpose
        A_trans = {'N'};
        B_trans = {'N'};
        % matrix size
        F = struct('small', 1,'medium', 2, 'large', 10);
        dim_inner = struct('small', 1,'medium', 2, 'large', 4);
        dim_outer_1 = struct('small', 1,'medium', 2, 'large', 4);
        dim_outer_2 = struct('small', 1,'medium', 2, 'large', 4);
        function_to_test = {@mmxmultc,@mmxmultc_nothreads};
    end
    methods (Test)
        function test_multiplication(testCase,function_to_test,A_type,B_type,A_trans,B_trans,F,dim_inner,dim_outer_1,dim_outer_2)
%             fprintf('%s %dx%d times %dx%x %d freqs\n',[A_type ' ' B_type ' ' A_trans B_trans],dim_outer_1,dim_inner,dim_inner,dim_outer_2,F );
            % generate the sizes of the matrices
            if strcmpi(A_trans,'T')
                if strcmpi(B_trans,'T')
                    % AtBt
                    Asize = [dim_inner dim_outer_1 F];
                    Bsize = [dim_outer_2 dim_inner F];
                else
                    % AtB
                    Asize = [dim_inner dim_outer_1 F];
                    Bsize = [dim_inner dim_outer_2 F];
                end
            else
                if strcmpi(B_trans,'T')
                    %ABt
                    Asize = [dim_outer_1 dim_inner F];
                    Bsize = [dim_outer_2 dim_inner F];
                else
                    % AB
                    Asize = [dim_outer_1 dim_inner F];
                    Bsize = [dim_inner dim_outer_2 F];
                end
            end
            % generate the matrices
            A = rand(Asize);
            if strcmpi(A_type,'COMPLEX')
                A = A+1i.*rand(Asize);
            end
            B = rand(Bsize);
            if strcmpi(B_type,'COMPLEX')
                B = B+1i.*rand(Bsize);
            end
            % apply the function
            try
                C = function_to_test(A,B,[A_trans B_trans]);
            catch err
                keyboard
            end
            % compute the same result with a simple for loop
            Ctest = zeros(dim_outer_1,dim_outer_2,F);
            for ff=1:F
                Atemp = A(:,:,ff);
                Btemp = B(:,:,ff);
                if A_trans=='T'
                    Atemp = Atemp.';
                end
                if B_trans=='T'
                    Btemp = Btemp.';
                end
                Ctest(:,:,ff) = Atemp*Btemp;
            end
            % check whether the two are equal
            testCase.verifyEqual(C,Ctest,'AbsTol',1e-15);
        end
    end
end