classdef converter
    %CONVERTER Summary of this class goes here
    %   Detailed explanation goes here
    
    methods (Static) 
        function [vecFromVech] = getVecFromVech(n)
            % -------------------------------------------------------------------------
            % computes the auxiliary matrix vecFromVech such that: vec(Sigma)=vecFromVech vech(Sigma)
            %
            % Inputs:
            % - n: number of variables
            % Outputs:
            % - vecFromVech: auxiliary matrix
            %
            % This version: February 24, 2015
            % -------------------------------------------------------------------------
            
            Aux = eye(n);
            last = zeros(n^2,1);
            last(n^2,1)=1;
            V(1).V = last;
            for j=2:n
                A = eye(j);
                if j<n
                    B = [zeros( n*(n-j)+n-j, j);eye(j)];
                    for m=2:j
                        B = [B; kron(Aux(:,n-j+1),A(m,:))];
                    end
                    clear A;
                    V(j).V = [B,V(j-1).V];
                    clear B
                else
                    B = eye(j);
                    for m=2:j
                        B = [B;kron(Aux(:,n-j+1),A(m,:))];
                    end
                    clear A;
                    V(j).V = [B,V(j-1).V];
                    clear B
                end
            end
            vecFromVech = V(n).V;
            
        end
        function [vechFromVec] = getVechFromVec(n)
            Identity = eye(n);
            vechFromVec = kron(Identity(1,:),Identity);
            for i=2:n
                vechFromVec = [vechFromVec; kron(Identity(i,:),Identity(i:end,:))];
            end

        end
    end
    
end

