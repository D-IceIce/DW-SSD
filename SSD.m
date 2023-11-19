classdef SSD
    properties
        lambda;
        bata;
        L = 0.05;
        tol = 1e-5;
        rho = 1.01;
        maxIter = 500;

        k
        M
        W
        
        S
        O
        B
        T
        N
        theta
    end
    
    methods
        function obj = process(obj)
            % basis
            obj.S = obj.basis();
            obj.lambda = obj.L/sqrt(size(obj.O,1)*size(obj.O,2));
            obj.bata = 100*obj.L/sqrt(size(obj.O,1)*size(obj.O,2));
            % rpca
            [obj.B, obj.T, obj.N, obj.theta, obj.W] = obj.rpca();
        end
        
        function S = basis(obj)
            kx = 25; ky = 25;
            S{1} = obj.bsplineBasis(size(obj.O,1),kx,5);
            S{2} = obj.bsplineBasis(size(obj.O,2),ky,5);
        end
        
        function S = bsplineBasis(obj, n, k, sd)
            if n == k 
                S = eye(n);
            else
                bd = sd-1;
                knots = [ones(1,bd) linspace(1,n,k) n * ones(1,bd)];
                nKnots = length(knots) - sd;
                kspline = spmak(knots,eye(nKnots));
                S=spval(kspline,1:n)';
            end
        end
        
        function [S, T, N, theta,W] = rpca(obj)
            T = zeros(size(obj.O));
            N = zeros(size(obj.O));
            Y1 = zeros(size(obj.O));
            D{1} = obj.diffMat([size(obj.S{1},2),size(obj.S{2},2)],'circular');
            D{2} = D{1}';
            Y2 = zeros([size(obj.S{1},2),size(obj.S{2},2)]);
            Y3 = zeros([size(obj.S{1},2),size(obj.S{2},2)]);
            Z1 = zeros([size(obj.S{1},2),size(obj.S{2},2)]);
            Z2 = zeros([size(obj.S{1},2),size(obj.S{2},2)]);
            mu = 1e-3;
            for ii = 1:obj.maxIter
                % update W
                W = obj.dynamicW(obj.M,N);
                % update theta
                theta = lyap((obj.S{1}'*obj.S{1})\(D{1}'*D{1}+D{2}'*D{2}),obj.S{2}'*obj.S{2},-(obj.S{1}'*obj.S{1})\(obj.S{1}'*(obj.O-T-N-Y1/mu)*obj.S{2}+D{1}'*(Z1+Y2/mu)+D{2}'*(Z2+Y3/mu)));
                % update Z
                Z1 = obj.solveL2(D{1}*theta+Y2/mu,1/mu);
                Z2 = obj.solveL2(D{2}*theta+Y3/mu,1/mu);
                % update T
                preT = sum(T(:) > 0);
                T = obj.solveL0(obj.O-obj.S{1}*theta*obj.S{2}'-N-Y1/mu,W*obj.lambda/mu);%
                currT = sum(T(:) > 0);
                % update N
                N = obj.solveL21(obj.O-obj.S{1}*theta*obj.S{2}'-T-Y1/mu,obj.bata/mu);
                % check converged
                st = obj.S{1}*theta*obj.S{2}'+T+N-obj.O;
                err = norm(st(:))/norm(obj.O(:));
                if err<obj.tol || (preT>0 && currT>0 && preT == currT)
                    break
                end
                % update mu & Y
                Y1 = Y1 + mu*st;
                Y2 = Y2 + mu*(Z1-D{1}*theta);
                Y3 = Y3 + mu*(Z2-D{2}*theta);
                mu = obj.rho*mu;
            end
            S = obj.S{1}*theta*obj.S{2}';
            % disp(['Iter:', num2str(ii)]);
        end
        
        function Z = solveL1(obj, A,tau)
            Z = max(0,A-tau)+min(0,A+tau);
            Z = max(Z,0);
        end
        
        function Z = solveL21(obj, A, tau)
            Z = A;
            for i = 1:size(A,2)
                tmp = norm(A(:, i), 'fro');
                if tmp > tau
                    coef = (tmp - tau) / tmp;
                    Z(:, i) = coef * A(:, i);
                end
            end
        end

        function W = dynamicW(obj,M,N)
            if size(M,3) >= obj.k
                W0 = obj.subfubc(sum(cat(3,M,N),3),3);
                W1 = mat2gray(W0);
                W = mat2gray(1./(W1+1));
            else
                W = ones(size(N));
            end
        end

        function out = subfubc(obj, img, len)
            m35 = imfilter(double(img), ones(len)/(len*len), 'symmetric');
            op = zeros(3*len, 3*len);
            [X, Y] = meshgrid(len/2:len:len*2.5, len/2:len:len*2.5);
            op(sub2ind(size(op), ceil(X(:)), ceil(Y(:)))) = 1;
            op(ceil(len*1.5), ceil(len*1.5)) = 0;
            m3b = imdilate(m35, op);
            out = (m35 - m3b).^2.*(m35 - m3b>0);
        end
        
        function Z = solveL2(obj, A, tau)
            Z = max(norm(A(:))-tau,0).*A/norm(A(:));
        end
        
        function Z = solveL0(obj, A,tau)
            Z = 1-exp(-A.^2./(2*tau.^2));
        end
        
        function D = diffMat(obj,obs,type)
            D = eye(obs);
            IX = sub2ind([obs obs],2:obs,1:obs-1);
            D(IX) = -1;

            if nargin==3 && strcmp('circular',type)
                D(1,end) = -1;
            elseif nargin==3 && strcmp('none',type) || nargin==2
                D = D(:,1:end-1);
            else
                error('Unrecognized Input');
            end
        end

    end
    
end