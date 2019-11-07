function M = liespecialeuclideanfactory()
% Returns a manifold structure to optimize over the special Euclidean group
% 
% function M = specialeuclideanfactory(n)
% function M = specialeuclideanfactory(n, k)
%
% The special Euclidean group (the manifold of rigid transformations):
% This is a product manifold of the rotations group SO(n) and the
% translation group R^n, copied k times. (See note below.)
%
% Points on the manifold are represented as structures X with two fields.
% X.R is a 3D array of size nxnxk such that each slice X.R(:, :, i)
% corresponds to a rotation matrix (orthogonal with determinant 1).
% X.t is a matrix of size nxk such that each column X.t(:, i) corresponds
% to a translation vector.
%
% Tangent vectors are represented as structures with the same fields. Note
% that rotational components of the tangent vectors are represented in the
% Lie algebra, i.e., each slice Xdot.R(:, :, i) is a skew-symmetric matrix.
% Use M.tangent2ambient(X, Xdot) to obtain a representation in the ambient
% space. This is often necessary when defining problem.ehess(X, Xdot).
%
% This is a description of SE(n)^k with the induced metric from the
% embedding space (R^nxn)^k x (R^n)^k, i.e., this manifold is a Riemannian
% submanifold of the embedding Euclidean space with the usual inner
% product.
%
% By default, k = 1.
%
% Note: this is a product geometry: it may not be the "appropriate"
% geometry to give to SE(n) for your application. In particular, this is
% not the Lie geometry of SE(n), because SE(n) is not a direct product of
% SO(n) and R^n: it is only a semidirect product. Following a comment by
% Martijn Zeestraten on the Manopt forum, see this file for more
% information about the Lie geometry:
%   http://ethaneade.com/lie.pdf
%
% See rotationsfactory and euclideanfactory for details.
%
% See also: rotationsfactory euclideanfactory

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Sep. 23, 2014.
% Contributors: 
% Change log:

    n = 3;
    M.nmatrep = 4;
    k = 1;
    
    M.id = eye(M.nmatrep);
    M.inv =@(X) [X(1:3,1:3)',-X(1:3,1:3)'*X(1:3,4); zeros(1,3),1];
    
    M.name = @() sprintf('SE(3) Lie Group');
    
    M.dim = @() 6;
    
    M.inner = @(x, d1, d2) d1(:).'*d2(:);
    
    M.norm = @(x, d) norm(d(:));
    
    M.typicaldist = @() pi*sqrt(n*k);
    
    M.proj = @projliealgebra;
    function projection = projliealgebra(X,H)
        R = X(1:3,1:3);
        RtS = R'*H(1:3,1:3);
        projection  = [0.5*(RtS-RtS'),R'*H(1:3,4);zeros(1,4)];
    end
    
    M.tangent = @(X,H) projliealgebra(M.id,H);
    
    M.tangent2ambient_is_identity = false;
    M.tangent2ambient = @(X, U) multiprod(X, U);
	
	M.egrad2rgrad = M.proj;
	
% 	M.ehess2rhess = @ehess2rhess;
% 	function Rhess = ehess2rhess(X, Egrad, Ehess, H)
%         % Reminder : H contains skew-symmeric matrices. The actual
%         % direction that the point X is moved along is X*H.
% 		Xt = multitransp(X);
% 		XtEgrad = multiprod(Xt, Egrad);
%         symXtEgrad = multisym(XtEgrad);
% 		XtEhess = multiprod(Xt, Ehess);
% 		Rhess = multiskew( XtEhess - multiprod(H, symXtEgrad) );
%     end
    
    % This QR-based retraction is only a first-order approximation
    % of the exponential map, not a second-order one.
    M.retr_qr = @retraction_qr;
    function Y = retraction_qr(X, U, t)
        Y = zeros(M.nmatrep,M.nmatrep);
        if nargin < 3
            Y(1:3,1:3) = qr_unique(X(1:3,1:3) + X(1:3,1:3)*U(1:3,1:3));
            Y(:,4) = X(:,4)+U(:,4);
        else
            Y(1:3,1:3) = qr_unique(X(1:3,1:3) + t*X(1:3,1:3)*U(1:3,1:3));
            Y(:,4) = X(:,4)+t*U(:,4);
        end
    end

%     M.invretr_qr = @inverse_retraction_qr;
%     function S = inverse_retraction_qr(X, Y)
%         
%         % Assume k = 1 in this explanation:
%         % If Y = Retr_X(XS) where S is a skew-symmetric matrix, then
%         %  X(I+S) = YR
%         % for some upper triangular matrix R. Multiply with X' on the left:
%         %  I + S = (X'Y) R    (*)
%         % Since S is skew-symmetric, add the transpose of this equation:
%         %  2I + 0 = (X'Y) R + R' (X'Y)'
%         % We can solve for R by calling solve_for_triu, then solve for S
%         % using equation (*).
%         R = zeros(n, n, k);
%         XtY = multiprod(multitransp(X), Y);
%         H = 2*eye(n);
%         for kk = 1 : k
%             R(:, :, kk) = solve_for_triu(XtY(:, :, kk), H);
%         end
%         % In exact arithmetic, taking the skew-symmetric part has the
%         % effect of subtracting the identity from each slice; in inexact
%         % arithmetic, taking the skew-symmetric part is beneficial to
%         % further enforce tangency.
%         S = multiskew(multiprod(XtY, R));
%         
%     end
    
%     % A second-order retraction is implemented here. To force its use,
%     % after creating the factory M, execute M.retr = M.retr_polar.
%     M.retr_polar = @retraction_polar;
%     function Y = retraction_polar(X, U, t)
%         if nargin == 3
%             tU = t*U;
%         else
%             tU = U;
%         end
%         Y = X + multiprod(X, tU);
%         for kk = 1 : k
%             [Uk, ~, Vk] = svd(Y(:, :, kk));
%             Y(:, :, kk0) = Uk*Vk';
%         end
%     end

%     M.invretr_polar = @inverse_retraction_polar;
%     function S = inverse_retraction_polar(X, Y)
%         
%         % Assume k = 1 in this explanation:
%         % If Y = Retr_X(XS) where S is a skew-symmetric matrix, then
%         %  X(I+S) = YM
%         % for some symmetric matrix M. Multiply with X' on the left:
%         %  I + S = (X'Y) M    (*)
%         % Since S is skew-symmetric, add the transpose of this equation:
%         %  2I + 0 = (X'Y) M + M (X'Y)'
%         % We can solve for M by calling sylvester, then solve for S
%         % using equation (*).
%         MM = zeros(n, n, k);
%         XtY = multiprod(multitransp(X), Y);
%         H = 2*eye(n);
%         for kk = 1 : k
%             MM(:, :, kk) = sylvester_nochecks(XtY(:, :, kk), XtY(:, :, kk)', H);
%         end
%         % In exact arithmetic, taking the skew-symmetric part has the
%         % effect of subtracting the identity from each slice; in inexact
%         % arithmetic, taking the skew-symmetric part is beneficial to
%         % further enforce tangency.
%         S = multiskew(multiprod(XtY, MM));
%         
%     end

    % By default, use QR retraction
    M.retr = M.retr_qr;
%     M.invretr = M.invretr_qr;
    
    % For backward compatibility:
%     M.retr2 = M.retr_polar;
    
    M.exp = @exponential;
    function Y = exponential(X, U, t)
        if nargin == 3
            exptU = t*U;
        else
            exptU = U;
        end
        for kk = 1 : k
            exptU(:, :, kk) = expm(exptU(:, :, kk));
        end
        Y = multiprod(X, exptU);
    end
    
    M.log = @logarithm;
    function U = logarithm(X, Y)
		U = M.inv(X)*Y;
        for kk = 1 : k
            % The result of logm should be real in theory, but it is
            % numerically useful to force it.
            U(:, :, kk) = real(logm(U(:, :, kk)));
        end
        % Ensure the tangent vector is in the Lie algebra.
        U = M.tangent(M.id,U);
    end


    M.hash = @(X) ['z' hashmd5(X(:))];
    
    M.rand =@randompose;
    function randpose = randompose(maxtranslation)
        if nargin<1
            maxtranslation = 1;
        end
        translation = rand(3,1);
        randpose = [randrot(n, k),maxtranslation*translation/norm(translation);zeros(1,3),1];
    end
    
    M.randvec = @randomvec;
    function V = randomvec(X,maxtranslation)
        if nargin<2
            maxtranslation = 1;
        end
        U = randskew(n, k);
        nrmU = sqrt(U(:).'*U(:));
        U = U / nrmU;
        translation = rand(3,1);
        V = [U,maxtranslation*translation/norm(translation);zeros(1,4)];
        V = V / M.norm(M.id,V);
    end
    
    M.lincomb = @matrixlincomb;
    
    M.zerovec = @(x) zeros(M.nmatrep, M.nmatrep, k);
    
    M.transp = @(x1, x2, d) d;
    M.isotransp = M.transp; % the transport is isometric
    
    M.pairmean = @pairmean;
    function Y = pairmean(X1, X2)
        V = M.log(X1, X2);
        Y = M.exp(X1, .5*V);
    end
    
    M.dist = @(x, y) M.norm(x, M.log(x, y));
    
    M.vec = @(x, u_mat) u_mat(:);
    M.mat = @(x, u_vec) reshape(u_vec, [M.nmatrep, M.nmatrep, k]);
    M.vecmatareisometries = @() true;
    
    %% custom field
    M.injectivityradius = 1;
    
    
    M.basisvec =@getbasisvectors;    
    function e = getbasisvectors(i)
        e = zeros(M.nmatrep,M.nmatrep);
        switch i
            case 1
                e(2,3) = -1;
                e(3,2) =  1;
            case 2
                e(1,3) =  1;
                e(3,1) = -1;
            case 3
                e(1,2) = -1;
                e(2,1) =  1;
            case 4
                e(1,4) = 1;
            case 5
                e(2,4) = 1;
            case 6
                e(3,4) = 1;
            otherwise
                error('Cannot access basis vector %d because there only %d of them!',i, M.dim);
        end
    end

    M.gramianmatrix = computegramian();
    function gramian = computegramian()
        gramian = eye(M.dim(),M.dim()); 
        gramian(1:3,1:3) = 2*gramian(1:3,1:3);
    end

    M.tovec =@tovectorcomponents;
    function veccomps = tovectorcomponents(vecmat)
        veccomps = zeros(M.dim(),1);
        for i = 1:M.dim()
            veccomps(i) = M.inner(M.id,M.basisvec(i),vecmat)/ M.inner(M.id,M.basisvec(i),M.basisvec(i));
        end
    end

    M.tomat =@tomatrixrepresentation;
    function matrep = tomatrixrepresentation(vecrep)
        matrep = M.zerovec();
        for i = 1:M.dim()
            matrep = matrep + vecrep(i)*M.basisvec(i);
        end
    end

    M.Ad =@adjointmatrepgroup;
    function Admatr = adjointmatrepgroup(S)
        R = S(1:3,1:3);
        t = S(1:3,4);
        Admatr = [R,zeros(3,3);CrossProdMatrix(t)*R,R];
    end

    M.adjointmap =@(map) M.gramianmatrix\(map')*M.gramianmatrix;

    M.dlog =@derlog;
    function dirder = derlog(x,y)
        
        s = M.log(x,y);
        omeu = M.tovec(s);
        ome = omeu(1:3);
        u = omeu(4:6);

        theta2 = trace(s(1:3,1:3)'*s(1:3,1:3))/2;
        theta = sqrt(theta2);
        if theta2>1e-4
            aTheta = sin(theta) / theta;
            bTheta = (1-cos(theta)) / theta2;
            cTheta = (1-aTheta) / theta2;
            eTheta = (bTheta-0.5*aTheta)/(1-cos(theta));
        else
            aTheta = 1 - theta.*theta/factorial(3) + theta.^4/(factorial(5));
            bTheta = 0.5 - theta.*theta/24 + theta.^4/(factorial(6));
            cTheta = 1/(factorial(3)) - theta.*theta/(factorial(5)) + theta.^4/(factorial(7));
            eTheta = (bTheta-2*cTheta)/(2*aTheta);
        end
        
        dlogrot = (eye(3) -0.5*s(1:3,1:3) + eTheta * s(1:3,1:3) * s(1:3,1:3) );
        W = (cTheta-bTheta)*eye(3) + (aTheta-2*bTheta)/(theta2)*s(1:3,1:3) + (bTheta-3*cTheta)/(theta2)*(ome*ome');
        B = bTheta*CrossProdMatrix(u)+cTheta*(ome*(u')+u*(ome'))+(ome'*u)*W;
        dirder = [dlogrot,zeros(3,3);-dlogrot*B*dlogrot,dlogrot]*M.Ad(M.inv(x)*y);
    end

end
