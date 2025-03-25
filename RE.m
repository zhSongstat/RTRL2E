function [Xk1,A,B] = RE(T,Omega,Omegac,lambda,iR,K,tol)
X = T;
[n1,n2,n3] = size(X);
X(Omegac) = 0;
flag = 0;
while flag == 0
    for k = 1:K
        Xk = X;
        X_hat = fft(Xk,[],3);
        for i = 1:n3
            [U0,~,V0] = svd(X_hat(:,:,i),'econ');
            U_hati = U0(:,1:iR);
            V_hati = V0(:,1:iR);
            t = diag(U_hati.'*X_hat(:,:,i)*V_hati);
            S_hati = my_prox_l1(diag(t),lambda);
            X_hat2(:,:,i) = U_hati*S_hati*V_hati.';
        end
        X = ifft(X_hat2,[],3);X(Omega) = T(Omega);
        Xk1 = real(X);
        if norm(Xk1(:)-Xk(:))/norm(Xk1(:)) < tol
            flag = 1;
            
        elseif k==K
            flag = 1;
        end
    end
    
end
Ss = iR;
Xest = fft(Xk1,[],3);
for i = 1 : n3
    [U,S,V] = svds(Xest(:,:,i),Ss);
    A(:,:,i) = U*S;
    B(:,:,i) = V.';
end
A = ifft(A,[],3);
B = ifft(B,[],3);
end