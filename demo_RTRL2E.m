function Lk1 = demo_RTRL2E(X,L,U,V,W,etamax,beta3,mu3,rho,lambda,Omega,R,itmax,eps)
N = sum(W(:));
[~,~,n3] = size(X);
alpha = 0.01; 
beta = 0.01;
flag = 0;
%% initial 
eta = log(0.01);
I = eye(R);
while flag == 0
    for k = 1:itmax
        Lk = L;

        L_hat = fft(L,[],3);
        U_hat = fft(U,[],3);
        V_hat = fft(V,[],3);
        for i = 1:n3
            U_hat(:,:,i) = (mu3.*L_hat(:,:,i))*(V_hat(:,:,i).')*pinv(beta3.*I+mu3.*V_hat(:,:,i)*V_hat(:,:,i).');
            V_hat(:,:,i) = pinv(mu3.*U_hat(:,:,i).'*U_hat(:,:,i)+beta3.*I)*U_hat(:,:,i).'*(mu3.*L_hat(:,:,i));
        end
        L = ifft(L_hat,[],3);
        U = ifft(U_hat,[],3);
        V = ifft(V_hat,[],3);
        
        P = -exp(3*eta).*sqrt(2/pi).*W.*exp(-0.5*exp(2*eta).*((X-L).^2)).*(X-L);
        UV = tprod(U,V);UV(Omega) = X(Omega);
        D_L = P + mu3.*(L - UV)+rho.*L;
        L = L - alpha.*D_L;
        
        E0 = W.*exp(-0.5.*exp(2*eta).*((X-L).^2));
        E1 = W.*exp(-0.5.*exp(2*eta).*((X-L).^2)).*(X-L).^2;
        D_eta = 0.5*N*exp(eta)./sqrt(pi) - exp(eta)*sqrt(2/pi)*sum(E0(:)) + exp(3*eta)*sqrt(2/pi)*sum(E1(:)) + lambda;
        eta = eta - alpha.*D_eta;
        D_lambda = eta-etamax;
        lambda = max(lambda+beta.*D_lambda,0);
        Lk1 = L;
        if norm(Lk1(:)-Lk(:))/norm(Lk(:)) < eps
            flag = 1;
        elseif k == itmax
            flag = 1;
        end
    end
end
end