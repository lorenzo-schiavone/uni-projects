function [x,iter ,resvec ,flag] = mygmres(A,b,tol ,maxit,x0)
% initialization
r0 = b - A*x0;
k = 0;
beta=norm(r0);
resvec = zeros(maxit);
resvec(1)=beta;
V = zeros(size(A,1), maxit);
V(:,1) = r0/beta;
exit_test = tol*norm(b);
flag = 0;
H = zeros(maxit,maxit);

if (resvec(1) < exit_test)
    x=x0;
else
while (resvec(k+1) > exit_test) && (k < maxit)
    k = k+1;
    v_new= A*V(:,k); % v_{k+1}
    % arnoldi method
    h = zeros(k,1);
    for j =1:k
        h(j)=v_new'*V(:,j); %h(j) are the component of the new column of H
        v_new = v_new - h(j)*V(:,j);
    end
    
    H(1:k,k)=h; % H is a square matrix here
    h_new = norm(v_new);
    if h_new < 1e-14 % happy breakdown 
        disp('happy breakdown')
        flag = -1;
        [Q,R]=qr(H(1:k,1:k));
        break
    else
    V(:,k+1)=v_new/h_new;
    H(k+1, k) = h_new;
    end
    [Q,R]=qr(H(1:k+1,1:k));
    resvec(k+1)=abs(beta*Q(1,k+1));
end

rhs = beta* Q(1,1:k)';
y = R(1:k,:) \ (rhs);
x=x0 + V(:,1:k)*y;
end

iter = k;
resvec=resvec(1:(k+1));
end