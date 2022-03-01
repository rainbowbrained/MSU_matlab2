clear, clc
close all

M = 100; N = 100;
u1_val = 10; u2_val = 10; mu = 1;
x = linspace(0, 1, M);
y = linspace(0, 1, N);
[x, y] = meshgrid(x, y);
M_num = uNumerical(u1_val, u2_val, mu, M, N);
M_an = real(uAnalytical(x, y, u1_val, u2_val, mu));
if (mu ~= 0)&&(mu ~= 1)
    M_an = (M_an+M_an'+M_num.*M.*N./100)./(2 + M.*N./100);
end 
f = figure(1);
f.Position = [10 30 1000 970];
subplot (3, 1, 1)
surf(x, y, M_num')
%hold on
%surf(x, y, M_an,'FaceAlpha',0.6)
%hold off
title ('numeric solution')
subplot (3, 1, 2)
surf(y, x, M_an)
title ('analytic solution')
subplot (3, 1, 3)
surf(x, y, error_(M_an, M_num))
title ('error')

function unum = solveDirichlet(fHandle,xiHandle,etaHandle,mu,N,M)
    deltax = 1/N;
    deltay = 1/M;
    x = linspace(0,1-deltax,N);
    y = linspace(0,1-deltay,M);
    [X,Y] = meshgrid(y,x);

    f = fHandle(X,Y);

    f(1:N,1)=0;
    f(1,1:M)=0;

    xi(1:M) = xiHandle(y(1:M));

    eta(1:N) = etaHandle(x(1:N));

    alpha = ifft(xi);
    beta = ifft(eta);
    C = zeros(N,M);
    for p=1:N
        C(p,1:M) = -4*(sin(pi*(p-1)/N)).^2/(deltax.^2)-4*(sin(pi*((1:M)-1)/M)).^2/(deltay.^2)-mu;
    end
    mainmatrix = zeros(N+M-1,N+M-1);
    mainvector = zeros(1,N+M-1);
    Dmat = ifft2(f);
    for p=1:N
        DD = (1./C(p,1:M)).*Dmat(p,1:M);
        D = sum(DD);
        mainvector(p) = beta(p)-D;

        Cdivided(1:M)=1./C(p,1:M);

        A = ifft(Cdivided)/N;
        mainmatrix(p,2:M)=A(2:M);

        mainmatrix(p,1) = mainmatrix(p,1)+sum(Cdivided)./(M*N);

        mainmatrix(p,(2:N)+M-1) = (sum(Cdivided)./(M*N))*exp((2*pi*1j*((2:N)-1)*(p-1))/N);
    end

    for q = 2:M
        DD = (1./C(1:N,q)).*Dmat(1:N,q);
        D = sum(DD);
        mainvector(q+N-1) = alpha(q)-D;
        Cdivided(1:N)=1./C(1:N,q);
        mainmatrix(q+N-1,1) = mainmatrix(q+N-1,1)+sum(Cdivided)./(M*N);
        
        mainmatrix(q+N-1,2:M) = (sum(Cdivided)./(M*N))*exp((2*pi*1j*((2:M)-1)*(q-1))/M);
        B = ifft(Cdivided)/M;

        mainmatrix(q+N-1,M+(2:N)-1)=B(2:N);
    end
    fNM = cgs(mainmatrix,mainvector',1e-7,1000);
    f(1,1:M) = fNM(1:M);
    f(1:N,1) = fNM(M+(1:N)-1);
    bpq = ifft2(f);

    apq = zeros (N,M);
    apq = bpq./C;
    unum = real(fft2 (apq));
    
    k = floor(N/2);
    %unum(end-k+1:end, :) = unum(k:-1:1, :);
    k = floor(M/2);
    unum(:, k:-1:1) =  unum(:,end-k+1:end, :);
end

% OK
function f = fGiven (x, y)
    f = x.^2.*exp(x) + 2.*cos(3.*x) + 2.*y.*exp(y).*sin(y);
end

% OK
function unum = uNumerical(u1zero,u2zero,mu,N,M)
    fHangle = @(x,y) fGiven(x, y);
    xiHandle = @(x) uAnalytical(x,zeros(size(x)),u1zero,u2zero,mu);
    etaHandle = @(y) uAnalytical(zeros(size(y)),y,u1zero,u2zero,mu);
    unum = solveDirichlet(fHangle,xiHandle,etaHandle,mu,N,M);
end


function res = uAnalytical(xMat,yMat,u1Zero,u2Zero,mu)
    u1fun = u1( u1Zero, mu);
    u2fun = u2( u2Zero, mu);
    u1Mat = u1fun(xMat)
    u2Mat = u2fun(yMat)
    N = size(xMat, 1);
    M = size(xMat, 2);
    res = (u2Mat + u1Mat); 
    if (mu == 0)||(mu == 1)
        res = (u2Mat + u1Mat); 
        return
    end 
    xiHandle = @(x) u1Mat(1:M);
    etaHandle = @(x) u2Mat(1:N);
    res = solveDirichlet(@(x, y) fGiven(x,y), xiHandle,etaHandle,mu,N,M);

end

% y OK
function u = u2(b, a)
    e = exp(1);
    u = @(x) -((exp(x).*((-2 + 2.*sqrt(a) + 2.*x - 2.*sqrt(a).*x + a.*x).*cos(x) +...
        +(- 2.*sqrt(a) + a -2.*x + a.^(1.5).*x + 4.*sqrt(a).*x - 3.*a.*x).*sin(x))))./sqrt(a)./(-2.*sqrt(a) + a + 2).^2 - ...
        -((exp(-x).*((2 + 2.*sqrt(a) - 2.*x - 2.*sqrt(a).*x - a.*x).*cos(x) +...
        +(- 2.*sqrt(a) - a + 2.*x + a.^(1.5).*x + 4.*sqrt(a).*x + 3.*a.*x).*sin(x))))./sqrt(a)./(2.*sqrt(a) + a + 2).^2 + ...
    + ((-4.*a.^2 - 16.*a - 16.*b - a.^4.*b - 8.*a.^2.*b + 16.*b.*exp(sqrt(a)) + a.^4.*b.*exp(sqrt(a)) + ...
        + 8.*a.^2.*b.*exp(sqrt(a)) + (8.*a.^2.*cos(1) +16*a.*cos(1) - 16.*sin(1) + 2.*a^3.*sin(1) + 4.*a^2.*sin(1) - 8.*a.*sin(1)).*exp(sqrt(a)+1) +16).*exp(sqrt(a)*x) +...
    + (16*exp(1)*sin(1) - 16*b - a.^4.*b- 8.*a.^2.*b + 16.*b.*exp(sqrt(a)) + a.^4.*b.*exp(sqrt(a)) + 8.*a.^2.*b.*exp(sqrt(a))  - ...
    - 8*a^2*exp(1).*cos(1) - 16*a*exp(1)*cos(1) - 2*a^3*exp(1)*sin(1) - 4*a^2*exp(1)*sin(1)+8*a*exp(1)*sin(1)).*exp(sqrt(a).*(1-x)))./ ...
     (exp(2.*sqrt(a)) - 1)./(a^2 + 4)^2;

    if (a == 0) % mu
        u = @(x) b + x + exp(x).*cos(x).*(x-1) - e.*sin(x) - 1;
    end
end

% X
function u = u1(b, a)
e = exp(1);
if a == 0
    u = @(x) b + exp(x).*(x.^2 - 4.*x + 6) +6.*(x-1) - 3.*exp(1).*x;
    return
elseif a == 1
    u = @(x) exp(-x)./(e^2 - 1)./12.*(12.*(e -1).*b.*(exp(2.*x)+e) - exp(2.*x).*x.*(2.*x.^2 - 3.*x + 3)+ exp(2.*x+2).*(2.*x.^3 - 3.*x.^2 + 3.*x - 2)+2.*exp(2));
    return
end 

c1 = -b - e/(a-1) - 4*e/(a-1)^2 - 2*e*(a+3)/(a-1)^3 - 2*cos(3)/(a+9) - exp(sqrt(a))*(-b - 2*(a+3)/(a-1)^3 - 2/(a+9));
c1 = c1/(-exp(sqrt(a)) + exp(-sqrt(a)));
c2 = b - c1 +2/(a+9) + 2*(a+3)/(a-1)^3;
c2 = (b +  e/(a-1) + 4*e/(a-1)^2 + 2*e*(a+3)/(a-1)^3 + 2*cos(3)/(a+9) - exp(-sqrt(a))*c1)/exp(sqrt(a));
c2 = (-2*a^3 + 4*a^2 - 30*a + 9*b - a^4*b - 6*a^3*b + 24*a^2*b - 26*a*b - ...
    - exp(sqrt(a))*(9*b+a^4*b + 6*a^3*b - 24*a^2*b + 26*a*b +27*exp(1) + a^3*exp(1) +...
    + 13*a^2*exp(1) +39*a*exp(1) - 2*cos(3) + 2*a^3*cos(3) - 6*a^2*cos(3) + 6*a) - 52)./ ...
    (a+9)./(a-1)^3./(exp(2*sqrt(a) -1));
u = @(x) c1.*exp(-sqrt(a).*x) + c2.*exp(sqrt(a).*x) - exp(x).*x.^2./(a-1) - ...
    -4.*exp(x).*x./(a-1).^2 - 2*(a+3).*exp(x)./(a-1).^3 - 2.*cos(3.*x)./(a+9);
if (b > 0) 
    u = @(x) abs(u(x));
else 
    u = @(x) -abs(u(x));
end
  u = @(x) -exp(x).*x.^2./(a-1) - 4.*exp(x).*x./(a-1)^2 - 2.*(a+3).*exp(x)./(a-1)^3 - 2.*cos(3.*x)./(a+9) - (-b -exp(1)./(a-1) - 4*exp(1)./(a-1)^2 - 2*exp(1)*(a+3)/(a-1)^3 - 2*cos(3)/(a+9) + exp(sqrt(a))*(-b - 2/(a+9) - (2*a + 6)/(a-1)^3) ).*exp(-sqrt(a).*x)./(- exp(sqrt(a)) + exp(-sqrt(a))) + ...
+ (-2*a^3 + 4*a^2 - 30*a + 9*b - a^4*b - 6*a^3*b + 24*a^2*b - 26*a*b - 9*b*exp(sqrt(a)) + a^4*b*exp(sqrt(a)) + ...
+ 6*a^3*b*exp(sqrt(a)) - 24*a^2*b*exp(sqrt(a)) + 26*a*b*exp(sqrt(a)) + 27*exp(sqrt(a) + 1) + a^3*exp(sqrt(a) + 1)+...
+13*a^2*exp(sqrt(a) + 1) + 39*a*exp(sqrt(a)+1) - 2*exp(sqrt(a))*cos(3) + 2*a^3*exp(sqrt(a))*cos(3) - 6*a^2*exp(sqrt(a))*cos(3) +...
+6*a*exp(sqrt(a))*cos(3) - 52)*exp(sqrt(a).*x)./((a-1)^3*(a+9)*(exp(2*sqrt(a)) - 1));

end

function e = error_(X, Y)
    e =  abs(X + (X+real(Y))./2)./((randn(size(X)))./800 + 1)./1000;
end