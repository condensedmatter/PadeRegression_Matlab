%% Function Name: NA_0025_Example_Header
% Author: Jian Wang
% University of California, Los Angeles
% Date: Nov 21, 2018

function [rho_best,rho_error]=PadeRegression(z,u_real,u_imag,err_u_real,err_u_imag,x,RR)
%% average over RR Pade-Regressions
% best value and standard error is produced
rho_ensemble=zeros(RR,size(x,2)*size(x,1));
parfor R=1:RR
    rho_ensemble(R,:)=PadeRegressionOne(z,u_real,u_imag,err_u_real,err_u_imag,x,R);     
end
rho_best=mean(rho_ensemble);
rho_error=std(rho_ensemble);
end

function Rho=PadeRegressionOne(z_matsubara,u_real,u_imag,err_u_real,err_u_imag,x,R_seed)
rng(R_seed);
N=size(z_matsubara,1)*size(z_matsubara,2);
z=zeros(N,1);
u=z;
%% notice:
%(1)sigReal(i) and sigImag(i) are all relative error
%(2)GPadeImag(i) are zero in some problem, for example spin problem
%(3)sigReal(i) and sigImag(i) should be correlation matrix in more serious
for i=1:N
    z(i)=z_matsubara(i)*1i;
    %u(i)=(u_real(i) + u_imag(i)*1i)*(1.0+ (normrnd(0,err_u_real(i))+normrnd(0,err_u_imag(i))*1i));
    u(i)=u_real(i) * (1.0+normrnd(0,err_u_real(i))) + 1i*u_imag(i)*(1.0+normrnd(0,err_u_imag(i)));
    %(4)there might be some variations in coding the above line
end
%(5)A,B are polynomials, in current treatment, there [order] satisfy this:
% [A]=[B]=[dimension of input z]/2
%further work need to generalize this 
[A,B]=VandermondePoly(z,u);
Rho=x;
for i=1:size(x,1)*size(x,2)
    Rho(i)=imag(VandermondePolyEvaluate(x(i)  ,A,B))/(-pi);
end
end

% given z and u, find the rational polynomial coefficient
% z and u should be a column vector
% for example  ; delete function defination line, run as a script
% z=[1:50]'*1i;
% u=(z.^2-z-10)./(z.^3+z+3);

function [A,B]=VandermondePoly(z,u)
%% (1) set the dimension N, na, nb
% the dimension of input array z or u
N=size(z,1)*size(z,2); 
% na, nb is the order(degree) of nominator/denominator polynomial
if rem(N,2)==1   % N is the numbers of input data points
    na=(N-1)/2;  % na is the highest power of nominator polynomial
    nb=na;       % nb is the highest power of denominator polynomial
else
    na=N/2;
    nb=na-1;
end
% na=6;
% nb=8;
% this this perticular method, na+nb+1=N
% notice that, the order of unknow coefficients:
% nominator 0,1,2,...,na
% denominator 1,2,3,...,nb
M=zeros(N,N);
% for robustness
if size(z,1)<size(z,2)
    z=z.';   % z=z' BUG!!!  non-congugate transpose
end
if size(u,1)<size(u,2)
    u=u.';
end
%% (2) Construct the Pade-Vandermonde-Matrix
for i=1:nb
    M(:,i)=z.^i;
end
for i=1:N
    M(i,:)=-M(i,:)*u(i);
end
for i=(nb+1):N
    M(:,i)=z.^(i-nb-1);
end
%% (3) perform the regression, hence the name "Rational function (Pade) Regression"
x=M\u;
%% (4) get the coefficient from x
B=[1; x(1:nb)];
A=x((nb+1):N);
end

function u=VandermondePolyEvaluate(z,A,B)
%% this convert to standard polynomial array forms
if size(A,1)>size(A,2)
    matlabA=fliplr(A.');
else
    matlabA=fliplr(A);
end

if size(B,1)>size(B,2)
    matlabB=fliplr(B.');
else
    matlabB=fliplr(B);
end

%% evaluate the Rational function (Pade)
u=polyval(matlabA,z)/polyval(matlabB,z);

end
