clc;clear;
%% data name

analytic_para='analytic_para.txt';
dataname_para='dataname_para.txt';
dataname_best='28_best1.dat';
dataname_error='28_error1.dat';

%% read analytic_para.txt
fid = fopen(analytic_para);
tline = fgetl(fid);
RR=str2num(tline);
tline = fgetl(fid);
xN=str2num(tline);
x=zeros(1,xN);
i=1;
while(i<=xN)
    tline = fgetl(fid);
    x(i)=str2double(tline);
    i=i+1;
end
fclose(fid);

%% read dataname_para.txt
fid = fopen(dataname_para);
tline = fgetl(fid);
d=str2num(tline); 
Ns=zeros(1,d);
i=1;
while(i<=d)
    tline = fgetl(fid);
    Ns(i)=str2double(tline);
    i=i+1;
end
Ks=zeros(d,1);
i=1;
while(i<=d)
    tline = fgetl(fid);
    Ks(i)=str2double(tline);
    i=i+1;
end
tline = fgetl(fid);
A=str2double(tline);
fclose(fid);

%% message
disp('*********************************** reading parameters ')
disp('system information')
disp(d)
disp(Ns ) %row vector
disp(Ks)  %
disp(A)
disp(RR)
disp(xN) %column vector


%% read data

M_best = dlmread(dataname_best);
M_error = dlmread(dataname_error);

m_best = reshape(M_best,Ns);
m_error = reshape(M_error,Ns);

m_sig=m_error./m_best;
disp('largest relative error is')
error_limit=max(max(abs(m_sig)));
disp(error_limit) ;

% now m_best and m_sig are high-D tensor to be analytic continued.
% d=2 case, imagesc(m_best)

% the goal is to generate similar stuff, with N0 replaced
% (notice that Matlab index starts with 1)
% Nss is the output version of Ns
Nss=Ns;
Nss(1)=xN;

output_best=zeros(Nss);
output_error=zeros(Nss);


%% perform the AC

z=(0:(Ns(1)/2-1))*pi/Ns(1);


%d=2
for i=1:Ns(2)
    disp(i);
    u_real=m_best(1:(Ns(1)/2),i);
    err_u_real=m_sig(1:(Ns(1)/2),i)*0.5;
    u_imag=zeros(1,Ns(1)/2);
    err_u_imag=u_imag;
    [output_best(:,i),output_error(:,i)]=PadeRegression(z,u_real,u_imag,err_u_real,err_u_imag,x,RR);
end

%[rho_best,rho_error]=PadeRegression(z,u_real,u_imag,err_u_real,err_u_imag,x,RR);


%% check output relative error
figure(1);
p=1;
imagesc( flipud(    -( output_best )      ) );
colorbar();
%imagesc( flipud(   ( (  output_best.^p ) + (  output_error.^p )  ).^(1/p)   ) );
%colorbar();

figure(3);
imagesc( flipud(   log( abs(  output_error./output_best )  )   ) );
colorbar();
 
ee=  reshape(abs( output_error./ output_best ), 200*16,1);
figure(2);
hist(log10(ee));
