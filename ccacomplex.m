clear all
close all
clc;
d=3;%no.of correlated components
r=10;%no.of runs
M=[20 50 100 200];%sample sizes
n=20;%dimension
m=20;%m=n
T=[0.9 0.7 0.5 zeros(1,min([m,n])-d)];
z=1;%for loop to save sungular values in one matrix(sg)
sg=zeros(r,n);
E=[1 0 0 0.9 0 0;0 1 0 0 0 0.7;0 0 1 0 0.5 0;0.9 0 0 1 0 0;0 0 0.5 0 1 0;0 0.7 0 0 0 1];%sigma
for a=1:length(M);
  for l=1:r;%loop of 1000 for averaging singualar values
  mu=zeros(1,length(E));
  S=complex(mvnrnd(mu,E,M(a)),mvnrnd(mu,E,M(a)));
  S=S.';
  %disp('the mixing vector a is')
  a1 =complex(zeros(n,length(E)*0.5)+randn(n,length(E)*0.5)*0.3,(zeros(n,length(E)*0.5)+(randn(n,length(E)*0.5)*0.3)));
  %disp(a1);
  %disp('the noise vector n is')
  n1 = complex(zeros(n,M(a))+randn(n,M(a))*0.1,(zeros(n,M(a))+randn(n,M(a))*0.1));
  %disp(n1);
  s1=zeros(length(E)/2,M(a));
  %disp('the first source ');
    for i=1 : length(E)/2
    s1(i,:) = S(i,:);
    end
  %disp('the system x is');
  x=(a1*s1)+n1;
  %disp(x);
  disp('the mixing vector a is')
  a2 =complex(randn(n,length(E)/2)*0.25,(randn(n,length(E)/2)*0.25));
  %disp(a2);
  %disp('the noise vector n is')
  n2 = complex(zeros(n,M(a))+randn(n,M(a))*0.2,(zeros(n,M(a))+randn(n,M(a))*0.2));
  %disp(n2);
  s2=zeros(length(E)/2,M(a));
    for j=((length(E))/2)+1:length(E) 
    k=j-length(E)/2;
    s2(k,:) = S(j,:);
    end
 %disp(s2);
 %disp('the system y is');
 y=(a2*s2)+n2;
 rxx=(1/M(a))*((x)*(x'));
 ryy=(1/M(a))*((y)*(y'));
 rxy=(1/M(a))*((x)*(y'));
 A=(sqrtm(inv(rxx)))*(rxy)*(sqrtm(inv(ryy)));
 disp(A);
 [F,K,G]=svd(A);
PF=0.001;
s=0;
k=1;
pr=1;
 for s=0:m-1
    ts=chi2inv(1-PF,2*(m-s)*(n-s))
  for i=s+1:m
    pr=pr*(1-K(i,i)^2);
%     disp('pr value');
%     disp(pr)   
  end      
ct=-(2*M(a)-(m+n+1))*log(pr);
pr=1;
disp('ct values')
disp(ct);
if ct<ts
    disp('the number of correlated components');
    ds=s;
    disp(ds);
    break    
end
 end
 %disp(K);
    for t=1:n
    sg(z,t)=K(t,t);
    end
    z=z+1;%to go to next row
     %disp(l);
     %disp(M(a));
end
end
%storing Individual singualr values for M
m1=zeros(r,n);%To store the values of M=20 for 1000 runs
m2=zeros(r,n);
m3=zeros(r,n);%to store the values of M=50 for 1000 runs
m4=zeros(r,n);
for q=1:r
    m1(q,:)=sg(q,:);
end
u=1;
for w=r+1:2*r
    m2(u,:)=sg(w,:);
    u=u+1;
end
h=1;
for g=2*r+1:3*r
    m3(h,:)=sg(g,:);
    h=h+1;
end
   j=1;
for v=3*r+1:4*r
    m4(j,:)=sg(v,:);
    j=j+1;
end
%calculating means
M1=zeros(1,n);
M2=zeros(1,n);
M3=zeros(1,n);
M4=zeros(1,n);%To store the mean of singualr values averaged 1000 times
for i=1:n
M1(1,i)=mean(m1(:,i));
M2(1,i)=mean(m2(:,i));
M3(1,i)=mean(m3(:,i));
M4(1,i)=mean(m4(:,i));
end
N=[M1;M2;M3;M4];
i=1:m;
plot(i,T,'k*')
%disp(i);
hold on
plot(i,M1,'g>');
plot(i,M2,'mx');
plot(i,M3,'bo');
plot(i,M4,'rs');
p=min(n,m);
hold off















    
    




