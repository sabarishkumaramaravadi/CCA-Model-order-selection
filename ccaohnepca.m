clear all
close all
clc
M=1000;
mu=zeros(1,4);
E=[1 0 0.9 0;0 1 0 0.8;0.9 0 1 0;0 0.8 0 1];
S=mvnrnd(mu,E,M);
S=S.';
std1=2;%[2 0; 0 1];
disp('the mixing vector a is')
a1 =zeros(5,2)+randn(5,2);
disp(a1);
disp('the noise vector n is')
n1 = zeros(5,M)+randn(5,M)*0.4;
disp(n1);
s1=zeros(2,M);
disp('the first source ');
for i=1 : 2
    s1(i,:) = S(i,:);
end
disp('the system x is');
x=(a1*s1)+n1;
disp(x);
szx=length(x)
stda2=3;%[1 0; 0 2];
disp('the mixing vector a is')
a2 =randn(5,2)*stda2;
disp(a2);
disp('the noise vector n is')
n2 = zeros(5,M)+randn(5,M)*0.2;
disp(n2);
s2=zeros(2,M);
for i=3 : 4
    j=i-2;
    s2(j,:) = S(i,:);
end
disp(s2);
disp('the system y is');
y=(a2*s2)+n2;

rxx=(1/M)*((x)*(x.'));
ryy=(1/M)*((y)*(y.'));
rxy=(1/M)*((x)*(y.'));
A=(sqrtm(inv(rxx)))*(rxy)*(sqrtm(inv(ryy)));
disp(A);
[F,K,G]=svd(A);
disp(K);