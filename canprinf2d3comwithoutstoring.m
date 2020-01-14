clear all
close all
clc;
z=1;
C =['r^','o','x','s','*'];
d=3;%no.of correlated components
r=100;%no.of runs
M=100;%sample sizes
n=10;%dimension
m=10;%m=n
T=[0.9 0.8 0.7 zeros(1,min([m,n])-d)];
R=[3 4 5 7 9 10];
E=[1 zeros(1,8) 0.9;0 1 zeros(1,6) 0.8 0;0 0 1 zeros(1,4) 0.7 zeros(1,2) ;zeros(1,3) 2 zeros(1,5) 0;zeros(1,4) 2 zeros(1,4) 0;zeros(1,5) 2 zeros(1,4);zeros(1,6) 2 zeros(1,3);0 0 0.7 zeros(1,4) 1 0 0;0 0.8 zeros(1,6) 1 0;0.9 zeros(1,8) 1];
   
   for b=1:length(R)
      for l=1:r;%loop of 1000 for averaging singualar values
  mu=zeros(1,length(E));
  S=complex(mvnrnd(mu,E,M),mvnrnd(mu,E,M)); %mixing signal
  S=S.';
  a1 =complex(zeros(n,length(E)*0.5)+randn(n,length(E)*0.5)*0.3,(zeros(n,length(E)*0.5)+(randn(n,length(E)*0.5)*0.3)));
  n1 = complex(zeros(n,M)+randn(n,M)*0.1,(zeros(n,M)+randn(n,M)*0.1)); %noise
  s1=zeros(length(E)/2,M);
    for i=1 : length(E)/2%disp('the first source ');
    s1(i,:) = S(i,:);
    end
    
  %Model
  x=(a1*s1)+n1;
  
  [U,W,V]=svd(x,'econ');
  rank = R(b);
  xrx=(U(:,1:rank)')*x;
  disp(xrx);
  
  a2 =complex(randn(n,length(E)/2)*0.25,(randn(n,length(E)/2)*0.25));
  n2 = complex(zeros(n,M)+randn(n,M)*0.2,(zeros(n,M)+randn(n,M)*0.2));
  %disp(n2);
  s2=zeros(length(E)/2,M);
    for j=((length(E))/2)+1:length(E) 
    k=j-length(E)/2;
    s2(k,:) = S(j,:);
    end
    
 y=(a2*s2)+n2;
 
 [U,W,V]=svd(y,'econ');
 
 yry=U(:,1:rank)'*y;
 rxx=(1/M)*((xrx)*(xrx'));
 ryy=(1/M)*((yry)*(yry'));
 rxy=(1/M)*((xrx)*(yry'));
 A=(sqrtm(inv(rxx)))*(rxy)*(sqrtm(inv(ryy)));
 disp(A);
[F,K,G]=svd(A);

for t=1:rank
   ss(l,t)= K(t,t);
end


      end
      sm=zeros(1,rank);
      
sm=mean(ss);
PF=0.001;
s=0;
k=1;
pr=1;
L=1;
for s=0:rank-1
    ts=chi2inv(1-PF,2*(rank-s)*(rank-s)); %chisquare distribute-Detector
  for i=s+1:rank
    pr=pr*(1-sm(1,i)^2);  
  end 
  for j=1:s+1
     L= L*K(i,i)^2;
  end
  
  
ct=-(2*M-(2*(rank)+1)+L)*log(pr); %ITC-Detector
h(b,l)=ct;
pr=1;
L=1;
disp('ct values')
disp(ct);
if ct<ts
    disp('the number of correlated components');
    ds=s;
    disp(ds);
    break  
elseif s==rank-1
    ds=s;
end
    
end
      j=1:rank;
      
      
 plot(j,sm,C(z));
hold on
z=z+1;
   end
   hisogram(h);
   legend('show')
   xlabel('index-i');
   ylabel('sample-Ki');
i=1:length(T);
plot(i,T,'k*');
hold off





