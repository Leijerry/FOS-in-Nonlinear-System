clc
clear all
%System=1;
%System=2;
System=3;
 
if System==1
    x=0.1*randn(1,3000);
    plot(x);
    y=Noise_less_system_1(x);
    figure;
    plot(y);
    k=[2,5,7,6,4,3,5,3,4,5];
    l=[1,4,3,7,1,4,6,2,1,2];
    Order=[2,2,2,2,2,2,2,2,2,2];
end;
if System==2
    x=rand(1,3000);
    plot(x);
    y=Noise_less_system_2(x);
    figure;
    plot(y);
    k=[2,2,3,3,4,4];
    l=[1,2,3,4,4,5];
    Order=[2,2,2,2,2,2];
end;
if System==3
    x=0.1*randn(1,3000);
    y=Noise_less_system_3(x);
    k=[2,2,3,3,4,4];
    l=[1,2,3,4,4,5];
    Order=[2,2,2,2,2,2];
end;
percent=[0 25 50 100];
mse_absolute_3rd_set=zeros(size(k,2),4);
mse_relative_3rd_set=zeros(size(k,2),4);
for i=1:size(percent,2)
    for model=1:size(k,2)
        disp(['For model', num2str(model),'(when P is ',num2str(percent(i)),'%)']);
        tic
        a(model)=FOS(x,y,percent(i),k(model),l(model),Order(model));
        toc
        mse_absolute_1st_set(model,i)=a(model).mse_absolute_1st_set;
        mse_relative_1st_set(model,i)=a(model).mse_relative_1st_set;
        mse_absolute_2nd_set(model,i)=a(model).mse_absolute_2nd_set;
        mse_relative_2nd_set(model,i)=a(model).mse_relative_2nd_set;
    end
    [q(i),p(i)]=min(mse_absolute_2nd_set(:,i));
    b(i)=a(p(i));
    modelout(a(p(i)).Selected_candidates,a(p(i)).a,a(p(i)).K,a(p(i)).L,percent(i),a(p(i)).Order);
    v3(i,1:1000)=Modloutput(a(p(i)).a,a(p(i)).Selected_candidates,x(2001:3000),a(p(i)).K,a(p(i)).L,a(p(i)).N,a(p(i)).Order);
    mse_absolute_3rd_set(p(i),i)=a(p(i)).mse_absolute_3rd_set;
    mse_relative_3rd_set(p(i),i)=a(p(i)).mse_relative_3rd_set;
    if i==1
        figure;
        P_lot=plot(2001:3000,y(1,2001:3000),'--',2001:3000,v3(i,1:1000));
        set(P_lot(1),'Color','blue','linewidth',2);
        set(P_lot(2),'Color','red','linewidth',1);
        figure;
        P_lot=plot(2001:2100,y(1,2001:2100),'--',2001:2100,v3(i,001:100));
        set(P_lot(1),'Color','blue','linewidth',2);
        set(P_lot(2),'Color','red','linewidth',1);
        figure;
        P_lot=plot(2450:2550,y(1,2450:2550),'--',2450:2550,v3(i,450:550));
        set(P_lot(1),'Color','blue','linewidth',2);
        set(P_lot(2),'Color','red','linewidth',1);
        figure;
        P_lot=plot(2900:3000,y(1,2900:3000),'--',2900:3000,v3(i,900:1000));
        set(P_lot(1),'Color','blue','linewidth',2);
        set(P_lot(2),'Color','red','linewidth',1);
    else
        figure;
        P_lot=plot(2001:3000,y(1,2001:3000),2001:3000,v3(i,1:1000));
        set(P_lot(1),'Color','blue','linewidth',1);
        set(P_lot(2),'Color','red','linewidth',1);
        P_lot=plot(2001:2100,y(1,2001:2100),2001:2100,v3(i,001:100));
        set(P_lot(1),'Color','blue','linewidth',1);
        set(P_lot(2),'Color','red','linewidth',1);
        figure;
        P_lot=plot(2450:2550,y(1,2450:2550),2450:2550,v3(i,450:550));
        set(P_lot(1),'Color','blue','linewidth',1);
        set(P_lot(2),'Color','red','linewidth',1);
        figure;
        P_lot=plot(2900:3000,y(1,2900:3000),2900:3000,v3(i,900:1000));
        set(P_lot(1),'Color','blue','linewidth',1);
        set(P_lot(2),'Color','red','linewidth',1);
    end
end;
result(1:size(k,2),1)=[1:size(k,2)]';
result(1:size(k,2),2)=k';
result(1:size(k,2),3)=l';
result(1:size(k,2),4)=Order';
 
rr=4;
for ll=1:1:size(percent,2)
    rr=rr+1;
    result(:,rr)=100*percent(ll)/(percent(ll)+100);
    rr=rr+1;
    result(:,rr)=mse_absolute_1st_set(:,ll);
    rr=rr+1;
    result(:,rr)=mse_relative_1st_set(:,ll);
    rr=rr+1;
    result(:,rr)=mse_absolute_2nd_set(:,ll);
    rr=rr+1;
    result(:,rr)=mse_relative_2nd_set(:,ll);
    rr=rr+1;
    result(:,rr)=mse_absolute_3rd_set(:,ll);
    rr=rr+1;
    result(:,rr)=mse_relative_3rd_set(:,ll);
end;
result1=result';
function y = Noise_less_system_1(x)
y = zeros ( 1,3000);
for i = 5:3000
    y(n) = 1+0.23*x(n) + 0.16*x(n-1) + 0.13*x(n-2)*x(n-3) + 0.08*(x(n-1)).^2 - 0.26*(y(n)).^2 - 0.08*y(n)*y(n-1) + 0.23*y(n-2)*y(n-3); 
end
function y = Noise_less_system_2(x)
y = zeros (1,3000);
for i = 5:3000
    y(i) = 0.098 + 0.413*x(i-2) + 0.597*y(i-1) - 0.78*x(i-3)*x(i-1) + 1.12*x(i-1)*y(i-2) - 0.7*y(i-1)*y(i-2); 
end

function y = Noise_less_system_3(x)
y = zeros ( 1,3000);
for i = 5:3000
    y(i) =0.1 + 0.4*x(i) - 0.3*y(i)  - 0.9*x(i-1)*y(i-2) + 0.76*y(i-1)*y(i-2);
end

classdef FOS
        properties ( SetAccess = public )
            Selected_candidates = 0;
            z=0;
            K=0;
            L=0;
            N=0;
            a=0;
            Order=0;
            mse_absolute_1st_set = 0;
            mse_relative_1st_set = 0;
            mse_absolute_2nd_set = 0;
            mse_relative_2nd_set = 0;
            mse_absolute_3rd_set = 0;
            mse_relative_3rd_set = 0;
    end
    
    methods
        function obj = FOS(x,y,Noise_percentage,K,L,Order)
            
 [obj.z,obj.Order,obj.Selected_candidates,obj.K,obj.L,obj.N,obj.a,obj.mse_absolute_1st_set,obj.mse_relative_1st_set,obj.mse_absolute_2nd_set,obj.mse_relative_2nd_set,obj.mse_absolute_3rd_set,obj.mse_relative_3rd_set]= fast_orthogonal_search(x,y,Noise_percentage,K,L,Order);
    end
    end
end


function [z,Order,Selected_candidates,K,L,N,a,mse_absolute_1st_set,mse_relative_1st_set,mse_absolute_2nd_set,mse_relative_2nd_set,mse_absolute_3rd_set,mse_relative_3rd_set]=fast_orthogonal_search(x,y,Noise_percentage,K,L,Order)
z=add_noise(y,Noise_percentage);
N= max ([K;L]);
M=factorial(K+L+1+Order)/(factorial(K+L+1)*factorial(Order));
Candidates= transpose(1:M);
E(1) =1;
G(1)= sum(z(N+1:1000))/(1000-N);
Selected_candidates(1) = Candidates(1,1);
A=0;
M_S_E = sum(z(N+1:1000).^2)/(1000-N)-G(1)^2*E(1);
for m=2:M
    A_lpha = 0;
    position=size(Selected_candidates,2)+1;
    for mm=2:M
        if Candidates(mm,1)>-1
            P_m=getcandidate3(x,z,mm,N,1000,K,L,Order);
            d(mm,1) = sum(P_m)/(1000-N);
            for r = 1: position-1
                A_lpha(mm,r) = d(mm,r)/E(1,r);
                if (r+1 < position)
                    P_rplus1 = getcandidate3(x,z,Selected_candidates(r+1),N,1000,K,L,Order);
                    d(mm,r+1) = sum(P_m.*P_rplus1)/(1000-N)-d(mm,1:r)*A(r+1,1:r)';
                else
                    d(mm,position) = sum(P_m.*P_m)/(1000-N)-A_lpha(mm,1:position-1).^2* (E(1,1:position-1)');
                    c=sum(P_m.*z(N+1:1000))/(1000-N)-A_lpha(mm,1:position-1)*(E(1,1:position-1).*G(1,1:position-1))';
                    g(mm)=c/d(mm,position);
                end
            end
        end
    end
    Q = g.^2.*d(:,position)';
    [Q_M,sorted]=sort(Q,'descend');
    next_candidate=sorted(1,1);
    if(Q_M(1,1) > (4/(1000-N))*M_S_E)
        G(1,position) = g(next_candidate);
        E(1,position)=d(next_candidate,position);
        Selected_candidates(position)=next_candidate;
        A(position,1:position-1) = A_lpha(next_candidate,1:position-1);
        M_S_E = M_S_E-Q_M(1,1);
    end
    Candidates(next_candidate,:)=-1;
end
a=a_cof(G,A);
v1=Modloutput(a,Selected_candidates,x(1,1:1000),K,L,N,Order);
mse_absolute_1st_set=(sum((v1(N+1:1000)-z(N+1:1000)).^2))/(1000-N);
mse_relative_1st_set=(100*sum((v1(N+1:1000)-z(N+1:1000)).^2))/((1000-N)*var(z(N+1:1000)));
v2=Modloutput(a,Selected_candidates,x(1,1001:2000),K,L,N,Order);
mse_absolute_2nd_set=(sum((v2(N+1:1000)-z(N+1001:2000)).^2))/(1000-N);
mse_relative_2nd_set=(100*sum((v2(N+1:1000)-z(N+1001:2000)).^2))/((1000-N)*var(z(N+1001:2000)));
v3=Modloutput(a,Selected_candidates,x(1,2001:3000),K,L,N,Order);
mse_absolute_3rd_set=(sum((v3(N+1:1000)-z(N+2001:3000)).^2))/(1000-N);
mse_relative_3rd_set=(100*sum((v3(N+1:1000)-z(N+2001:3000)).^2))/((1000-N)*var(z(N+2001:3000)));

function p = getcandidate3(x,z,t,N,N1,K,L,Order)
p=1;
s=listOrder(K,L,Order);
i=s(t,:);
j=1;
while (j<= Order*2 && i(j) ~= -1)
    if(i(j)==120)
        p = p.*x(N+1-i(j+1) : N1-i(j+1));
    else
         p = p.*z(N+1-i(j+1) : N1-i(j+1));
    end
j=j+2;
end

function z=add_noise(y,P)
mu=0;
v= P*var(y)/100;
sigma= (v)^.5;
z=zeros(1,2000);
noise=sigma*randn(1,3000)+mu;
for n=1:3000
    z(1,n)=y(1,n)+noise(1,n);
end;

function s = listOrder(K,L,Order)
M=factorial(K+L+1+Order)/(factorial(K+L+1)*factorial(Order));
s=-1*ones(floor(M),2*Order);
r=1;
for i=1:Order
    s(r+1:r+number_of_candidates_in_nthorder(K+L+1,i),1:(i*2))=candidatesof_nthorder(K,L,i);
    r=r+number_of_candidates_in_nthorder(K+L+1,i);
end;



function M=number_of_candidates_in_nthorder(N,Order)
M=factorial(N+Order-1)/(factorial(N-1)*factorial(Order));

function s=candidatesof_nthorder(K,L,Order)
N=K+L+1;
M=number_of_candidates_in_nthorder(N,Order);
s=-1*ones(M,2*Order);
r=1;
orderminus1=Order-1;
if Order==1
    s=x_and_y_terms(K,L);
   else
    a=x_and_y_terms(K,L);
    b=candidatesof_nthorder(K,L,Order-1);
    for i=1:N
        for o=(1+floor(number_of_candidates_in_nthorder(N,orderminus1))- floor(number_of_candidates_in_nthorder(N+1-i,orderminus1))): floor(number_of_candidates_in_nthorder(N,orderminus1))
            s(r,1:2)=a(i,1:2);
            s(r,3:2*Order)=b(o,:);
            r=r+1;
        end;
    end;
end

function s=x_and_y_terms(K,L)
r=0;
for i=1:K
    r=r+1;
    s(r,[1,2])=[121,i];
end;
for i=1:L+1
    r=r+1;
    s(r,[1,2])=[120,i-1];
end;

function a=a_cof(G,A)
s=size(G);
a(1:s(2))=0;
for m=1:s(2)
    v(m)=1;
    for i=m:s(2)
        if(i>m)
            v(i)=0;
            v(i)=v(i)-A(i,m:i-1)*v(m:i-1)';
        end
        a(m)=a(m)+G(i)*v(i);
    end
end


function v=Modloutput(a,P,x,K,L,N,Order)
M=size(a,2);
s=listOrder(K,L,Order);
v=zeros(1,1000);
for i=N+1:1000
    for j=1:M
        p=1;
        h=1;
        while (h<=Order*2 && s(P(j),h) ~= -1)
            if (s(P(j),h) ==120)
                p=p.*x(i-s(P(j),h+1));
            else
                 p=p.*v(i-s(P(j),h+1));
            end
            h=h+2;
        end
        v(i)=v(i)+a(j)*p;
    end
end

function modelout(sp,a,k,l,percentage,Order)
disp('-------------------------------------------------------------------------------------------------');
disp(['         ','Equation of the Identified System (when P is ',num2str(percentage),'%)  ']);
disp('-------------------------------------------------------------------------------------------------');
M=size(a,2);
c=P_candidates_strings(k,l,Order);
disp(['y[n]=',num2str(a(1,1))]);
for n=2:M
    disp(['     ','+','(',num2str(a(1,n)),')',char(c(sp(1,n),:))]);
end;
disp('-----------------------------------------------------------------------');

function s=P_candidates_strings(K,L,Order)
j=factorial(K+L+1+Order)/(factorial(K+L+1)*factorial(Order));
s=32*ones(j,7*Order);
 
r=1;
s(1,1)=49;
for i=1:Order
    s(r+1:r+number_of_candidates_in_nthorder(K+L+1,i),1:(i*7))=P_candidates_stringsin_nthorder(K,L,i);
    r=r+number_of_candidates_in_nthorder(K+L+1,i);
end;
function s= P_candidates_stringsin_nthorder(k,l,Order)
N=k+l+1;
M=number_of_candidates_in_nthorder(N,Order);
s=32*ones(M,2*Order);
r=1;
orderminus1=Order-1;
 
if Order==1
    s=P_candidates_strings_xandyterms(k,l);
else
    a=P_candidates_strings_xandyterms(k,l);
    b=P_candidates_stringsin_nthorder(k,l,orderminus1);
    for i=1:N
        for o=(1+number_of_candidates_in_nthorder(N,orderminus1)-number_of_candidates_in_nthorder(N+1-i,orderminus1)):number_of_candidates_in_nthorder(N,orderminus1)
            s(r,1:7)=a(i,1:7);
            s(r,8:7+size(b,2))=b(o,:);
            r=r+1;
        end;
    end;
end

function s=P_candidates_strings_xandyterms(k,l)
j=k+(l+1);
s=32*ones(j,7);
r=0;
for i=1:k
    r=r+1;
    if floor (i/10)>0
        s(r,5)=48+floor(i/10);
    else
        s(r,5)=32;
    end;
    s(r,[1,2,3,4,6,7])=[121,91,110,45,(48+i-(floor(i/10)*10)),93];
end;
r=r+1;
s(r,[1,2,5,7])=[120,91,110,93];
 
for i=1:l
    r=r+1;
    if floor (i/10)>0
        s(r,5)=48+floor(i/10);
    else
        s(r,5)=32;
    end;
    s(r,[1,2,3,4,6,7])=[120,91,110,45,(48+i-(floor(i/10)*10)),93];
    
end;
