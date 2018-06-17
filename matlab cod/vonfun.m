
%function von_result=vonfun(signal,epsilong0,t_size);
function [von_result,lamda]=vonfun(signal,epsilong0,t_size);
%%计算矩阵A
I=eye(t_size);
A=zeros(t_size-3,t_size);
for i=1:t_size-3
    A(i,i)=-1;
    A(i,i+1)=3;
    A(i,i+2)=-3;
    A(i,i+3)=1;
end
format long;
A_s=sparse(A);
i=1;
epsilong=epsilong0;
%epsilong=1;
lamda=1/epsilong;
%%计算(I(n)+λ^2??A^T (n)A(n) )??Y^'=Y
temp1=inv(I+lamda*A_s'*A_s);
von_result=temp1*signal;
%%计算最佳平滑因子
P1=epsilong*eye(t_size);
P2=eye(t_size-3);
P1_s=sparse(P1);
P2_s=sparse(P2);
P=[P1_s,zeros(t_size,(t_size-3));zeros((t_size-3),t_size),P2_s];
P_s=sparse(P);
v1=von_result-signal;
v2=A_s*von_result;
B1=eye(t_size);
B1_s=sparse(B1);
B2=A;
B2_s=sparse(B2);
B=[B1;B2];
B_s=sparse(B);
H1=B1_s'*P1_s*B1_s;
H2=B2_s'*P2_s*B2_s;
H=B_s'*P_s*B_s;
H1_s=sparse(H1);
H2_s=sparse(H2);
H_s=sparse(H);
sigma1=v1'*P1*v1/(t_size-trace(inv(H_s)*H1_s));
sigma2=v2'*P2*v2/((t_size-3)-trace(inv(H_s)*H2_s));
P1_t=sigma2/sigma1*P1;
k=0;
while(abs(1-k)>0.0001)
k=sigma2/sigma1;
epsilong=P1_t(1,1);
lamda=1/epsilong;
temp1=inv(I+lamda*A_s'*A_s);
von_result=temp1*signal;
v1=von_result-signal;
v2=A_s*von_result;
P1=epsilong*eye(t_size);
P1_s=sparse(P1);
P=[P1_s,zeros(t_size,(t_size-3));zeros((t_size-3),t_size),P2_s];
P_s=sparse(P);
H1=B1_s'*P1_s*B1_s;
H=B_s'*P*B_s;
H1_s=sparse(H1);
H_s=sparse(H);
sigma1=v1'*P1*v1/(t_size-trace(inv(H_s)*H1_s));
P1_t=sigma2/sigma1*P1;
i=i+1;
end
 
%%计算最终平滑序列
lamda=1/epsilong;
temp1=inv(I+lamda*A_s'*A_s);
von_result=temp1*signal;