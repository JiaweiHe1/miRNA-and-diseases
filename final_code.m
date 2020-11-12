%NRWRH:
function order = NRWRH( St,Sd,A,Eta,Lam,r )
size1=size(St);
size2=size(Sd);
size3=size(A);
he_t=repmat(sum(St,2),1,size1(2));
he_d=repmat(sum(Sd,2),1,size2(2));
he_a1=repmat(sum(A,2),1,size3(2));
he_a2=repmat(sum(A,1),size3(1),1);
he_at=repmat(sum(A,2),1,size1(2));
he_ad=repmat(sum(A,1),size2(1),1);
 
Mtt=(ones(size1)-min(Lam,he_at)).*St./he_t;
Mdd=(ones(size2)-min(Lam,he_ad)).*Sd./he_d;
Mtd=Lam*max(0,A./he_a1);
Mdt=Lam*max(0,A./he_a2);
Mdt=Mdt';
 
v=eye(size3(2));
u=max(0,A./he_a2);
p=[(1-Eta)*u;Eta*v];
M=[Mtt Mtd;Mdt Mdd];
size4=size(p);
c4=size4(2);
final=zeros(size4(1),1);
for i=1:c4;
    p0=p(:,i);
    p1=p0;
    p2=(1-r)*M'*p1+r*p0;
    L1norm=max(abs(p2-p1));
    while(L1norm>=10^(-10))
        p1=p2;
        p2=(1-r)*M'*p1+r*p0;
        L1norm=max(abs(p2-p1));
    end
    final=[final p2];
end
final(:,1)=[];
final=final(1:size1(1),:);
final=final-99999*A;
[x,y]=sort(final,'descend');
x1=x(1,:);
y1=y(1,:);
order=[y1;x1];
end


revised:
function W = revisedW( S )
A=sign(S);
M=99999;
S0=S+eye(size(S));
S1=S-eye(size(S));
size1=size(S);
S2=repmat(sum(S),size1(1),1);
p0=S./S2;
 
W=zeros(size1(1),1);
 
for i=1:size1(2)
    
    H=(1-S(:,i))*(1-S(:,i))'.*S0;
    f=zeros(size1(1),1);
    b=p0(:,i)'*S1(:,i);
    a=S1(:,i)';
    Aeq=ones(size1(1),1)';
    beq=1;
    lb=zeros(size1(1),1);
    ub=M*A(:,i);
    w=quadprog(H,f,-a,-b,Aeq,beq,lb,ub);
    W=[W w];
end
W(:,1)=[];
 
 
end


trueorder:
function final =trueorder( drug,target,order )
order=order';
size1=size(order);
zero=zeros(size1(1),1);
order=[zero order];
for i=1:size1(1)
    order(i,2)=target(order(i,2));
    order(i,1)=drug(i,1);
end
final=order;
end

%整合相似性矩阵:
DM=xlsread('miRNA-disease??????×?±à??.xlsx');
c=size(DM);
sizeDMh=c(1);
D=xlsread('disease??×?±à??.xlsx');
a=size(D);
sizeD=a(1);
 M=xlsread('miRNA??×?±à??.xlsx');
b=size(M);
sizeM=b(1);
FS=textread('miRNA?????à???????ó.txt');
SSD=textread('?????????à???????ó.txt');
 
IP=zeros(sizeM,sizeD);
for h=1:sizeDMh
    m=DM(h,1);
    d=DM(h,2);
    IP(m,d)=1;
end
 
FD=zeros(sizeD,1);
for d=1:sizeD
    FD(d,1)=FD(d,1)+norm(IP(:,d));
end
FD=FD./sizeD;
 
FM=zeros(sizeM,1);
for m=1:sizeM
    FM(m)=FM(m)+norm(IP(m,:));
end
FM=FM./sizeM;
 
KD=zeros(sizeD,sizeD);
for i=1:sizeD
    for j=1:sizeD
        KD(i,j)=exp(-norm(IP(:,i)-IP(:,j))/FD(i));
        KD(i,j)=max(KD(i,j),0);
    end
end
 
KM=zeros(sizeM,sizeM);
for i=1:sizeM
    for j=1:sizeM
        KM(i,j)=exp(-norm(IP(i,:)-IP(j,:))/FM(i));
        KM(i,j)=max(KM(i,j),0);
    end
end
 
SD=zeros(sizeD,sizeD);
for i=1:sizeD
    for j=1:sizeD
        if SSD(i,j)==0
            SD(i,j)=KD(i,j);
        else SD(i,j)=SSD(i,j);
        end
    end
end
 
SM=zeros(sizeM,sizeM);
for i=1:sizeM
    for j=1:sizeM
        if FS(i,j)==0
            SM(i,j)=KM(i,j);
        else SM(i,j)=FS(i,j);
        end
    end
end
