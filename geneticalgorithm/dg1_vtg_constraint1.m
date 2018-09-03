function[Pl,val,CVD,losses]=dg1_vtg_constraint1(l,n)
c=0;
R=16;
vmin=zeros(R,1);
vmax=zeros(R,1);
CVD=zeros(R,1);
v=zeros(15,1);
losses=zeros(R,1);
[K L]=size(n);
val=zeros(K,5);
Pl=zeros(K,1);
Z1_2=1.35309+1.32349i;
Z2_3=1.17024+1.14464i;
Z3_4=0.84111+0.82271i;
Z4_5=1.52348+1.02760i;
Z4_14=2.23081+1.50470i;
Z4_15=1.19702+0.80740i;
Z3_11=1.79553+1.21110i;
Z11_12=2.44845+1.65150i;
Z12_13=2.01317+1.35790i;
Z2_9=2.01317+1.35790i;
Z9_10=1.68671+1.13770i;
Z6_8=1.25143+0.84410i;
Z6_7=1.08820+0.73400i;
Z2_6=2.55727+1.72490i;
S_load=[(44100+44990i);
  (70000+71414i);
  (140000+142828i);
  (44100+44990i);
  (70000+71414i);
  (140000+142828i);
  (140000+142828i);
  (70000+71414i);
  (44100+44990i);
  (140000+142828i);
  (70000+71414i);
  (44100+44990i);
  (70000+71414i);
  (140000+142828i)];

pf=0.85; % power factor
for h=1:K
    k=n(h);
    x=0;
S=zeros(14,1);
I=zeros(14,1);
S_gen=zeros(14,1);
S_gen(l(h)-1)=k;
S_p=S_load;
S=S_p-S_gen;

for i=1:15
    V(i)=11000+0i;
end

while x==0
    c=c+1;
for i=2:15
    I(i)=conj(S(i-1)/V(i));
end

%backward propagation
I4_5=I(5);
I4_14=I(14);
I4_15=I(15);
I3_4=I(4)+I(15)+I(14)+I(5);
I12_13=I(13);
I11_12=I(12)+I(13);
I3_11=I(11)+I11_12;
I2_3=I(3)+I3_11+I3_4;
I6_8=I(8);
I6_7=I(7);
I2_6=I(6)+I(7)+I(8);
I9_10=I(10);
I2_9=I(9)+I(10);
I1_2=I(2)+I2_9+I2_6+I2_3;

%forward propagation
V(2)=V(1)-((I1_2)*Z1_2);
V(3)=V(2)-((I2_3)*Z2_3);
V(4)=V(3)-((I3_4)*Z3_4);
V(5)=V(4)-((I4_5)*Z4_5);
V(6)=V(2)-((I2_6)*Z2_6);
V(8)=V(6)-((I6_8)*Z6_8);
V(7)=V(6)-((I6_7)*Z6_7);
V(9)=V(2)-((I2_9)*Z2_9);
V(10)=V(9)-((I9_10)*Z9_10);
V(11)=V(3)-((I3_11)*Z3_11);
V(12)=V(11)-((I11_12)*Z11_12);
V(13)=V(12)-((I12_13)*Z12_13);
V(14)=V(4)-((I4_14)*Z4_14);
V(15)=V(4)-((I4_15)*Z4_15);
          

for i=2:15
    Snew(i-1)=V(i)*conj(I(i));   % power mismatch method
    A(i-1)=S(i-1)-Snew(i-1);
    end
if(max(A)>0.00001)
   x=0;  
    else
    x=1;
end
end
%  I1_2
[th12,I12]=cart2pol(real(I1_2),imag(I1_2));
% I2_3
 [th23,I23]=cart2pol(real(I2_3),imag(I2_3));
% I3_4
 [th34,I34]=cart2pol(real(I3_4),imag(I3_4));
% I4_5
[th45,I45]=cart2pol(real(I4_5),imag(I4_5));
% I4_14
[th414,I414]=cart2pol(real(I4_14),imag(I4_14));
% I4_15
[th415,I415]=cart2pol(real(I4_15),imag(I4_15));
% I3_11
[th311,I311]=cart2pol(real(I3_11),imag(I3_11));
% I11_12
[th1112,I1112]=cart2pol(real(I11_12),imag(I11_12));
% I12_13
[th1213,I1213]=cart2pol(real(I12_13),imag(I12_13));
% I2_9
[th29,I29]=cart2pol(real(I2_9),imag(I2_9));
% I9_10
[th910,I910]=cart2pol(real(I9_10),imag(I9_10));
% I6_8
[th68,I68]=cart2pol(real(I6_8),imag(I6_8));
% I6_7
[th67,I67]=cart2pol(real(I6_7),imag(I6_7));
% I2_6
[th26,I26]=cart2pol(real(I2_6),imag(I2_6));
% %%POWER loss in each branch

 P(1)=((I12))^2*(real(Z1_2));
 P(2)=((I23))^2*(real(Z2_3));
 P(3)=((I34))^2*(real(Z3_4));
 P(4)=((I45))^2*(real(Z4_5));
 P(5)=((I414))^2*(real(Z4_14));
 P(6)=((I415))^2*(real(Z4_15));
 P(7)=((I311))^2*(real(Z3_11));
 P(8)=((I1112))^2*(real(Z11_12));
 P(9)=((I1213))^2*(real(Z12_13));
 P(10)=((I29))^2*(real(Z2_9));
 P(11)=((I910))^2*(real(Z9_10));
 P(12)=((I68))^2*(real(Z6_8));
 P(13)=((I67))^2*(real(Z6_7));
 P(14)=((I26))^2*(real(Z2_6));
the=zeros(15,1);
v=zeros(15,1);
for i=1:15
     [the(i),v(i)]=cart2pol(real(V(i)),imag(V(i)));
end
%  v
norm_vtg=v/11000;
 vmin(h,1)=min(v);
 vmax(h,1)=max(v);
 c=0;
 for i=1:15
     c=c+abs(1-norm_vtg(i));
 end
 CVD(h,1)=c;
P_loss=0;
 for i=2:15
 P_loss=P_loss+P(i-1); %total power loss
 end
 Pl(h)=P_loss;
 val(h,:)=[l(h) n(h) vmin(h) vmax(h) Pl(h)];
end
losses=Pl()./62162;

