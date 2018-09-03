clear all
clc
disp('CODE FOR FUNCTION WITH 2 OBJECTIVES: LOSSES AND CVD')
disp('---------------------------------------------------')
R=16;
losses=zeros(R,1);
god=zeros(40,1);
CVD=zeros(R,1);
don=zeros(R,1);
vmin=zeros(R,1);
vmax=zeros(R,1);
cost=zeros(R,1);
lol=zeros(R,5);
Po=zeros(R,1);
I=zeros(14,1);
I1=zeros(14,1);
P_b=zeros(14,1);
P_b1=zeros(14,1);
pf=0.85;
a=zeros(40,1);
nvar=2;
nbits=4;
N=nvar*nbits; 
Nt=N;% number of bits in a chromosome
M=R;
popsize=M;% number of chromosomes must be even
last=10; % number of generations
sel=0.5; % selection rate
ind=zeros(M,1);
maxit=0;
M2=2*ceil(sel*M/2);
keep=M2;% number of chromosomes kept
mutrate=0.1; % mutation rate
nmuts=ceil(mutrate*N*(M-1)); % number of mutations
% creates M random chromosomes with N bits
pop=round(rand(M,N)); % initial population 
bi=zeros(R,1);
 bi=bi2de(pop(:,(1:4)));
 
  for i=1:M
     bin=bi2de(pop(:,(5:8)));
  end
   for i=1:M
 while bi(i,1)==1 | bi(i,1)==0 | bi(i,1)==15 
       pop(i,(1:4))=round(rand(1,4));
         bi(i)=bi2de(pop(i,(1:4)));
 end
   end
      for i=1:R
 while  bin(i,1)<4 | bin(i,1)>12
      pop(i,(5:8))=round(rand(1,4));
      bin(i)=bi2de(pop(i,(5:8)));
 end
  end
 for i=1:R
     Po(i)=((10^5)*bin(i));
 end
 while maxit<20
lower=10450;
upper=11550;

     n=zeros(R,1);
Q=zeros(R,1);
Q=(tan(acos(pf)))*Po;
n=Po+j.*Q;
for i=1:M
 while bi(i,1)==1 | bi(i,1)==0 | bi(i,1)==15 
       pop(i,(1:4))=round(rand(1,4));
         bi(i)=bi2de(pop(i,(1:4)));
 end
end   
for i=1:M
 while bi(i,1)==1 | bi(i,1)==0 | bi(i,1)==15 
       pop(i,(1:4))=round(rand(1,4));
         bi(i)=bi2de(pop(i,(1:4)));
 end
     end
 for i=1:R
 while bin(i,1)<4 | bin(i,1)>12
      pop(i,(5:8))=round(rand(1,4));
      bin(i)=bi2de(pop(i,(5:8)));
 end
 end
 for i=1:R
     Po(i)=((10^5)*bin(i));
 end
n=zeros(R,1);
Q=zeros(R,1);
Q=(tan(acos(pf)))*Po;
n=Po+j.*Q;
 for i=1:R
    while lol(i,5)> 62126
         pop(i,(1:4))=round(rand(1,4));
             bi(i,1)=bi2de(pop(i,(1:4)));
              while bi(i,1)==1 | bi(i,1)==0 | bi(i,1)==15 
                    pop(i,(1:4))=round(rand(1,4));
             bi(i,1)=bi2de(pop(i,(1:4)));
              end
             [cost(i,1),lol(i,1:5)]=dg1_vtg_constraint1(bi(i,1),n(i,1));
    
    end
 end

for i=1:M
 while bi(i,1)==1 | bi(i,1)==0 | bi(i,1)==15 
       pop(i,(1:4))=round(rand(1,4));
         bi(i)=bi2de(pop(i,(1:4)));
 end
end
    
       
[cost,lol,CVD,losses]=dg1_vtg_constraint1(bi,n); % calculates population
     don=zeros(R,1);
     don=zeros(R,1);
    for i=1:R
        while lol(i,3)<lower | lol(i,4)>upper
            pop2(i,1:4)=round(rand(1,4));
            bi(i,1)=bi2de(pop2(i,1:4));
     while bi(i,1)==1 | bi(i,1)==0 | bi(i,1)==15 
             pop(i,(1:4))=round(rand(1,4));
             bi(i,1)=bi2de(pop(i,(1:4)));
     end
              [cost(i,1),lol(i,1:5),CVD,losses]=dg1_vtg_constraint1(bi(i,1),n(i,1));
        end             
     end      
          for i=1:M
 while bi(i,1)==1 | bi(i,1)==0 | bi(i,1)==15 
       pop(i,(1:4))=round(rand(1,4));
         bi(i)=bi2de(pop(i,(1:4)));
 end
          end
    
[cost,lol,CVD,losses]=dg1_vtg_constraint1(bi,n);
ff=(0.5.*CVD)+(0.5.*losses);                                 % cost using ff                       
 [cost1,ind]=sort(ff);% min cost in element 1
  par=bi(ind,:);
  bi=par;
  non=n(ind,:);
  n=non;
  lol=lol(ind,:);
  CVD=CVD(ind,:);
  figure(maxit+1)
     for i=1:R
     plot(par(i),real(n(i)),'*r')
     xlabel('Bus no.')
     ylabel('Size of DG')
     hold on
     end
 pop1=pop(ind,:); % sorts population with
% %                   % lowest cost first
pop=pop1;
  M1=ceil((popsize-keep)/2); % number of matings
prob=flipud([1:keep]'/sum([1:keep]));% weights
                              % chromosomes based
                               % upon position in
                                % list
odds=[0 cumsum(prob(1:keep))']; % probability

pick1=rand(1,M1); % mate #1
pick2=rand(1,M1); % mate #2
% ma and pa contain the indicies of the chromosomes
% that will mate
ic=1;
while ic<=M1
for id=2:keep+1
if pick1(ic)<=odds(id) & pick1(ic)>odds(id-1)
ma(ic)=id-1;
end % if
if pick2(ic)<=odds(id) & pick2(ic)>odds(id-1)
pa(ic)=id-1;
end % if
end % for
ic=ic+1;
end % while
ma1=ma';
pa1=pa';
for ix=1:(keep/2); % index of mate #1
xp=ceil(rand(1,M1)*(Nt-1)); % crossover point
% pop
pop(keep+ix,1:xp)=pop(ma1(ix),1:xp);
pop(keep+ix,xp+1:Nt)=pop(pa1(ix),xp+1:Nt);
%                first offspring
pop(keep+ix+1,1:xp)=pop(pa1(ix),1:xp);
pop(keep+ix+1,xp+1:Nt)= pop(ma1(ix),xp+1:Nt);
end
%                 second offspring
%_____________________________________

% Mutate the population
nmut=ceil((popsize-1)*Nt*mutrate); % total number
%                                    of mutations
mrow=ceil(rand(1,nmut)*(popsize-1))+1; % row to mutate
mcol=ceil(rand(1,nmut)*Nt); % column to mutate
for ii=1:nmut
pop(mrow(ii),mcol(ii))=abs(pop(mrow(ii),mcol(ii))-1);
% toggles bits
end % ii
% pop=pop1;
bi=bi2de(pop(:,1:4));
bin=bi2de(pop(:,5:8));
for i=1:R
     Po(i)=((10^5)*bin(i));
 end
n=zeros(R,1);
Q=zeros(R,1);
Q=(tan(acos(pf)))*Po;
n=Po+j.*Q;
 
 maxit=maxit+1;

          [god(maxit),lol1,CVD,losses]=dg1_vtg_constraint1(bi(1),n(1));
      [maxit  cost1(1)]
      a(maxit)=cost1(1);
      figure(41)
         plot(maxit,cost1(1),'*k')
   hold on
   
 xo=0;
   end
 [loss,lol1,CVD,losses]=dg1_vtg_constraint1(bi(ind(1)),n(1));
 
disp(['Location at which DG should be connected : ',num2str(bi(ind(1)))])
disp(['Size of DG : ',num2str(n(1))])
disp(['Vmin=',num2str(lol1(1,3))])
disp(['Vmax=',num2str(lol1(1,4))])
[V,kk,I,P_b]=mega(bi(ind(1)),n(1));
figure(42)

disp(['Losses with DG : ',num2str(kk)])
[V1,kk1,I1,P_b1]=mega(bi(ind(1)),0);
disp(['Losses without DG : ',num2str(kk1)])
disp('Voltage profile without DG')
V1'
disp('Voltage profille with DG')
V'
disp('Current without DG   Current with DG  Power loss without DG    Power loss with DG')
[I1 I P_b1 P_b]
figure(42)
for i=1:15
 plot(i,V(i),'>g')
 hold on
end
for i=1:15
  plot(i,V1(i),'>m')
 hold on
end
lower=10450;
upper=11550;
figure(42)
for i=0:0.01:15
plot(i,upper,'-m')
hold on
plot(i,lower,'-r')
hold on
end
xlabel('Bus no.')
ylabel('Voltage')


 