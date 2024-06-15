function [output]=F1(fileID,M1,data2,alfa)


% data1=[
% -0.057500066	0.064477221	0.057632061	0.097959364	0.089923049	0.069479881	0.050435079	0.008150246	-0.048231539	0.084271043	0.112001106	0.047885	0.084858	0.049397	0.060845	0.033231	0.07502	0.02959
% 
% ];
data1=data2;
P0=data1;

%alfa=1;
m=length(P0);

P=zeros(m,1);
 
 %k=k1;
M=M1;
  
 for i=1:m
     P(i)=P0(i);
 end

%finalF=zeros(m,1);








MD=1+M;
%MD=10;

A=zeros(MD,MD);
C=zeros(MD,1);

NONIU=100;
niu1=linspace(0.01,1,NONIU);
     ccc=zeros(M,NONIU,m);
     for i1=1:NONIU
         i1;
         alfa=niu1(i1);
        for t=1:M
         for i=1:m
             ccc(t,i1,i)=CC(t,niu1(i1),i,alfa);
         end
     end
     end
min=100000;
for i1=1:NONIU
    niu=niu1(i1);
 
    
    %%%%%First ROW%%%%%%%%%%%%%%%%
    A(1,1)=m;
   
    icol=1;
      for t=1:M
   
      
            icol=icol+1;
            A(1,icol)=0;
            for i=1:m
               
              A(1,icol)=A(1,icol)+ccc(t,i1,i)  ;
               
            end
         
      end
    
      
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    

    
    %%%%%%%%% first scolomn%%%%%%%%%%%%%%%%
       irow=1;
        for d=1:M
   
       
            irow=irow+1;
            A(irow,1)=0;
            for i=1:m
              
              A(irow,1)=A(irow,1)+ccc(d,i1,i)  ;
               
            end
            
        end
  
    
      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  
  
    
    %%%%%%%%%%%%%%%%%%%%%Matrix%%%%%%%%%%
    
        irow=1;
         for t=1:M
      %t
        
            irow=irow+1;
            icol=1;
         for d=1:M    
    
      
            icol=icol+1;
            A(irow,icol)=0;
            for i=1:m
              
              A(irow,icol)=A(irow,icol)+ccc(t,i1,i)*ccc(d,i1,i) ;
              
            end
           
     
        end
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    

   %%%%%%%%%%%%%%%%%% C matrix%%%%%%%%%%%%%%%%%%%%%
   C(1)=0;
   for i=1:m
   C(1)=C(1)+P(i);    
   end
   
  
   
   icol=1;
    for t=1:M
   
      
           icol=icol+1;
          C(icol)=0;
   for i=1:m
      
   C(icol)=C(icol)+P(i)*ccc(t,i1,i);  
     
   end  
      
     
    end
   
   
     
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
% B=pcg(A'*A,A'*C,1e-10);
% B=(A'*A)\(A'*C);
% B= linsolve(A,C);

B=A\C;
% B=pinv(A)*C;
F=zeros(m,1);

for i=1:m
   F(i)=B(1);
  
   irow=1;
      for t=1:M
  
    
   irow=irow+1;
  
    F(i)=F(i)+ccc(t,i1,i)*B(irow);
 
    
      end
   
      
  
end



eps=0;

for i=1:m
if P(i)==0
 % eps=eps+abs((P(i)-F(i))/(P(i)+0.001));  
else
eps=eps+abs((P(i)-F(i))/P(i)+0.00001);
end
end
eps=eps*100/m;

if(eps<min)
    min=eps;
    niumin=niu;
    minerror=min;
   % finalF=F;
    Bfinal=B;
end

end

minerror;
niumin;

   x=linspace(1,m,m);
    Psm=zeros(m,1);
   F11=zeros(m,1);
   for i=1:m
    
   
   Psm(i)=P0(i);
   end


   
   
   for i=1:m
     
   F11(i)=Bfinal(1);
   irow=1;
      for t=1:M
  
    
   irow=irow+1;
  
    F11(i)=F11(i)+CC(t,niumin,i,niumin)*Bfinal(irow);
 
       end
     
   Ff1=F11;
     
   
 

end
 %figure(1001)
   % plot(x,F11)
    % hold on;
   % plot(x,Psm)



fprintf(fileID,'%12.8f\n',niumin);
output=Bfinal;
%Bfinal(1)
end
