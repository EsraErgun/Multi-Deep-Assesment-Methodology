function [niumin1,MAPE,MAPE1,A,B,C,Ffitting]=fitting(fileID,fdata1,l1,M1,alfa1,m1,an)

%alfa=alfa1;
% fdata=zeros(13,1);
% for i=1:13
%   fdata(i)=sin(i*pi/12)*cos(i*pi/12);  
% end
% fdata=[
% -0.057500066	0.064477221	0.057632061	0.097959364	0.089923049	0.069479881	0.050435079	0.008150246	-0.048231539	0.084271043	0.112001106	0.047885	0.084858	0.049397	0.060845	0.033231	0.07502	0.02959
% 
% ];

fdata=fdata1;
P=zeros(length(fdata),1);

P0=fdata;
P=P0;
j=length(fdata);

l=l1;
M=M1;
m=m1;
alfar=ones(m,1);

 fileID1 = fopen('alfar.dat','r');
 
 for i=1:m
 
     line = fgets(fileID1);
   alfar(i)= sscanf(line,'%f'); 
 
 end
 fclose(fileID1);
  alfar;


finalF=zeros(j,1);

%alfa=zeros(m,1);

MD=1+(l)*m;
%MD=10;

NONIU=100;
niu5=linspace(0.01,1,NONIU);

DDD=zeros(l,NONIU,j-l,m);
%alfa=0.01;
%FFF=zeros(l,2*M,NONIU,j-l,m);
 for i1=1:NONIU
   
for r=1:m   
   % alfa=niu5(i1);
 for k=1:l
   for i=l+1:j
   % DDD(k,i1,i-l,r)=0;
    DDD(k,i1,i-l,r)=an(r,1);
 for n=1:M   
   
 
  DDD(k,i1,i-l,r)=  DDD(k,i1,i-l,r)+ Dkrn(k,n,niu5(i1),i,alfar(r))*an(r,n+1); 
 
 end
 
     
 end
 end
end
 end  
   
min=1000000000000000000000000000;
for i1=1:NONIU
    niu=niu5(i1);
  A=zeros(MD,MD);
C=zeros(MD,1); 
    A(1,1)=j-l; 
    %%%%%First ROW%%%%%%%%%%%%%%%%
   
    icol=1;
   
            for r=1:m
                for k=1:l
                    
                 %   for n=1:M
               icol=icol+1;  
            A(1,icol)=0;
            for i=l+1:j
              
              A(1,icol)=A(1,icol)+DDD(k,i1,i-l,r);
        
            end
           % end
            end
            end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %%%%%%%%%First colomn%%%%%%%%%%%%%%%%5
  
     irow=1;
   
      
            for r=1:m
                for k=1:l
                  %  for n=1:M
               irow=irow+1;  
            A(irow,1)=0;
            for i=l+1:j
              
              A(irow,1)=A(irow,1)+DDD(k,i1,i-l,r);
        
            end
            
           % end
            end
            end        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  
   C(1)=0;
   for i=l+1:j
   C(1)=C(1)+P(i);    
   end
   
    %%%%%%%%%%%%%%%%%%%%%danarcheni%%%%%%%%%%
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%row 1->1+l*(2*M)*m %%%%%%%%%%%%%%%%%%%%%%%
     row=1;
    for rrow=1:m
        for krow=1:l
         % for nrow=1:M 
          row=row+1;
          col=1;
          
          for rcol=1:m
            for kcol=1:l
             % for ncol=1:M  
                 col=col+1  ;
                for i=l+1:j
              
                 A(row,col)=A(row,col)+DDD(kcol,i1,i-l,rcol)*DDD(krow,i1,i-l,rrow);
             
                end         
          end
              
          end  
          
          
             C(row)=0;
   for i=l+1:j
   C(row)=C(row)+P(i)*DDD(krow,i1,i-l,rrow);    
   end  
        
        end
    end
           
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
           
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% B=pcg(A'*A,A'*C,1e-10);
% B=(A'*A)\(A'*C);
% B=pinv(A)*C;
B=A\C;
% B= linsolve(A,C);
%  A\eye(size(A))


F=zeros(j,1);

for i=l+1:j
   F(i)=B(1);
 
   irow=1;   
 
   for r=1:m
     for k=1:l  
      %  for n=1:M   
       irow=irow+1;  
      
    F(i)=F(i)+DDD(k,i1,i-l,r)*B(irow);
       
       
       % end
     end
   end
  
end

for i=1:l
   F(i)=P(i);
end
eps=0;
for i=l+1:j
if P(i)==0
   % eps=eps+abs( (P(i)-F(i))/(P(i)+0.0001) );
else
eps=eps+abs( (P(i)-(F(i)) )/abs(P(i)));

end

end
%eps=eps*100/(j-k-1);
eps=eps*100/(j-l);
%size(B)

if(eps<min)
    min=eps;
    niumin=niu;
    minerror=min;
  %  finalF=F;
    Bfinal=B;
    mini1=i1;
end

end

 minerror; %%here we print the optimal error
 niumin; %% here we print the optimal niu
% 
   x=linspace(1,j,j);
    Psm=zeros(j,1);
    for i=1:j
        Psm(i)=P0(i);
        end
   F11=zeros(j,1);
  

   for i=1:l
 
    F11(i)=P(i);
   end

   for i=l+1:j
     
       F11(i)=Bfinal(1);
 
   irow=1;
     
    for r=1:m
     for k=1:l  
       % for n=1:M   
       irow=irow+1;  
      
    F11(i)=F11(i)+DDD(k,mini1,i-l,r)*Bfinal(irow);
               
       % end
     end
    end
    
end

 for i=1:1+l*m
 fprintf(fileID,'%40.20f\n',Bfinal(i));
  end

%AAA=Bfinal;
niumin1=niumin;
% figure(5)
%    plot(x,Psm) %here we draw the modeled fucntion
%     hold on;
%    plot(x,F11) % here we draw the initial function
%    grid on
%    legend('P','F')
   Ffitting=F11';
   MAPE=0;
for i=l+1:j
    MAPE=MAPE+abs(Psm(i)-F11(i))/abs(Psm(i))*100/(j-l);
   MAPE1(i)=abs(Psm(i)-F11(i))/abs(Psm(i))*100/(j-l);

end
%  xlswrite('rate_modeling_pinv.xlsx',MAPE, 'C1:C1')
end
 