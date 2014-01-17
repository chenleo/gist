
%%This program uses the positive and negative samples to train Fisher model, and predicts the piRNA.
%%you need put three files "positive173090.txt,negative193321.txt,prediction.txt" into the current work folder, and run "predicting.m".
%%%This program consists of 3 parts: part 1)transform the FASTA sequences into 1364-D vectors. The FASTA sequences include positive training set
%%"positive173090.txt"(transformed from Line 12 to Line 177), negative training set "negative193321.txt°±(transformed Line 180-281),
%%%%sequences to be predicted "prediction.txt"(transformed from Line 290 to Line 410).
%%%%                                 part 2)train the Fisher Discriminant Model based on 120000 pairs of
%%%%randomly selected vectors from positive and negative sets (Line 412 to 560)
%%%                                   part 3)predict, output the predicted piRNA as "predictedpiRNA.txt" (Line 565 to END)                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%% part 1: transform seuqences into vectors%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%read positive sequences%%%%%%%%%%%%%%%%%%%%%%
Pxulie=cell(1,173090);  %%%% positive sequences
Pshuo=cell(1,173090);   %%%% the illustration of sequences

fid1=fopen('positive173090.txt','r');  
j=0; 
while feof(fid1)==0               % judge if reach the end of the file, otherwise, do the loop.
    
   f=fgetl(fid1);                % read the title line
   ss=size(f);
   s=ss(1,2);
            if (s>0)&&(f(1)=='>')  
               ff1=fgetl(fid1); 
               j=j+1;
               Pshuo{1,j}=f; 
               Pxulie{1,j}=ff1;
            end 
end
sta1=fclose(fid1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%generate K-mer strings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str2=cell(1,16);  %%%%
mu='AGCT';
geshu=0; %%%% 
          for i=1:4
               for j=1:4               
                     geshu=geshu+1;
                     str2{1,geshu}=strcat(mu(1,i), mu(1,j));
             
               end
         end 

str3=cell(1,64);  %%%%
mu='AGCT';
geshu=0; %%%% 
        for i=1:4
              for j=1:4
                    for k=1:4
                        
                            geshu=geshu+1;
                           str3{1,geshu}=strcat(mu(1,i), mu(1,j),mu(1,k)); 
                    end
              end 
       end


str4=cell(1,256);  %%%%
mu='AGCT';
geshu=0;          %%%%
for i=1:4
      for j=1:4
         for k=1:4
             for m=1:4
                 
             
                   geshu=geshu+1;
                   str4{geshu}=strcat(mu(1,i), mu(1,j),mu(1,k),mu(1,m));
                 
             end
         end
      end 
end


str5=cell(1,1024);  
mu='AGCT';
geshu=0; 
for i=1:4
      for j=1:4
         for k=1:4
             for m=1:4
                 for n=1:4
             
                           geshu=geshu+1;
                           str5{geshu}=strcat(mu(1,i), mu(1,j),mu(1,k),mu(1,m),mu(1,n));
                 end
             end
         end
      end 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fidxie = fopen('vectorpositive173090.txt','w');
xun=173090; 
 

for sui=1:xun  
    Plianall=zeros(1,1364);    
    ff=Pxulie{1,sui}; 
    ss=size(ff);
    s=ss(1,2);      %%%length
    
                %%%%%%%%1-mer frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
                 bases=basecount(ff);%%%%%%bases.A, bases.G, bases.C, bases.T 
                  
                 Plianall(1,1)=(bases.A)/s;
                 Plianall(1,2)=(bases.G)/s;
                 Plianall(1,3)=(bases.C)/s;
                 Plianall(1,4)=(bases.T)/s;
                %%%%%%%%2-mer%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                for j=1:2:s-1  
                   for k=1:16 
                      if ff(1,j:j+1)==str2{1,k}  
                           Plianall(1,k+4)=Plianall(1,k+4)+1;
                      end
                   end
                end  
                
                Plianall(1,5:20)=Plianall(1,5:20)/(s/2);
                
    
             %%%%%%%3-mer frequency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                for j=1:3:s-2  %%%
                   for k=1:64 
                      if ff(1,j:j+2)==str3{1,k}  
                           Plianall(1,k+20)=Plianall(1,k+20)+1;
                      end
                   end
                end  %%%
                
                 Plianall(1,21:84)=Plianall(1,21:84)/(s/3);
                
                %%%%%%%%4-mer frequency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                for j=1:4:s-3  
                   for k=1:256
                      if ff(1,j:j+3)==str4{1,k}  %%%%%%%%%%%%%%
                           Plianall(1,k+84)=Plianall(1,k+84)+1;  %%%%
                      end
                   end
                end  %%%
                   Plianall(1,85:340)=Plianall(1,85:340)/(s/4);
                   
                   
                 %%%%%%%%5-mer frequency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                for j=1:5:s-4  %%%%
                   for k=1:1024
                      if ff(1,j:j+4)==str5{1,k}  %%%%%%%%%%%%%%
                           Plianall(1,k+340)=Plianall(1,k+340)+1;  %%%%
                      end
                   end
                end  %%%
           
                
                  Plianall(1,341:1364)=Plianall(1,341:1364)/(s/5);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 zhuan=Plianall(1,:);  %%%  1X1364
 
 
 for bao=1:1364
     fprintf(fidxie,'%19.18f ',zhuan(1,bao)); 
 end 

fprintf(fidxie,'\n'); 
    
end 
staxie=fclose(fidxie);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the end of positive %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%read negative sequences%%%%%%%%%%%%%%%%%%%%%%
Nxulie=cell(1,193321);  %%%% negative sequences
Nshuo=cell(1,193321);   %%%% the illustration of sequences

fid1=fopen('negative193321.txt','r');  
j=0; 
while feof(fid1)==0               % judge if reach the end of the file, otherwise, do the loop.
    
   f=fgetl(fid1);                % the title line
   ss=size(f);
   s=ss(1,2);
            if (s>0)&&(f(1)=='>')  
               ff1=fgetl(fid1); 
               j=j+1;
               Nshuo{1,j}=f; 
               Nxulie{1,j}=ff1;
            end 
end
sta1=fclose(fid1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fidxie = fopen('vectornegative193321.txt','w');
xun=193321; 
 

for sui=1:xun  
     Plianall=zeros(1,1364);   
    ff=Nxulie{1,sui}; 
    ss=size(ff);
    s=ss(1,2);      %%%length
    
                %%%%%%%%1-mer frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
                 bases=basecount(ff);%%%%%%bases.A, bases.G, bases.C, bases.T 
                  
                 Plianall(1,1)=(bases.A)/s;
                 Plianall(1,2)=(bases.G)/s;
                 Plianall(1,3)=(bases.C)/s;
                 Plianall(1,4)=(bases.T)/s;
                %%%%%%%%2-mer%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                for j=1:2:s-1  
                   for k=1:16 
                      if ff(1,j:j+1)==str2{1,k}  
                           Plianall(1,k+4)=Plianall(1,k+4)+1;
                      end
                   end
                end  
                
                Plianall(1,5:20)=Plianall(1,5:20)/(s/2);
                
    
             %%%%%%%3-mer frequency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                for j=1:3:s-2  %%%
                   for k=1:64 
                      if ff(1,j:j+2)==str3{1,k}  
                           Plianall(1,k+20)=Plianall(1,k+20)+1;
                      end
                   end
                end  %%%
                
                 Plianall(1,21:84)=Plianall(1,21:84)/(s/3);
                
                %%%%%%%%4-mer frequency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                for j=1:4:s-3  
                   for k=1:256
                      if ff(1,j:j+3)==str4{1,k}  %%%%%%%%%%%%%%
                           Plianall(1,k+84)=Plianall(1,k+84)+1;  %%%%
                      end
                   end
                end  %%%
                   Plianall(1,85:340)=Plianall(1,85:340)/(s/4);
                   
                   
                 %%%%%%%%5-mer frequency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                for j=1:5:s-4  %%%%
                   for k=1:1024
                      if ff(1,j:j+4)==str5{1,k}  %%%%%%%%%%%%%%
                           Plianall(1,k+340)=Plianall(1,k+340)+1;  %%%%
                      end
                   end
                end  %%%
           
                
                  Plianall(1,341:1364)=Plianall(1,341:1364)/(s/5);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 zhuan=Plianall(1,:);  %%%  1X1364
 
 
 for bao=1:1364
     fprintf(fidxie,'%19.18f ',zhuan(1,bao)); 
 end 

fprintf(fidxie,'\n'); 
    
end 
staxie=fclose(fidxie);  %%% 


%%%%%%%%%%%%





%%%%%%%%%%%%%%%read prediction sequences%%%%%%%%%%%%%%%%%%%%%%


fid1=fopen('prediction.txt','r');  
j=0; %%% how many sequences?
while feof(fid1)==0               % judge if reach the end of the file, otherwise, do the loop.
    
   f=fgetl(fid1);                % read the next line to the title line
   ss=size(f);
   s=ss(1,2);
            if (s>0)&&(f(1)=='>')  
               ff1=fgetl(fid1); 
               j=j+1;
           end 
end
sta1=fclose(fid1);

jzhong=j; %% there are izhong sequences to be predicted

prelie=cell(1,jzhong);  %%%% prediction sequences
preshuo=cell(1,jzhong);   %%%% the illustration of sequences

fid1=fopen('prediction.txt','r');  
j=0; 
while feof(fid1)==0               % judge if reach the end of the file, otherwise, do the loop.
    
   f=fgetl(fid1);                % the title line
   ss=size(f);
   s=ss(1,2);
            if (s>0)&&(f(1)=='>')  
               ff1=fgetl(fid1); 
               j=j+1;
               preshuo{1,j}=f; 
               prelie{1,j}=ff1;
            end 
end
sta1=fclose(fid1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fidxie = fopen('vectorprediction.txt','w');
xun=jzhong; 


for sui=1:xun  
    Plianall=zeros(1,1364);    
    ff=prelie{1,sui}; 
    ss=size(ff);
    s=ss(1,2);      %%%length
    
                %%%%%%%%1-mer frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
                 bases=basecount(ff);%%%%%%bases.A, bases.G, bases.C, bases.T 
                  
                 Plianall(1,1)=(bases.A)/s;
                 Plianall(1,2)=(bases.G)/s;
                 Plianall(1,3)=(bases.C)/s;
                 Plianall(1,4)=(bases.T)/s;
                %%%%%%%%2-mer frequency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                for j=1:2:s-1  
                   for k=1:16 
                      if ff(1,j:j+1)==str2{1,k}  
                           Plianall(1,k+4)=Plianall(1,k+4)+1;
                      end
                   end
                end  
                
                Plianall(1,5:20)=Plianall(1,5:20)/(s/2);
                
    
             %%%%%%%3-mer frequency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                for j=1:3:s-2  %%%
                   for k=1:64 
                      if ff(1,j:j+2)==str3{1,k}  
                           Plianall(1,k+20)=Plianall(1,k+20)+1;
                      end
                   end
                end  %%%
                
                 Plianall(1,21:84)=Plianall(1,21:84)/(s/3);
                
                %%%%%%%%4-mer frequency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                for j=1:4:s-3  
                   for k=1:256
                      if ff(1,j:j+3)==str4{1,k}  %%%%%%%%%%%%%%
                           Plianall(1,k+84)=Plianall(1,k+84)+1;  %%%%
                      end
                   end
                end  %%%
                   Plianall(1,85:340)=Plianall(1,85:340)/(s/4);
                   
                   
                 %%%%%%%%5-mer frequency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                for j=1:5:s-4  %%%%
                   for k=1:1024
                      if ff(1,j:j+4)==str5{1,k}  %%%%%%%%%%%%%%
                           Plianall(1,k+340)=Plianall(1,k+340)+1;  %%%%
                      end
                   end
                end  %%%
           
                
                  Plianall(1,341:1364)=Plianall(1,341:1364)/(s/5);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 zhuan=Plianall(1,:);  %%%  1X1364
 
 
 for bao=1:1364
     fprintf(fidxie,'%19.18f ',zhuan(1,bao)); 
 end 

fprintf(fidxie,'\n'); 
    
end 
staxie=fclose(fidxie);  %%% 


%%%%%%%%part 2 : train the Fisher Discriminant Model based on 120000 pairs of randomly selected vectors from positive and negative sets



%%%following sentences are reading data from piRNA and non-piRNA vectors£¨random select£¨get mean£¨covariance matrix sw
xun=120000;  %%% training set consists of 120000 positives£¨120000 negatives; 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(1) calculating mean Pa, Na for positive and negative samples %%%%%%%%%%%%%%%%%%%
fidmoreP = fopen('vectorpositive173090.txt','r');  

Pa=zeros(1,1364); 
Psui=randperm(173090); %%%random the all 173090 vectors of positive sample
pding=Psui(1,1:xun);  %%%get the first 120000 vectors from randomized 173090
Pding=sort(pding); %%%sorting them 

dd=1;

for kk=1:173090    
      f=fgetl(fidmoreP);   
   if kk==Pding(1,dd) %%%%
      yun=sscanf(f,'%f'); 
      Pa=Pa+yun';         %%%adding
      dd=dd+1;      
      if dd==xun+1  %%%see if we have enough vectors
         break 
      end
   end 
end
Pa=Pa/xun;  %%% get mean for positive set
stamoreP=fclose(fidmoreP);



fidmoreN = fopen('vectornegative193321.txt','r');  
Na=zeros(1,1364); %%
Nsui=randperm(193321); 
nding=Nsui(1,1:xun); 
Nding=sort(nding); 

dd=1;

for kk=1:193321    
      f=fgetl(fidmoreN);       
      
   if kk==Nding(1,dd) 
      yun=sscanf(f,'%f');
      Na=Na+yun';
      dd=dd+1;
      if dd==xun+1 
         break 
      end
   end 
   
end
Na=Na/xun;  
stamoreN=fclose(fidmoreN);
%%%%%%%%%%%%%%%%%%%%%%% £®1£©end of calculating means for  positive and negative sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%£®2£©covariance matrix sw %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fidmoreP = fopen('vectorpositive173090.txt','r');  
dd=1;
Psw=zeros(1364,1364); %%%because sw=Psw+Nsw;so, first we get Psw for positive samples

for kk=1:173090    
      f=fgetl(fidmoreP);                
   if kk==Pding(1,dd) 
      yun=sscanf(f,'%f');
      Psw=Psw+(yun'-Pa)'*(yun'-Pa);%%%calculating Psw,note that Pa is the mean which has been got above
      dd=dd+1;
      if dd==xun+1  
         break 
      end
   end 
end
stamoreP=fclose(fidmoreP);



fidmoreN = fopen('vectornegative193321.txt','r');  %%%
dd=1;
Nsw=zeros(1364,1364);

for kk=1:193321  
      f=fgetl(fidmoreN);  
      
   if kk==Nding(1,dd) 
      yun=sscanf(f,'%f');
      Nsw=Nsw+(yun'-Na)'*(yun'-Na);%%%%%%calculating Nsw,note that Na is the mean which has been got above
      dd=dd+1;
      if dd==xun+1  
         break 
      end
   end   
end
stamoreN=fclose(fidmoreN);
sw=Psw+Nsw; 
%%%%%%%£®2) the end of calculating sw %%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%£®3£©calculating fisher discriminant vector fN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fN1=inv(sw)*(Pa-Na)'; 
fN=fN1/norm(fN1); %%%get fisher discriminant vector fN,which is the w in my article
%%%%%%%%%%%%%%%%(3)the end of calculating fN %%%%%%%%%%%%%%%%%%%%%s%%%%%%%%%


 
%%%%%%%%%%%%%%%£®4£©calcualting fN*X for 120000 pairs of positive and negative sample%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx3=zeros(1,xun);  %%% ready to record the fN*X for positive samples
yy3=zeros(1,xun);  %%%ready to record the fN*X for negative samples

fidmoreP = fopen('vectorpositive173090.txt','r');  
dd=1;
for kk=1:173090   
      f=fgetl(fidmoreP);                
   if kk==Pding(1,dd) 
      yun=sscanf(f,'%f');
      xx3(1,dd)=yun'*fN;  %%% record the fN*X for positive samples
      dd=dd+1;
      if dd==xun+1 
         break 
      end
   end 
end
stamoreP=fclose(fidmoreP);


fidmoreN = fopen('vectornegative193321.txt','r');  
dd=1;
for kk=1:193321   
      f=fgetl(fidmoreN);      
   if kk==Nding(1,dd)       
      yun=sscanf(f,'%f');
      yy3(1,dd)=yun'*fN;              %%%record the fN*X for negative samples
      dd=dd+1;
      if dd==xun+1 
         break 
      end
   end    
end
stamoreN=fclose(fidmoreN);

%%%%%%%%%%%%%%%%%%%%the end of calcualting fN*X for 120000 pairs of positive and negative sample%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%part 3 : predict%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xxx3=zeros(1,jzhong); %%% get ready to record the fN*X for all sequences to be predicted  
zzz3=zeros(1,jzhong);%%%


fidceP = fopen('vectorprediction.txt','r');  
for kk=1:jzhong    
      f=fgetl(fidceP);        
      yun=sscanf(f,'%f');
      xxx3(1,kk)=yun'*fN;  %%% getting fN*X for all sequences to be predicted  
end
staceP=fclose(fidceP);

%%%%%%%%%%%%%%%%%%%%%%%%
Pxulie=cell(1,jzhong);  %%%% store sequences to be predicted
Pshuo=cell(1,jzhong);   %%%% store the illustration
Plength=zeros(1,jzhong);  %%% store the sequence length

fid1=fopen('prediction.txt','r');  %%%%  
j=0; 

while feof(fid1)==0               % judge if reach the end of the file, otherwise, do the loop.
    
   f=fgetl(fid1);                % read  the title line
   ss=size(f);
   s=ss(1,2);
            if (s>0)&&(f(1)=='>')                
                ff1=fgetl(fid1); 
                j=j+1;
                Pshuo{1,j}=f; 
                Pxulie{1,j}=ff1;
                Plength(1,j)=length(ff1);
            end 
            
end
sta1=fclose(fid1);
%%%%%%%%%%%%%%%%%%%%output the predicted sequences %%%%%%%%%%%%%%%%%%%%%%%%%%%
Nmean=mean(yy3);  %%%% mean of 120000 negative sequences
Nstd=std(yy3);    %%%% std of 120000 negative sequences 
fid8=fopen('predictedpiRNA.txt','w'); %%% for output

% positivece=0; %%%%the number of predicted piRNAs

tiao=2; %%%precision>90%£¨coverage> 60% 
for i=1:jzhong
   
      if   (xxx3(1,i)> Nmean+tiao*( Nstd )) &&(Plength(1,i)>=25) %%% the sequences which fN*X larger than Nmean+2*( Nstd ) and no shorter than 25 are piRNAs
%         positivece=positivece+1;
        
        count1=fprintf(fid8,'%s\n',Pshuo{1,i});         
        count2=fprintf(fid8,'%s\n',Pxulie{1,i}); 
           
        
      end    
end






            













 
 

