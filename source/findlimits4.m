function m=findlimits4(imf)
[row, col]=size(imf);
countk=1;
k=row*ones(1,col);
l=ones(1,col);
for j=1:col
    p=0;
    for i=1:row-1
        if imf(i,j)==255 && imf(i+1,j)==255
            k(countk)=i;
            p=1;
            break;
        end
    end
    
    if p==1          
        q=0;
        for i=row:-1:k(countk)+6
            if imf(i,j)==255 && imf(i-1,j)==255
                l(countk)=i;
                q=1;
                break;
            end
        end
        countk=countk+1;
  
    
    if p==1 && q==0 && countk>2
        l(countk-1)=l(countk-2);
    end
    
    end
        
end
m=l-k;
len=max(m);
mini=min(m);
c=find(m<0);
d=find(m>0);
posmean=mean(m(d));
m(c)=posmean;
end
