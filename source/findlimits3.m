function [l, k, m, len]=findlimits3(imf)
[row, col]=size(imf);
k=row*ones(1,col);
l=ones(1,col);
for j=1:col
    p=0;
    for i=1:row
        if imf(i,j)~=0
            k(j)=i;
            p=1;
            break;
        end
    end
        if (j>1 && p==0)
            k(j)=k(j-1);
        end
           

    
    for i=1:row
        q=0;
        if imf(i,j)==255
            l(j)=i;
            q=1;
            break;
        end
    end
        if j>1 && q==0
            l(j)=l(j-1);
        end    
m(j)=l(j)-k(j);
end
len=max(m);
mini=min(m);
c=find(m<0);
d=find(m>0);
posmean=mean(m(d));
m(c)=posmean;
end
