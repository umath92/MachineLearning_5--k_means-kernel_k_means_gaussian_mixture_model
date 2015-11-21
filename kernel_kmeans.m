function kernel_kmeans()

circle=load('hw5_circle.mat');
assignin('base', 'hw5_circle',circle.points);
circle=double(circle.points);

kernel_k_2(circle,.055);
end


function kernel_k_2(blob,ss)

k=2;

index=ones(size(blob,1),1);
index(size(blob,1))=2;
index(35)=2;
index(57)=2;
index(93)=2;
index(105)=2;
index_old=inf(size(blob,1),1);

assignin('base', 'index',index);



for i=1:size(blob,1)
    for j=1:size(blob,1)
        K(i,j)=exp(-sum((blob(i,:)-blob(j,:)).^2)/ss);   
    end
end

assignin('base', 'K',K);

%%%%%%%%%%%%%%%%%
d1=Inf(size(blob,1),k);
d2=Inf(size(blob,1),k);
dist2=Inf;

%while(1)
iter=0;
while(iter<15)
    iter=iter+1;
    
    %d1=Inf(size(blob,1),k);
    dist1=0;
    index_old=index;
    
    for row=1:size(blob,1)
        for class=1:k
            intSum=internalSum(index,blob,class,row,K);
            l=lastSum(index,blob,class,K);
            d1(row,class)=K(row,row)-2*(intSum)+l;
        end
    end
    
    dist1=sum(sum(d1));
    if(dist1<dist2)
        [~,I]=min(d1');
    else
        [~,I]=min(d2');
    end
    index=I';
    dist2=dist1;
    d2=d1;
end

% plot stuff
    figure
    a=find(index==1);
    for i=1:size(a,1)
        x(i)=blob(a(i),1);
        y(i)=blob(a(i),2);
    end
    if(size(a,1)>=1)
        scatter(x,y,'g')
    end
    hold on

    b=find(index==2);
    x=[];
    y=[];
    for i=1:size(b,1)
        x(i)=blob(b(i),1);
        y(i)=blob(b(i),2);
    end
    if(size(b,1)>=1)
        scatter(x,y,'y+')
    end
    hold off


end

function [final]=internalSum(index,blob,class,xn,K)

sum=0;
den=0;
for row=1:size(blob,1)
    if (index(row)==class)
      sum=sum+K(row,xn);
      den=den+1;
    end
end

final=sum/den;
end

function [final] = lastSum(index,blob,class,K)
sum=0;
den=0;
for row=1:size(blob,1)
    if(index(row)==class)
        den=den+1;
    end
end

for row=1:size(blob,1)
    if (index(row)==class)
        for row2=row:size(blob,1)
            if (index(row2)==class)
              sum=sum+K(row,row2);
            end
        end
    end
end
final=sum/(den*den);

end
