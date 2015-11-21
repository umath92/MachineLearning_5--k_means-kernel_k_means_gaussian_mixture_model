function kmeans()


blob=load('hw5_blob.mat');
assignin('base', 'hw5_blob',blob.points);
blob=double(blob.points);


circle=load('hw5_circle.mat');
assignin('base', 'hw5_circle',circle.points);
circle=double(circle.points);


% Blob data

k_2(blob);
k_3(blob);
k_5(blob);

% Circle data
k_2(circle);
k_3(circle);
k_5(circle);

end

function k_2(blob)

k=2;
u=[Inf,Inf;Inf,Inf];
r=size(blob,1);
%u_new=[-1,1.4;0,0];

% Random Centers %
num=randi(r,1);
qq=blob(num,:);
u_new=[;0,0];
num=randi(r,1);
u_new=[qq;blob(num,:)];
% Random Centers %

index=ones(size(blob,1),1);
while(norm(u_new-u)>0)

    
    u=u_new;
    
    
    % Step 1
    for index_x=1:size(blob,1)
        mn=Inf;
        for index_k=1:k
           current=norm(blob(index_x,:)-u(index_k,:));
           if(current<mn)
               mn=current;
               index(index_x)=index_k;
           end
        end
    end
    assignin('base', 'index',index);

    % Step 2
    sum=[0,0;0,0];
    count=[0;0];
    for i=1:size(index,1)
        if(index(i)==1)
            sum(1,:)=sum(1,:)+blob(i,:);
            count(1)=count(1)+1;
        elseif (index(i)==2)
            sum(2,:)=sum(2,:)+blob(i,:);
            count(2)=count(2)+1;
        end
    end
    u_new=[sum(1,:)/count(1);sum(2,:)/count(2)];
end
    % plot stuff
    figure;
    a=find(index==1);
    for i=1:size(a,1)
        x(i)=blob(a(i),1);
        y(i)=blob(a(i),2);
    end
    if(size(a,1)>=1)
        scatter(x,y,'g');
    end
    hold on;
    
    b=find(index==2);
    x=[];
    y=[];
    for i=1:size(b,1)
        x(i)=blob(b(i),1);
        y(i)=blob(b(i),2);
    end
    if(size(b,1)>=1)
        scatter(x,y,'y+');
    end
    hold off;
    
    % pre plot stuff

end

function k_3(blob)

k=3;
u=[Inf,Inf;Inf,Inf;Inf,Inf];
%u_new=[-1,1.4;1,0.6;-0.5,1];

r=size(blob,1);

% Random Centers %
num=randi(r,1);
qq=blob(num,:);
u_new=[;0,0;0,0];
num=randi(r,1);
qqq=blob(num,:);
num=randi(r,1);
u_new=[qq;qqq;blob(num,:)];


% Random Centers %


index=ones(size(blob,1),1);
while(norm(u_new-u)>0)
    
    
    
    u=u_new;
    
    
    % Step 1
    for index_x=1:size(blob,1)
        mn=Inf;
        for index_k=1:k
           current=norm(blob(index_x,:)-u(index_k,:));
           if(current<mn)
               mn=current;
               index(index_x)=index_k;
           end
        end
    end
    assignin('base', 'index',index);

    % Step 2
    sum=[0,0;0,0;0,0];
    count=[0;0;0];
    for i=1:size(index,1)
        if(index(i)==1)
            sum(1,:)=sum(1,:)+blob(i,:);
            count(1)=count(1)+1;
        elseif (index(i)==2)
            sum(2,:)=sum(2,:)+blob(i,:);
            count(2)=count(2)+1;
        else
            sum(3,:)=sum(3,:)+blob(i,:);
            count(3)=count(3)+1;
        end
    end
    
    if(count(1)>0)
        p=sum(1,:)/(count(1));
    else
        p=u(1,:);
    end
    if(count(2)>0)
        q=sum(2,:)/(count(2));
    else
        q=u(2,:);
    end
    if(count(3)>0)
        r=sum(3,:)/(count(3));
    else
        r=u(3,:);
    end
    u_new=[p;q;r];
    
end

% plot stuff
    figure;
    a=find(index==1);
    for i=1:size(a,1)
        x(i)=blob(a(i),1);
        y(i)=blob(a(i),2);
    end
    if(size(a,1)>=1)
        scatter(x,y,'g');
    end
    hold on;
    
    b=find(index==2);
    x=[];
    y=[];
    for i=1:size(b,1)
        x(i)=blob(b(i),1);
        y(i)=blob(b(i),2);
    end
    if(size(b,1)>=1)
        scatter(x,y,'y+');
    end
    hold on;
    
    c=find(index==3);
    x=[];
    y=[];
    for i=1:size(c,1)
        x(i)=blob(c(i),1);
        y(i)=blob(c(i),2);
    end
    if(size(c,1)>=1)
        scatter(x,y,'c*');
    end
    hold off;
    
    % pre plot stuff
    
end

function k_5(blob)


k=5;
u=[Inf,Inf;Inf,Inf;Inf,Inf;Inf,Inf;Inf,Inf];

r=size(blob,1);

% Random Centers %
num=randi(r,1);
qq=blob(num,:);
num=randi(r,1);
qqq=blob(num,:);
num=randi(r,1);
qqqq=blob(num,:);
num=randi(r,1);
qqqqq=blob(num,:);
num=randi(r,1);
u_new=[qq;qqq;qqqq;qqqqq;blob(num,:)];


% Random Centers %
index=ones(size(blob,1),1);
while(norm(u_new-u)>0)

    
    u=u_new;
    
    
    % Step 1
    for index_x=1:size(blob,1)
        mn=Inf;
        for index_k=1:k
           current=norm(blob(index_x,:)-u(index_k,:));
           if(current<mn)
               mn=current;
               index(index_x)=index_k;
           end
        end
    end
    assignin('base', 'index',index);

    % Step 2
    sum=[0,0;0,0;0,0;0,0;0,0];
    count=[0;0;0;0;0];
    for i=1:size(index,1)
        if(index(i)==1)
            sum(1,:)=sum(1,:)+blob(i,:);
            count(1)=count(1)+1;
        elseif (index(i)==2)
            sum(2,:)=sum(2,:)+blob(i,:);
            count(2)=count(2)+1;
        elseif (index(i)==3)
            sum(3,:)=sum(3,:)+blob(i,:);
            count(3)=count(3)+1;
        elseif (index(i)==4)
            sum(4,:)=sum(4,:)+blob(i,:);
            count(4)=count(4)+1;
        else
            sum(5,:)=sum(5,:)+blob(i,:);
            count(5)=count(5)+1;
        end
    end
    if(count(1)>0)
        p=sum(1,:)/(count(1));
    else
        p=u(1,:);
    end
    if(count(2)>0)
        q=sum(2,:)/(count(2));
    else
        q=u(2,:);
    end
    if(count(3)>0)
        r=sum(3,:)/(count(3));
    else
        r=u(3,:);
    end
    if(count(4)>0)
        s=sum(4,:)/(count(4));
    else
        s=u(4,:);
    end
    if(count(5)>0)
        t=sum(5,:)/(count(5));
    else
        t=u(5,:);
    end
    u_new=[p;q;r;s;t];
end

    % plot stuff
    figure;
    a=find(index==1);
    for i=1:size(a,1)
        x(i)=blob(a(i),1);
        y(i)=blob(a(i),2);
    end
    if(size(a,1)>=1)
        scatter(x,y,'g');
    end
    hold on;
    
    b=find(index==2);
    x=[];
    y=[];
    for i=1:size(b,1)
        x(i)=blob(b(i),1);
        y(i)=blob(b(i),2);
    end
    if(size(b,1)>=1)
        scatter(x,y,'y+');
    end
    hold on;
    
    c=find(index==3);
    x=[];
    y=[];
    for i=1:size(c,1)
        x(i)=blob(c(i),1);
        y(i)=blob(c(i),2);
    end
    if(size(c,1)>=1)
        scatter(x,y,'c*');
    end
    hold on;
    
    c=find(index==4);
    x=[];
    y=[];
    for i=1:size(c,1)
        x(i)=blob(c(i),1);
        y(i)=blob(c(i),2);
    end
    if(size(c,1)>=1)
        scatter(x,y,'r+');
    end
    hold on;
    
    c=find(index==5);
    x=[];
    y=[];
    for i=1:size(c,1)
        x(i)=blob(c(i),1);
        y(i)=blob(c(i),2);
    end
    if(size(c,1)>=1)
        scatter(x,y,'b*');
    end
    hold off;
    % pre plot stuff

end

