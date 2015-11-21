function GMM()

blob=load('hw5_blob.mat');
assignin('base', 'hw5_blob',blob.points);
blob=double(blob.points);

gmm(blob);
gmm2(blob);
end

function gmm(circle)


%%%%%%%%% Filling w(i,j)-- E step --- %%%%%%%%%%%%%%%%%%
l=1;
colorstring = 'kbgry';
figure
while(l<=5)
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Set the classes randomly %%%%%%%%%%%%%%%%%
    r=size(circle,1);
    keySet=zeros(r,1);
    valueSet=zeros(r,1);
    mapObj = containers.Map(keySet,valueSet);
    q=0;
    index=zeros(size(circle,1),1);
    while(q<=(r/3))
        num=randi(r,1);
        if isKey(mapObj,num)~=1
            mapObj(num) = 1;
            index(num)=1; 
            q=q+1;
        end
    end
    flag=0;
    for i=1:r
        if isKey(mapObj,i)~=1
            if(mod(flag,2)==0)
                index(i)=2;
            else
                index(i)=3;
            end
        end
        flag=flag+1;
    end
    assignin('base', 'index',index);
    %%%%%%%% Done setting the class randomly %%%%%%%%%%%%%%%%%

    %%% Initializing u and co-variance and phi %%%%%%%%%%
    x1=[];
    x2=[];
    x3=[];

    for i=1:r
        if(index(i)==1)
            x1=vertcat(x1,circle(i,:));
        elseif index(i)==2
            x2=vertcat(x2,circle(i,:));
        else
            x3=vertcat(x3,circle(i,:));
        end
    end

    u1=mean(x1);
    u2=mean(x2);
    u3=mean(x3);


    sigma1=cov(x1);
    sigma2=cov(x2);
    sigma3=cov(x3);

    phi1=sum(index==1)/size(index,1);
    phi2=sum(index==2)/size(index,1);
    phi3=sum(index==3)/size(index,1);
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
    
    
    
    
    count=0;
    while(count<=90)
        count=count+1;
        w=zeros(r,3);
        const=power((2*pi),1.5);


        for i=1:r
            for j=1:3
                if j==1
                    w(i,j)=phi1*exp((-1/2)*(circle(i,:)-u1(:,:))*inv(sigma1(:,:))*transpose(circle(i,:)-u1(:,:)))*(1/(const*power(det(sigma1(:,:)),0.5)));
                    summ(i)=w(i,j);
                elseif j==2
                    w(i,j)=phi2*exp((-1/2)*(circle(i,:)-u2(:,:))*inv(sigma2(:,:))*transpose(circle(i,:)-u2(:,:)))*(1/(const*power(det(sigma2(:,:)),0.5)));
                    summ(i)=summ(i)+w(i,j);
                else
                    w(i,j)=phi3*exp((-1/2)*(circle(i,:)-u3(:,:))*inv(sigma3(:,:))*transpose(circle(i,:)-u3(:,:)))*(1/(const*power(det(sigma3(:,:)),0.5)));
                    summ(i)=summ(i)+w(i,j);
                end
            end
        end

        for i=1:r
            for j=1:3
                w(i,j)=w(i,j)/summ(i);
            end
        end


        %%%%%%%%% Estimating parameters -- M step --- %%%%%%%%%%%%%%%%%%
        
        qq=sum(w(:,1))+sum(w(:,2))+sum(w(:,3));
        phi1=sum(w(:,1))/qq;
        phi2=sum(w(:,2))/qq;
        phi3=sum(w(:,3))/qq;

        s1=0;
        s2=0;
        s3=0;
        t1=0;
        t2=0;
        t3=0;

        for i=1:r
            for j=1:3
                if j==1
                    s1=s1+(w(i,j)*circle(i,:));
                    t1=t1+w(i,j);
                end
                if j==2
                    s2=s2+(w(i,j)*circle(i,:));
                    t2=t2+w(i,j);
                end
                if j==3
                    s3=s3+(w(i,j)*circle(i,:));
                    t3=t3+w(i,j);
                end 
            end
        end

        ss1=0;
        ss2=0;
        ss3=0;
        tt1=0;
        tt2=0;
        tt3=0;

        for i=1:r
            for j=1:3
                if j==1
                    ss1=ss1+(w(i,j)*transpose(circle(i,:)-u1)*(circle(i,:)-u1));
                    tt1=tt1+w(i,j);
                end
                if j==2
                    ss2=ss2+(w(i,j)*transpose(circle(i,:)-u2)*(circle(i,:)-u2));
                    tt2=tt2+w(i,j);
                end
                if j==3
                    ss3=ss3+(w(i,j)*transpose(circle(i,:)-u3)*(circle(i,:)-u3));
                    tt3=tt3+w(i,j);
                end
            end
        end

        sigma1=ss1/tt1;
        sigma2=ss2/tt2;
        sigma3=ss3/tt3;

        u1=s1/t1;
        u2=s2/t2;
        u3=s3/t3;

        %%% Compute Log-Likelihood %%%%%
        s=0;
        for j=1:3
            for i=1:r
                if j==1
                    s=s+w(i,j)*(log(phi1)+log(exp((-1/2)*(circle(i,:)-u1(:,:))*inv(sigma1(:,:))*transpose(circle(i,:)-u1(:,:)))*(1/(const*power(det(sigma1(:,:)),0.5)))));
                end
                if j==2
                    s=s+w(i,j)*(log(phi2)+log(exp((-1/2)*(circle(i,:)-u2(:,:))*inv(sigma2(:,:))*transpose(circle(i,:)-u2(:,:)))*(1/(const*power(det(sigma2(:,:)),0.5)))));
                end
                if j==3
                    s=s+w(i,j)*(log(phi3)+log(exp((-1/2)*(circle(i,:)-u3(:,:))*inv(sigma3(:,:))*transpose(circle(i,:)-u3(:,:)))*(1/(const*power(det(sigma3(:,:)),0.5)))));
                end
            end
        end
        sss(count)=s;
        assignin('base', 'w',w);

        %%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%

        for i=1:r
            m=-inf;
            for j=1:3
                if(w(i,j)>m)
                    index(i)=j;
                    m=w(i,j);
                end
            end
        end
    end
    %sss
    plot(sss,'Color', colorstring(l));
    hold on;

    
    l=l+1;
end
    








end


function gmm2(circle)


%%%%%%%%% Filling w(i,j)-- E step --- %%%%%%%%%%%%%%%%%%
l=1;
colorstring = 'kbgry';
while(l<=1)
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Set the classes randomly %%%%%%%%%%%%%%%%%
    r=size(circle,1);
    keySet=zeros(r,1);
    valueSet=zeros(r,1);
    mapObj = containers.Map(keySet,valueSet);
    q=0;
    index=zeros(size(circle,1),1);
    while(q<=(r/3))
        num=randi(r,1);
        if isKey(mapObj,num)~=1
            mapObj(num) = 1;
            index(num)=1; 
            q=q+1;
        end
    end
    flag=0;
    for i=1:r
        if isKey(mapObj,i)~=1
            if(mod(flag,2)==0)
                index(i)=2;
            else
                index(i)=3;
            end
        end
        flag=flag+1;
    end
    assignin('base', 'index',index);
    %%%%%%%% Done setting the class randomly %%%%%%%%%%%%%%%%%

    %%% Initializing u and co-variance and phi %%%%%%%%%%
    x1=[];
    x2=[];
    x3=[];

    for i=1:r
        if(index(i)==1)
            x1=vertcat(x1,circle(i,:));
        elseif index(i)==2
            x2=vertcat(x2,circle(i,:));
        else
            x3=vertcat(x3,circle(i,:));
        end
    end

    u1=[-0.5,1.4];
    u2=[-0.5,1];
    u3=[1,0.2];


    sigma1=cov(x1);
    sigma2=cov(x2);
    sigma3=cov(x3);

    phi1=sum(index==1)/size(index,1);
    phi2=sum(index==2)/size(index,1);
    phi3=sum(index==3)/size(index,1);
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
    
    
    
    
    count=0;
    while(count<=100)
        count=count+1;
        w=zeros(r,3);
        const=power((2*pi),1.5);


        for i=1:r
            for j=1:3
                if j==1
                    w(i,j)=phi1*exp((-1/2)*(circle(i,:)-u1(:,:))*inv(sigma1(:,:))*transpose(circle(i,:)-u1(:,:)))*(1/(const*power(det(sigma1(:,:)),0.5)));
                    summ(i)=w(i,j);
                elseif j==2
                    w(i,j)=phi2*exp((-1/2)*(circle(i,:)-u2(:,:))*inv(sigma2(:,:))*transpose(circle(i,:)-u2(:,:)))*(1/(const*power(det(sigma2(:,:)),0.5)));
                    summ(i)=summ(i)+w(i,j);
                else
                    w(i,j)=phi3*exp((-1/2)*(circle(i,:)-u3(:,:))*inv(sigma3(:,:))*transpose(circle(i,:)-u3(:,:)))*(1/(const*power(det(sigma3(:,:)),0.5)));
                    summ(i)=summ(i)+w(i,j);
                end
            end
        end

        for i=1:r
            for j=1:3
                w(i,j)=w(i,j)/summ(i);
            end
        end


        %%%%%%%%% Estimating parameters -- M step --- %%%%%%%%%%%%%%%%%%
        
        qq=sum(w(:,1))+sum(w(:,2))+sum(w(:,3));
        phi1=sum(w(:,1))/qq;
        phi2=sum(w(:,2))/qq;
        phi3=sum(w(:,3))/qq;

        s1=0;
        s2=0;
        s3=0;
        t1=0;
        t2=0;
        t3=0;

        for i=1:r
            for j=1:3
                if j==1
                    s1=s1+(w(i,j)*circle(i,:));
                    t1=t1+w(i,j);
                end
                if j==2
                    s2=s2+(w(i,j)*circle(i,:));
                    t2=t2+w(i,j);
                end
                if j==3
                    s3=s3+(w(i,j)*circle(i,:));
                    t3=t3+w(i,j);
                end 
            end
        end

        ss1=0;
        ss2=0;
        ss3=0;
        tt1=0;
        tt2=0;
        tt3=0;

        for i=1:r
            for j=1:3
                if j==1
                    ss1=ss1+(w(i,j)*transpose(circle(i,:)-u1)*(circle(i,:)-u1));
                    tt1=tt1+w(i,j);
                end
                if j==2
                    ss2=ss2+(w(i,j)*transpose(circle(i,:)-u2)*(circle(i,:)-u2));
                    tt2=tt2+w(i,j);
                end
                if j==3
                    ss3=ss3+(w(i,j)*transpose(circle(i,:)-u3)*(circle(i,:)-u3));
                    tt3=tt3+w(i,j);
                end
            end
        end

        sigma1=ss1/tt1;
        sigma2=ss2/tt2;
        sigma3=ss3/tt3;

        u1=s1/t1;
        u2=s2/t2;
        u3=s3/t3;

        %%% Compute Log-Likelihood %%%%%
        s=0;
        for j=1:3
            for i=1:r
                if j==1
                    s=s+w(i,j)*(log(phi1)+log(exp((-1/2)*(circle(i,:)-u1(:,:))*inv(sigma1(:,:))*transpose(circle(i,:)-u1(:,:)))*(1/(const*power(det(sigma1(:,:)),0.5)))));
                end
                if j==2
                    s=s+w(i,j)*(log(phi2)+log(exp((-1/2)*(circle(i,:)-u2(:,:))*inv(sigma2(:,:))*transpose(circle(i,:)-u2(:,:)))*(1/(const*power(det(sigma2(:,:)),0.5)))));
                end
                if j==3
                    s=s+w(i,j)*(log(phi3)+log(exp((-1/2)*(circle(i,:)-u3(:,:))*inv(sigma3(:,:))*transpose(circle(i,:)-u3(:,:)))*(1/(const*power(det(sigma3(:,:)),0.5)))));
                end
            end
        end
        sss(count)=s;
        assignin('base', 'w',w);

        %%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%

        for i=1:r
            m=-inf;
            for j=1:3
                if(w(i,j)>m)
                    index(i)=j;
                    m=w(i,j);
                end
            end
        end
    end
    %plot(sss,'Color', colorstring(l));
    %hold on;

    
    %index
    blob=circle;
    % plot stuff
    figure
    a=find(index==1);
    for i=1:size(a,1)
    x(i)=blob(a(i),1);
    y(i)=blob(a(i),2);
    end
    if(size(a,1)>=1)
    disp('1')
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
    disp('2')
    scatter(x,y,'y+')
    end
    hold on

    c=find(index==3);
    x=[];
    y=[];
    for i=1:size(c,1)
    x(i)=blob(c(i),1);
    y(i)=blob(c(i),2);
    end
    if(size(c,1)>=1)
    disp('3')
    scatter(x,y,'c')
    end
    hold off
    
    disp('mean:')
    u1
    u2
    u3
    disp('\n')
    
    disp('Co-variance matrix:')
    sigma1
    sigma2
    sigma3
    
    l=l+1;
end
    








end