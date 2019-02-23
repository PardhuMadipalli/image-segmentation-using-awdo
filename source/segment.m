%Author: Pardhu M
%Contact: pardhu.madipalli@outlook.com


close all;
clear all;
warning('off', 'all')
tic;
format long g;
enmethod='wdo';
result=zeros(7, 3);
f1=1;
for input=[1]            %index of the image, one can write it as input=[1:100] for 100 images with names as im(1).jpg to im(100).jpg 
    disp(num2str(input));
    name=num2str(input);
    input_name=['im(' name ')'];
    I1=((imread(['D:\Major_Project\code\TV1\' input_name '.jpg']))); %corresponds to the image im(1).jpg
    %I=I1;
    I=I1(:,60:(size(I1,2)-10));
    file_bound=['D:\folder_name\' input_name '_orig.png'];
    figure;
    imshow(I); 
    imwrite(I, file_bound, 'PNG')
    dis=I;
    
    %% Parameters for image enhancement
    I=uint8(I);
    orig=I;
    %% Image denoising
    med=Bayes(I); 
  

figure;
imshow(med); 
    %% Image Enhancement Nature inspired
    for lambda = [2]
        for gamma = [1000]
            bord=zeros(size(I));
            fname=[enmethod 'main'];
            fun=str2func(fname);
            EnhanceImage=fun(med, lambda, gamma);
            figure;
            imshow(EnhanceImage); %title('(c)');
            I=EnhanceImage;
            LAM=num2str(lambda);
            GAM=num2str(gamma);
            f=num2str(f1);
            [m n]=size(I);

            %% Image Dilation

            X=EnhanceImage>20;
            figure;imshow(X); % title('(d)');
            %imshow(X,[]);
            h = fspecial('gaussian',[15 15],0.01); 
            X=imfilter(X,h);
            Y=imfill(double(X),'holes');
            Z=bwareaopen(Y,400);
            %subplot(3,3,5);
            figure;
            imshow(Z); % title('(e)');
            SE=strel('disk',2);
            U=imdilate(Z,SE);
            figure;
            
            
            imshow(U); 
            %% Canny Edge Detection
            FinalEdge=edge(U,'canny');
            figure;
            
            
            imshow(FinalEdge);
            [inmi,in]=min(sum(I'));
            xi = in;

            yi = ceil(size(I,1)/2);
            near=zeros(m,n);

            for k=round(xi):1:m            
                l=round(yi);
                    if (FinalEdge(k,l)==1)
                    r2=k;
                    c2=l;
                    near(k,l)=255;
                    break;
                end
            end
            %% LineImage
            MaskImage=zeros(m,n);
            MaskImage(10:m-10,10:n-10)=1;
            LineImage=immultiply(FinalEdge,MaskImage);
            hel=zeros(m,n);
            LowerROI=EnhanceImage(r2-30:r2+30,1:n);
            file_bound=['D:\folder_name\' input_name '_ROI.png'];
            imwrite(LowerROI, file_bound, 'PNG')
            figure;
            
            imshow(LowerROI); 
    %% AWDO code
    im=LowerROI;

    
    [his, bin]=imhist(LowerROI);
    param.popsize = 50;            % population size.
    param.npar = 3;                % Dimension of the problem.
    par=num2str(param.npar);
    param.maxit = 70;              % Maximum number of iterations.
    maxV =     0.3;                % maximum allowed speed.
    dimMin = 1;                     % Lower dimension boundary.
    dimMax = numel(his);			% Upper dimension boundary.

    % Unlike the classical WDO, Adaptive WDO does not need the coefficients to
    % be predetermined. Coefficients; alpha, RT, g, and c, will be selected by
    % the CMAES algorithm:
    rec.arx = rand(4,param.popsize);   %consistent with the CMAES indexing
    %---------------------------------------------------------------

    % Initialize WDO population, position and velocity:
    % Randomize population in the range of [-1, 1]:
    % Please note that WDO/AWDO always works in the range of [-1, 1], but is mapped
    % to the upper/lower bounds of the problem before calling the pressure (cost,fitness) function.

    pos = 2*(rand(param.popsize,param.npar)-0.5);
    % Randomize velocity:
    vel = maxV * 2 * (rand(param.popsize,param.npar)-0.5);  
    %---------------------------------------------------------------

    % Call the pressure (cost,fitness) function
    % Evaluate initial population one member at a time: (simple Sphere Function is utilized here)
    for K=1:param.popsize,
        % map the AWDO 'pos' vector to problem itself using the upper/lower boundaries of each dimension.
        x = round((dimMax - dimMin) * ((pos(K,:)+1)./2) + dimMin);
        size(x);
        % call the pressure (cost,fitness) function:
        pres(K,:) = otsu(his, x);   %here insert your own cost function and make sure the mapping above carried out properly.
    end
    %----------------------------------------------------------------

    % Finding best air parcel in the initial population :
    [globalpres,indx] = max(pres);
    globalpos = pos(indx,:);
    maxpres(1) = max(pres);			% maximum pressure
    %-----------------------------------------------------------------

    % Rank the air parcels in descending order:
    %pres
    [sorted_pres rank_ind] = sort(pres, 'descend');
    %rank_ind
    %pos = pos(rank_ind,:)
    % Do not sort the air parcels! Sorting mixes up CMAES indexing.
    %pos = pos(rank_ind,:);
    keepglob=zeros(param.maxit,1);
    keepglob(1) = globalpres;
    %-----------------------------------------------------------------

    % Start iterations :
    iter = 1;   % iteration counter
    for ij = 2:param.maxit,
            % Update the velocity:
            for i=1:param.popsize
            % choose random dimensions:
            a = randperm(param.npar);        			
            % choose velocity based on random dimension:
                velot(i,:) = vel(i,a);			
                p=find(sorted_pres==pres(i),1); %p is rank here
                vel(i,:) = (1-rec.arx(1,i))*vel(i,:)-(rec.arx(2,i)*pos(i,:))+ ...
                        abs(1-1/p)*((globalpos-pos(i,:)).*rec.arx(3,i))+ ...
                        (rec.arx(4,i)*velot(i,:)/p);
            end

            % Check velocity:
            vel = min(vel, maxV);
            vel = max(vel, -maxV);
            % Update air parcel positions:
            pos = pos + vel;
            pos = min(pos, 1.0);
            pos = max(pos, -1.0); 

            % Call the pressure (cost,fitness) function
            % Evaluate initial population one member at a time: (simple Sphere Function is utilized here)
            for K=1:param.popsize	
                % map the AWDO 'pos' vector to problem itself using the upper/lower boundaries of each dimension.
                x = round((dimMax - dimMin) * ((pos(K,:)+1)./2) + dimMin);
                % call the pressure (cost,fitness) function:
                pres(K,:) = otsu(his, x);   %here insert your own cost function and make sure the mapping above carried out properly.
            end		

            % call CMAES with the coefficients and corresponding pressure, so
            % that CMAES can return the new set of coefficient values for next
            % iteration:
            [rec] = purecmaes_wdo(ij,rec,param.popsize,pres);

            %----------------------------------------------------
            % Finding best particle in population
            [maxpres,indx] = max(pres);
            maxpos = pos(indx,:);          % max location for this iteration
            %----------------------------------------------------
            % Rank the air parcels:
            [sorted_pres rank_ind] = sort(pres, 1, 'descend');
            % Do not sort the air parcels position, velocity and pressure!
            % Instead use "rank_ind" if needed.
            % pos = pos(rank_ind,:);
            % vel = vel(rank_ind,:);
            % pres = sorted_pres;  

            % Updating the global best:
            better = maxpres > globalpres;
            if better
                    globalpres = maxpres;		% global maximum pressure (cost,fitness) value
                    globalpos = maxpos;			% global maximum position vector, note that it is in the range of [-1, 1].
            end
            % Keep a record of the progress:
            keepglob(ij) = globalpres;
    %     	save WDOposition.txt pos -ascii -tabs;
    end
            %Save values to the final file.
            pressure = transpose(keepglob);
            keepsize=size(keepglob);
%             figure;
%             plot(1:param.maxit, keepglob);
%             title('pres vs iter num');
            finalpos = round((dimMax - dimMin) * ((globalpos+1)./2) + dimMin);
            sort(finalpos);
            keepglob(1);

    %         filenamestr = ['WDOpressure.txt'];
    %     	save(filenamestr, 'pressure', '-ascii' , '-tabs');

            % note that the 'globalpos' is the best solution vector found.
            % 'globalpos' is in the range of [-1 1] and needs to be scaled with upper/lower bounds of the specific problem that you are trying to solve.

    a=sort(finalpos);
    levels=a-[ones(1, param.npar)];

    imq=imquantize(im,levels);
    imq=imq-1;
    imf=uint8(imq*(255/param.npar));


    [l, k, m, wid]=findlimits3(imf);
    [row, col]=size(imf);

    %% showing the boundaries on original area
    for i=35:col-35
       dis(r2-30+k(i), i)=255;
       dis(r2-30+l(i), i)=255;
   end
%     for i=35:col-80
%        LowerROI(k(i), i)=255;
%        LowerROI(l(i), i)=255;
%     end

    for i=1:col
       bord(r2-30+k(i), i)=255;
       bord(r2-30+l(i), i)=255;
    end
    figure;
    imshow(dis)
    closed=testing(bord, f1);
    thickness=findlimits4(closed);

    file_bound=['D:\folder_name\' input_name '_bound.png'];
    imwrite(dis, file_bound, 'PNG');
    
    Mean=mean(thickness(35:col-35))*0.062;
    Median=median(thickness)*0.062;
    Median_str=num2str(Median);
    Mean_str=num2str(Mean);
    
    mean_tstr=mean(m(35:col-35))*0.062;
    median_tstr=(median(m)*0.062/f1);
    max_tstr=(max(m)*0.062/f1);
    
    result(input, :)=[input, Mean, mean_tstr]; 
   

         end
    end
 end
%     xlswrite('D:\folder_name\resultsss.xls', result);
%     toc; 
