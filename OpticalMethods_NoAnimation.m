%% Brian Yee
% ID 8365960762
% Wednesday Lab
% Skeleton code supplied by AME 341b

% For animations, uncomment all the plot functions - will increase runtime
% significantly

% AME-341b Special Experiment OPTICAL METHODS (OM). Specify the 'image_set'
% case and run.  Throughout this script there are several variables that
% need to be formulated into proper inputs; these are marked as "###".

clear; close all; clc;
setup_mode = 1;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  specify the working directory where the a/b image pair is stored
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
working_dir = 'C:\Users\4899b\Documents\USC\2021-22\AME 341b\Special Experiment OM\image_test2\';  % keep the "\" at the end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the filenames of the a/b image pair in working_dir.  This will create
% a structure where file names are images(1).name for img_a and,  
% images(2).name for img_b.
images = dir([working_dir,'img*.jpg']);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Load img_a and img_b and Pre-Process image pair:
%  flatten to 1D (grayscale) if the images are multidimensional (e.g., RGB)
%  This block uses the function 'flatten_rgb_image' which we provided; make
%  sure this function is in your Matlab path.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(images)
    % load image data into a temporary variable
    temp = imread([working_dir,images(i).name]);
    % show original image
    if setup_mode == 1
        figure(1);  subplot(2,2,i);
        imshow(temp);  title(images(i).name,'Interpreter','none'); hold on;
    end

	% Assign flattened (grayscale) image data to the structure "IMAGES"; 
    % this data will be used for interrogation.
    images(i).data = flatten_rgb_image(temp,1);
end
% End of image Pre-Processing setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Correlation Parameters for the image pair:
%  image box and search box dimensions (in pixels);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dimensions=[];
for i=1:length(images)
    dimensions(i,1:2)=size(images(i).data);
end

i1 = dimensions(1,:);
i2 = dimensions(2,:);

searchX = floor(i1(1)/4);
searchY = floor(i1(2)/2);
sX = floor(searchX/2);
sY = floor(searchY/2);
imageX = floor(searchX/4);
imageY = floor(searchY/4);
iX = floor(imageX/2);
iY = floor(imageY/2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  For each image box in img_a, calculate the displacement.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nx*ny = total number of displacement calculations (grid points).  This is
% a function of the image size and the Correlation Parameters from above.

nx = ((i1(1)-sX-1-sX-iX)/iX);
ny = ((i1(2)-sY-1-sY-iY)/iY);

% define quiver location arrays
Qx = [];
Qy = [];

%Each outer loops (p & q) creates a single displacement vector. This moves
%the entire "process" around the image

% graph settings
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
figure; set(gcf,'color','w');
%define displacement quiver arrays
quivX = [];
quivY = [];
for p = 1 + sX + iX : iX : i1(1)-sX
    % progress indicator...
    fprintf('%.2f%%...\n',100*(p)/(i1(1)-sX))
    for q = 1 + sY + iY: iY: i1(2)-sY
        fprintf('%.2f%%...\n',100*(p)/(i1(1)-sX) + (100/nx)*((q)/(i1(2)-sY)) )
        % pixel array A
        A = double( images( 1 ).data( ((p-iX):(p+iX)),((q-iY):(q+iY)) ) );  % specify array indices and convert to a double
        % NOTE: imshow does not like doubles, so imshow(uint8(A)) will display A nicely
        % plot img_a.jpg
%         subplot(2,3,1)
%         imshow(images(1).data)
%         r=rectangle('Position',[q-iY p-iX imageY imageX]);
%         drawnow
%         title('$img_a.jpg$')
%         xlabel('x [pixels]');
%         ylabel('y [pixels]');
%         % plot A in img_a.jpg
%         subplot(2,3,2)
%         imshow(uint8(A))
%         title('Image Box in $img_a.jpg$')
%         xlabel('x [pixels]');
%         ylabel('y [pixels]');
%         drawnow
        A_avg = mean(A,'all');
        % Find the displacement of A by correlating this pixel array with all 
        % possible destinations B(K,L) in search box S of img_b.
        %
        % Inner loops (i & j) are the process to find displacement for
        % a given A. 
        C=zeros(1+sX-iX,1);
        for i = -sX+iX: 1 :sX-iX  % x pixel shift within S
            for j = -sY+iY: 1 :sY-iY % y pixel shift within S
                
                % pixel array B      HINT: size(A) = size(B) < size(S)
                B = double( images( 2 ).data( ((p+i-iX):(p+i+iX)),((q+j-iY):(q+j+iY)) ) ); % specify array indices within S and convert to a double
                B_avg = mean(B,'all');
                % plot img_b.jpg
%                 subplot(2,3,4)
%                 imshow(images(2).data)
%                 rectangle('Position',[q-iY p-iX imageY imageX])
%                 rectangle('Position',[q-sY p-sX searchY searchX])
%                 r=rectangle('Position',[q+j-iY p+i-iX imageY imageX]);
%                 title('$img_b.jpg$')
%                 xlabel('x [pixels]');
%                 ylabel('y [pixels]');
%                 % plot B in img_b.jpg
%                 subplot(2,3,5)
%                 imshow(uint8(B))
%                 title('Image Box in $img_b.jpg$')
%                 xlabel('x [pixels]');
%                 ylabel('y [pixels]');
%                 drawnow
                % Calculate the correlation coefficient, C, for this pixel array.
                % Evaluate C at all possible locations (index shifts I,J).
                % The best correlation determines the displacement of A into img_b.
				%  Note: Double sum below effectively implements Double Riemann sum across k and l in lecture
                C(i + 1 + sX-iX , j + 1 + sY-iY ) = sum(sum( (A - A_avg).*(B - B_avg) ))/...
                          sqrt(sum(sum( (A - A_avg).^2 ))*sum(sum( (B - B_avg).^2 )));
            end % j
            % plot correlation coefficient as function of x and y
            % displacements of image box (i, j)
%             subplot(2,3,3)
%             surf(C);
%             xlabel('dx [pixels]');
%             ylabel('dy [pixels]');
%             zlabel('C');
%             title('Correlation Coefficient C in Search Box')
%             drawnow
        end % i
        [test3 test4] = find(ismember(C, max(C(:))));
        Qy = [Qy test3(1)-1-(sX-iX)];
        Qx = [Qx test4(1)-1-(sY-iY)];
        quivX = [quivX q];
        quivY = [quivY p];
        % plot quiver map, updated every time a search box has been
        % searched
%         subplot(2,3,6);
%         imshow(  temp  ); % show image
%         hold on;         % hold figure
%         u = Qx; v = Qy;
%         quiver( quivX,quivY,u,v ,'r','LineWidth',2);
%         xlabel('dx [mm]');
%         ylabel('dy [mm]');
%         title('Quiver Map of Displacement')
%         drawnow
    end % q
end % p

fprintf('Processing complete!\n')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot vectors on shifted image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; set(gcf,'color','w');
imshow(  temp  ); % show image 
hold on;         % hold figure
mmPx = 1/11;         % conversion factor to go from pixels to millimeters
u = Qx*mmPx; v = Qy*mmPx;       %multiply displacement vector arrays by conversion factors
quiver( quivX,quivY,u,v,'r','LineWidth',2);   % plot displacement vectors
% plot final quiver map
title('Quiver Map of Displacement');
xlabel('dx [mm]');
ylabel('dy [mm]');

% Fault line locator
faultY = [];
for i = 2:length(v)
    if u(i) - u(i-1) > 1
        faultY = [faultY quivY(i)];
    end
end

fprintf('Fault line is at y = %.2f +/- %.2f. \n',faultY(2)*mmPx,nx*mmPx)

        
