[biofilmradius2new, Rframelist2new, Raverage2new] = biofilm_PDMS(80,'stack19/overnight6_oldChannels_37C_30C_3_190319_TL_','190319_radius/radiusfilm',294,'half',1.4,2.) 
%[biofilmradius9new, Rframelist9new, Raverage9new] = biofilm_PDMS(80,'K2biofilm_090419/090419_bf/K2biofilm_MSggShocks_090419_TL','K2biofilm_090419/090419_radius2/radiusfilm',328,'half', 1.2,3)


function [biofilmradius,Rframelist, Raverage] = biofilm_PDMS(pillarlength,biofilmlocation,biofilmlocationBW,...
    numframes,circletype,radiusparam,radiusparam0)
% pillarLength is the length of the pillar in pixels, I removed this value
% so that it didn't falsely detect it as a newton's ring.
%Biofilmlocation is the location for the images for the biofilm. If you use
%biofilm images cropped in half, make sure the biofilm is centered and that
%the radial growth is facing east. biofilmlocationBW is the location for
%the images of the black and white version. We extract area values from
%these images which then gives us an estimation of the radius of the
%biofilm.
%In circletype you type in either "full" or "half" for the version of
%biofilm that you intend to investigate
%Radiusparam is a scalar that the radii gets mulitplied by, all values
%outside this radii gets removed. If the Newton's rings are stationary
%until the biofilm reaches a certain size, radiusparam0 is used to multiply
%the radii by a larger value and this value is then constant until the area
%of the biofilm catches up, it then resumes to use radiusparam as a scalar.
%If you have stationary newton's rings, radiusparam0 > radiusparam, if they
%are not stationary then radiusparam0 = radiusparam.

biofilmradius = [];
t = 1;
framelist = {};
framelist1 = {};
newtonringsperframe = {};
Raverage = [];
Rframelist = [];


radiilist = [];
framelist = [];

for a = 1 : numframes
    if a < 10
        a_ref = strcat('000',num2str(a));
    elseif a < 100
        a_ref = strcat('00',num2str(a));
    elseif a < 1000
        a_ref = strcat('0',num2str(a));
    else
        a_ref = num2str(a);
    end
    biofilm_location = strcat(biofilmlocation,a_ref,'.tif');
    biofilm_location2 = strcat(biofilmlocationBW,a_ref,'.tif');
    bf = imread(biofilm_location);
    bf2 = imread(biofilm_location2);
   
    polarIm = ImToPolar(double((adapthisteq(im2gray(bf))))/255,0,1,700,800);
    polarIm(1:pillarlength,:) = 0; %Remove pillar size from image.
    mysize = size(polarIm);
    w = mysize(2);
    xsize = 400;
    alpha = 0.4; %Thresshold value
    h = mysize(1);
    %Here we store the pixel averages for every horizontal line
    pavrlist = []; 
    %Here all the pixel averages for the whole frame is stored
    %It is stored for all y-coordinates.
    framepav = {};

    
    for y=1:h
        pv = 0;
        for x=1:xsize %I wanted to keep essential parts of the biofilm.
            if polarIm(y,x) > alpha %If above threshhold, then count as a potential "line".
                pv = 1 + pv;        %Add pixel values for each horizontal line   
            end
        end
        pavrlist = [pavrlist, pv/xsize]; %Divide by number of horizontal pixels for the biofilm to get the average
    end
    framepav = [framepav, pavrlist]; 
    framelist = [framelist; framepav];
    

    frame_a = cell2mat(framelist(t,1)); 
    beta = mean(frame_a)*1.; 
    lines = frame_a>beta; 
    lines_img = zeros(h,w);

    for q=1:h
        if lines(1,q) == 1
            lines_img(q,:) = 1;
        end
    end

    %Area conditions for the BW images. 
    if circletype == "half"
     rB = round(sqrt(2*bwarea(bwareafilt(logical(bf2),1))/pi)); %model the biofilm as a circle and figure out the radii.
     Rplaceholder = rB;
     if t == 1
         rB0 = rB * radiusparam0;
     end
    end
    if circletype == "full"
     rB = round(sqrt(bwarea(bwareafilt(logical(bf2),1))/pi)); %model the biofilm as a circle and figure out the radii.
     Rplaceholder = rB;

     
      if t == 1 
          %Initial condition for the radius parameter.
         rB0 = rB * radiusparam0;
     end
    end
    a
    sizebf2 = size(bf2);
    biofilmradius = [biofilmradius, rB];

    if rB < Rplaceholder*radiusparam
        %Whenever the biofilm reaches the same radii as the initial
        %condition, start using a dynamic radii.
        rB = Rplaceholder * radiusparam;
    end


    polarIm1 = polarIm;
    polarIm1 = imresize(polarIm1,[sizebf2(1),sizebf2(2)]);
    polarIm1(round(rB):end,:) = 0; %Remove all pixels outside rB.
    alpha = 0.3; %Thresshold value
    h = mysize(1);
    %Here we store the pixel averages for every horizontal line
    pavrlist1 = []; 
    %Here all the pixel averages for the whole frame is stored
    %It is stored for all y-coordinates.
    framepav1 = {};
    for y=1:h
        pv = 0;
        for x=1:xsize %I wanted to keep essential parts of the biofilm.
            if polarIm1(y,x) > alpha %If above threshhold, then count as a potential "line".
                pv = 1 + pv;        %Add pixel values for each horizontal line   
            end
        end
        pavrlist1 = [pavrlist1, pv/xsize]; %Divide by number of horizontal pixels for the biofilm to get the average
    end
    
    framepav1 = [framepav1, pavrlist]; 
    framelist1 = [framelist1; framepav1];
    beta = mean(frame_a)*1 ; %This beta works properly for the newton's rings.
    frame_b = cell2mat(framelist1(t,1)); 
    lines1 = frame_b<beta; 
    lines_img1 = zeros(h,w);
    radialdistances = []; %Store the radial distances here.
    for q=1:h
        if lines1(1,q) == 1
            lines_img1(q,:) = 1;
            radialdistances=[radialdistances, h-q+1]; %h-q+1 since the central region starts at the bottom of the polar transformation.
        end
    end

   f = fspecial('gaussian');
   img = adapthisteq(flip(filter2(f,polarIm1),1))>0.8;
   img = flip(img);
   montage({img,polarIm}) %Keep this to make sure that the algorithm works properly each time. This helps with choosing radius parameters.
   sizeimg = size(img);
   h = sizeimg(1);
   %Study at which diameter the rings occur
   ylist = zeros(1,h);
    for i=1:h
        if sum(img(i,:))>30
            ylist(i)=1;
        end
    end

    %Constructive interference condition.
    nlist = 1:16; %Put the number of rings to test here.
    nlist = sqrt((nlist-0.5).*(nlist+0.5).^-1); %Create a comparison list
    diameters = unique(find(ylist)); %Here we find the values of the diameters of the newton rings from the polar transformation.
    sizediam = size(diameters);
    epsilon = 0.01; %tolerance value for the inequality
    newtonrings = {};
    sizenlist = size(nlist);
    lambda = 500; %light wavelength
    Rvar = 0;
    Rnum = 0;
    
    for pa=1:sizediam(2)-1
       for qa=pa+1:sizediam(2)
           
         if pa~=qa
            quot = (diameters(pa)*diameters(pa))/(diameters(qa)*diameters(qa));
            for l=1:sizenlist(2)-1
                for h=l+1:sizenlist(2)
                Nlh = sqrt((l-0.5)/(h-0.5));
                    if abs(quot-Nlh) < epsilon
                       
                       radius = diameters(pa);
                       R = (radius*radius)/(lambda*(l-0.5));

                       %If you want to save the information about the
                       %diameters of the Newton's rings, uncomment these
                       %two lines below. The computation will become much slower.
                       %var = {a,[diameters(pa), diameters(qa)], R };
                       %newtonrings = [newtonrings; var]; 

                       Rvar = Rvar + R;
                       Rnum = Rnum + 1;
                    end
                end
            end
            end
         end
    end
   t = t + 1;
    newtonsringsperframe = [newtonringsperframe; newtonrings];
    Rvar = Rvar/Rnum; %Average value of the radius of curvature of PDMS.
    Raverage = [Raverage, Rvar];
    Rframelist = [Rframelist, a];
    scatter(Rframelist,Raverage,'filled') 
    
end
