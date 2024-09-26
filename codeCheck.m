% % Open a dialog box to select an image file
% [filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp;*.tif', 'Image files (*.jpg, *.png, *.bmp, *.tif)'; '*.*', 'All files'}, 'Select an Image File');
% 
% %Check if the user clicked 'Cancel' in the dialog
% if isequal(filename, 0) || isequal(pathname, 0)
%     disp('User canceled the operation.');
% else
%     % Load the selected image
%     imagePath = fullfile(pathname, filename);
%     img = imread(imagePath);
% 
%     % Display the image
%     imshow(img);
% end



global workimage;
global dispfig;
global report;
global HC;
global contour;
global skeleton;
global imgObjs;
global tIA;
global imgHistory;
global susstd;
global sapcount;
global bNoInterface;
% lines 18 to 37 defines the fields of chromosomes 
obj1 = struct('contour'    , 0, ... % the field defines the contour of the object
        'skeleton'   , 0, ... % the field defines the midial line of the object
        'original'   , 0, ... % the field defines the grayscale pixels of the object
        'filled'     , 0, ... % the field defines the binary data of the object
        'branchpts'  , 0, ... % the field defines the branchpoints of the object
        'centromere' , 0, ... % the field defines the contour of the object
        'centromereZ', 0, ... % the field defines the location of the centromeres of the object
        'numofcents' , 0, ... % the field defines the numbers of the centromeres of the object unused
        'organgle'   , 0, ...
        'shaynormed' , 0, ...
        'avgval'     , 0, ...
        'imsize'     , 0, ...
        'stdval'     , 0, ...
        'ancut'      , 0, ...
        'clusters'   , 0, ...
        'garbage'    , 0, ... % this tag is a binary tag 1 means the object does not worth further analysis
        'origin'     , 0, ... % origin ==objidx for original run object else origin = idx of originated from object
        'chained'    , 0, ... % array of chained objects originated fron same origin
        'prechange'  , 0, ... % empty if was never spliced m by n by number of changes 2d array otherwise
        'bDelete'    , 0, ... %object appeared in the original image then it was analysed
        'Drsn'       , []); % this tag is for text data that states why an object has been deleted
%if I need to I can ignore it by this
%flag

tRecursion = struct ('visited'  , [] , ...
    'indx'     , [] , ...
    'indy'     , [] , ...
    'branches' , [] , ...
    'branchnum', [] , ...
    'counter'  , [] );

imgObjs = struct('obj'            , [] , ...
    'allcontours'    , 0  , ...
    'allskeletons'   , 0  , ...
    'alloriginals'   , 0  , ...
    'allcentromeres' , 0  , ...
    'totalnumcents'  , 0  , ...
    'orginized'      , 0  , ...
    'orgallcents'    , 0  , ...
    'contflag'       , 0  , ...
    'skelflag'       , 0  , ...
    'centflag'       , 0  , ...
    'numobjects'     , 0   );

tSingleItem = struct('Name'   ,    []   , ...
                     'ID'     ,    []   , ...
                     'Value'  ,    []   );    

tSingleObject = struct('iImgID'           ,     0     , ...
                       'iObjID'           ,     0     , ...
                       'iPosCentromeres'  ,     []    , ... %x, y, id, of only the centers of the centromeres. 
                       'iNumCentromeres'  ,     0     , ... %the number of lines in iPosCentromeres
                       'iPosSkeleton'     ,     []    , ... %A sorted array of all the points of the skeleton of the object. Sorted by appearance. 
                       'iPosContour'      ,     []    , ...
                       'iPosCentromere'   ,     []);
                   
tMultiImgMap = struct('R'    ,     0    , ...
                      'U'    ,     0    , ...
                      'D'    ,     0    , ...
                      'L'    ,     0    );                   

tIA = struct( 'Folders'               , tSingleItem(1)      , ...              
              'MultiImgs'             , tSingleItem(1)      , ...
              'CrImg'                 , 0                   , ...
              'MultiImgMapFile'       , ''                  , ...
              'MultiImgMap'           , tMultiImgMap(1)     , ...
              'ReportFile'            , ''                  , ...
              'tRec'                  , tRecursion(1)       , ...
              'separated'             , []                  , ...
              'imstruct'              , []                  , ...
              'image'                 , []                  , ...
              'Switch_HighResCent4'   , ''                  , ...
              'iNumObjects'           , 0                   , ...
              'Objects'               , tSingleObject(1)    );


tIA.tRec.branches=[];
tIA.tRec.indx=[];
tIA.tRec.indy=[];
tIA.tRec.counter=1;
tIA.tRec.branchnum=1;
tIA.Switch_HighResCent4 = 'MaxP'; %   'Ziv'  or 'MaxP'
sapcount=1;

IA_Initialization;
init_params; % Ziv's note, my parameters are set in this function
chrom=uint8(img);
NumberedAreas=New_ROI(chrom); % define initial regions of intrest
imgObjs.numobjects = max(max(NumberedAreas)); % initial count of number of objects in the image

%20_10_08_MaxP: Ziv, Please make a vector of numbers for contour, that
%includes only the pixels that are participating in making the contour of
%each object. Same for skeleton, same for branch. It is unacceptable to use
%an image for each one of these features. Please use a more efficient
%database structure.
% concluded

construction=zeros(size(chrom)); %create empty matrices for future analysis
cont=construction;

%20_10_31_MaxP: And this is my contour accumulator. The other accumulators
%are changed by selection. Contour is the only one that aint. 
cont2=construction;

skel=construction;
branches=construction;
OLsus=zeros(1,imgObjs.numobjects);
% figure; imshow(uint8(chrom)); title('after halo removal')
j = 0;
q = 0;
for i=1:imgObjs.numobjects %analysis loop
    NA=NumberedAreas;
    NA(NA~=i)=0;
    NA(NA==i)=1; % work only on certin ROI
    temp= double(NA).*double(chrom); % get only the pixels with data within that ROI
    %tempROI=NumberedAreas(NumberedAreas==i);
    
    %20_09_25_MaxP: 
    %Ziv, This function works on the full image. You need to make it work only on the specific region otherwise it is crazily wasteful in terms of resources. Look how much space you are holding for a single object!
    %Specifically, edgeclean rotates through all objects and works on the
    %full image all the time while there is nothing in the area of the
    %object! 
    %I suggest finding each objects maximum and minimum and giving it a
    %bounding rectangle. The rest of the image should be avoided from
    %passing to edgeclean6. Later. For now we just connect things. 
    
    [thresh,edgeimage2]=edgeclean6(temp);% get binary image of contoured objects
    filled=imfill(edgeimage2,'holes'); % fill said contours
    cbranches=bwmorph(edgeimage2,'branchpoints'); % check if contours cut eachother
    [wbx,~]=find(cbranches);
    if ~isempty(wbx) % if contours cut each other seperate them

        cbranches=eightconnth(cbranches); % enlarge contours crossing area
        filled=filled-cbranches; % removes the contours crossing area from the object matrix
        for s=1:4  % clean the seperated objects
            filled=bwmorph(filled,'spur');
        end
        edges=bwmorph(filled,'remove'); % re apply contours to objects
%         if q==0
%            figure;imshow(edges); title('after branchpoints clean') 
%            q=q+1;
%         end
    else
        edges=edgeimage2;
    end
    labeledcontours=bwlabel(edges); % number the contours in the matrix
    zeroimg = zeros(size(labeledcontours));
    for k = 1:max(max(labeledcontours)) % for each contour save an independent object
        %         if (j == 107)
        %           q = 0;
        %         end
        j = j + 1; %objects counter
        imgObjs.obj(j).contour   = labeledcontours==k;
        %figure;imshow(imgObjs.obj(j).contour);
        imgObjs.obj(j).filled    = imfill(imgObjs.obj(j).contour,'holes');
        imgObjs.obj(j).original  = uint8(chrom).*uint8(imfill(imgObjs.obj(j).contour,'holes'));
        imgObjs.obj(j).skeleton  = bwmorph(imgObjs.obj(j).filled,'thin',inf);
        imgObjs.obj(j).branchpts = bwmorph(imgObjs.obj(j).skeleton,'branchpoints');
        imgObjs.obj(j).prechange=[];
        imgObjs.obj(j).origin=j;
        imgObjs.obj(j).chained=[];
        imgObjs.obj(j).ancut=[];
        imgObjs.obj(j).centromereZ=[];
        imgObjs.obj(j).clusters=[];
        imgObjs.obj(j).DF=0;
        imgObjs.obj(j).Drsn=[];
        cont2 = cont2 + imgObjs.obj(j).contour;
        skel=skel+imgObjs.obj(j).skeleton;
        branches=branches+imgObjs.obj(j).branchpts;
        construction=uint8(construction)+uint8(imgObjs.obj(j).original);
    end
    
    %temp(temp<thresh)=0;
    %filled_image; % raw chromosome
    cont=cont+edges;    

end
imgObjs.numobjects = j;


%MaxP 20_12_07 This piece of code was moved from the loop, to outside the
%loop. 


for iN = 1:imgObjs.numobjects
  imgObjs.obj(iN).bOpenShape = CheckOpenCurve(imgObjs.obj(iN).contour); %enter binary flag for if object's contour is a close shape
        %MaxP 20_10_31: And here I first see that the contour is present.
        %Question is how to define if it is a closed shape. Easy: Ziv says: The contour of a single object is here. So if this
        %single object is anywhere near the edge, cut it out for being an open shape. So add a pixel and see if it is at the edge. 
        %Obviously no shape can really be open, as they are all born out of
        %a binary fill. 
        %Now we know which is which and we do not allow open shapes. 
        %And now test that you do indeed detect an open shape this way. 
        %Max: 20_12_07: I see a mess in many figures: an object which is
        %clearly not closed is showing in the list. An object whose contour
        %is in one place and skeleton in another. A mess. It is weird here
        %that I am using i and j indices interchangeably. This is not
        %consistent. I will first put data in the j list, and only then
        %ask. 
  if (tIA.Parameters(3).Value == '1')
    figure; imshow(imgObjs.obj(iN).contour);
    hold on; imshow(imgObjs.obj(iN).skeleton);    
    close;
  end
        
  if (imgObjs.obj(iN).bOpenShape > 0)
    imgObjs.obj(iN).bDelete=1; 
  else
    imgObjs.obj(iN).bDelete=0;         
  end
end

%20_10_31_MaxP: And here I am referring to the accumulated cont, my cont2,
%and not Ziv's cont, which includes ALL objects. My cont2 includes only
%those not near the edge ("closed shapes"). 
%outim1=imfill(cont,'holes');
outim1=imfill(cont2,'holes');

outim=uint8(outim1).*construction;
% figure;imshow(skel);
% figure;imshow(cont);
%structimage.image=construction;
structimage.image=outim;

%20_10_31_MaxP: And here again, I am using only those contours I allow,
%same as only skels I allow etc. 
%structimage.contour=cont;
structimage.contour=cont2;

structimage.skeleton=skel;
structimage.brunchpoint=branches;
structimage.OLsus=OLsus;

%20_10_31_MaxP: And here again, I am using only those contours I allow,
%same as only skels I allow etc. 
%imgObjs.allcontours = cont;
imgObjs.allcontours = cont2;

imgObjs.allskeletons = skel;
imgObjs.alloriginals = construction;
imgObjs.allcentromeres=[];
imgObjs.imsize=size(cont);

% Assuming imgObjs is already defined in codeCheck.m
% Extract each field and assign to the base workspace


% figure; imshow(uint8(construction)); title('after basic cleaning function')
clear chrom;
clear cont;

%20_10_31_MaxP: Cleanup
clear cont2;

clear skel;
clear construction;

%from here
% 
% 
% 

High_Res_Cents(0) % this function searches for centromeres location and validats those areas
disp('cents 1 end');


highressep_Worker220521( 0); % this function tests if objects can be separated and separate them to improve centromere detection
disp('clean and sep 2 end');
High_Res_Cents(0);
% 
% 
% Check if contours are empty
if isempty(imgObjs.allcontours)
    msgbox('Please run FirstRun first');
else
    openGen;

    % Check if centromere flags or centromere data is empty
    if isempty(imgObjs.centflag) || isempty(imgObjs.allcentromeres)
        msgbox('Please run HighRes_Centromeres first');
    else
        % Toggle centromere display
        if imgObjs.centflag == 0
            % Combine original image with centromeres for display
            original_image = uint8(imgObjs.alloriginals);
            centromeres = imgObjs.allcentromeres;

            % Create an RGB image from the grayscale original image
            dispfig = cat(3, original_image, original_image, original_image);

            % Overlay red dots on centromere locations
            red_layer = dispfig(:, :, 1);
            green_layer = dispfig(:, :, 2);
            blue_layer = dispfig(:, :, 3);

            % Find the centromere locations and mark them in red
            [centY, centX] = find(centromeres > 0);
            for k = 1:length(centX)
                red_layer(centY(k), centX(k)) = 255;
                green_layer(centY(k), centX(k)) = 0;
                blue_layer(centY(k), centX(k)) = 0;
            end

            dispfig(:, :, 1) = red_layer;
            dispfig(:, :, 2) = green_layer;
            dispfig(:, :, 3) = blue_layer;

            imshow(dispfig); % Display the image with red centromeres
            imgObjs.centflag = 1;
        else
            % Display original image without centromeres
            dispfig = uint8(imgObjs.alloriginals);
            imshow(dispfig);
            imgObjs.centflag = 0;
        end
    end

    closeGen;
end



% disp('Starting centromere identification process') % inform the user of the progress
% warning('off','all');
% openGen; % open general data matrices
% [a,b]=size(imgObjs.alloriginals);
% imgObjs.allcentromeres=zeros(a,b);
% imgObjs.totalnumcents=0;
% if sapcount~=1 % check if the identification is being done after the object seperation
%     %create an empty txt file to save the data to
%     fid = fopen('Imagedata.txt','wt+');
%     fclose(fid);
% end
% for objnum=1:imgObjs.numobjects
% 
%     if imgObjs.obj(objnum).bDelete==0
%         ActualCentFunc2(objnum); % run the actual centromere analysis function
%        % ActualCentFunc(objnum); % for each non deleted object find the centromeres, costum
%     end
%     if mod(objnum,10)==0 % user progress indicator
%         str=[num2str(objnum),' done out of ',num2str(imgObjs.numobjects),' objects'];
%         disp(str);
%     end
% end
% if sapcount ~=1
%     %convet txt file to csv
%     CreatCSV
% end
% imgObjs.centflag=0;
% closeGen;
% 
% if isempty(imgObjs.allcontours)
%               f= msgbox('Please run FirstRun first')  
%             else
%                 openGen;
%                 if (isempty(imgObjs.centflag) | isempty(imgObjs.allcentromeres) )
%                     f= msgbox('Please run HighRes_Centromeres first')
%                 elseif imgObjs.centflag==0
%                     dispfig(:,:,1)=uint8(uint8(imgObjs.alloriginals))+uint8(255*imgObjs.allcentromeres);
%                     % imagesc(app.ImageAxes,uint8(dispfig));
%                     imgObjs.centflag=1;
%                 else
%                     dispfig(:,:,1)=uint8(imgObjs.alloriginals);
%                     imagesc(app.ImageAxes,uint8(dispfig));
%                     imgObjs.centflag=0;
%                 end
%                 closeGen;
% end
