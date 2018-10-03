function [dcell, BPix] = get_statistics(cormatcell)

%calculates relevant statistics.  Inputs are

%cormatcell:    cell array of correlation matrices from MDCCG.  The first
%               cell is the real correlation matrix

%establish array of correlation matrices
niter=length(cormatcell)-1;

arrcan = get_cormatarray(cormatcell);
dcell.original=arrcan.original;
dcell.mean=arrcan.mean;
arrcan=arrcan.array;

%get cutoff for individual pixels
coff1=assess_pixel_cutoff(arrcan(:,:,2:end));

%this second cutoff is for the inference on the size of the regions.
%if we leave the cutoff as coff1 then no regions will be picked up
%at all in the permuted datasets.  With higher bootstrapped samples, a
%differet cutoff for coff1 and coff2 isn't necessary.  This is simply done
%for illustrative purposes to reduce the computational burden of an example
%run out of the box.
coff2=quantile(arrcan(:,:,2:end), .95, 3);



%utilizes image processing toolbox.  

corcompcell=cell(1,niter+1);
BArea=cell(1,niter+1);

BPix=cell(1,niter+1);
Bsum=cell(1,niter+1);


%loop through all of the correlation matrices and caluclate the statistics

for knd=1:(niter+1)
    if knd==1
        coff=coff1;
    else
        coff=coff2;
    end
    indmat=arrcan(:,:,knd)>coff;

    corcompcell{1, knd}=indmat.*arrcan(:,:,knd) - indmat.*coff;

    BPix{1,knd}=regionprops(corcompcell{1,knd}>0, corcompcell{1,knd}, 'Area', 'PixelIdxList', 'PixelValues', 'FilledImage', 'SubarrayIdx', 'MajorAxisLength', 'MinorAxisLength', 'ConvexHull');

    BArea{1,knd}=cat(1,BPix{1,knd}.Area);


    fields = {'FilledImage', 'SubarrayIdx', 'MajorAxisLength', 'MinorAxisLength', 'ConvexHull'};
    tmp=rmfield(BPix{1, knd}, fields);

    Btmp=cellfun(@sum, struct2cell(tmp));

    Bsum{1,knd}=Btmp(3,:); %PixelIntensity sum

end
sizes=cell2mat(BArea(1,2:(niter+1))');
sums=cell2mat(Bsum(1,2:(niter+1)));





sizecutoff=assess_excursion_cutoff(sizes);
sumscutoff=assess_excursion_cutoff(sums);





Busesize=find(BArea{1,1}>sizecutoff);
Busesum=find(Bsum{1,1}>sumscutoff);

takesize=BPix{1,1}(Busesize);
takesum=BPix{1,1}(Busesum);



compsize=zeros(size(cormatcell{1,1}));
compsum=zeros(size(cormatcell{1,1}));


compsize(cat(1, takesize.PixelIdxList))=1;
compsum(cat(1, takesum.PixelIdxList))=1;
compsizeGrayscale=compsize.*corcompcell{1,1};
compsumGrayscale=compsum.*corcompcell{1,1};



dcell.compsumGrayscale=(abs(compsumGrayscale)>0).*dcell.original;
dcell.compsizeGrayscale=(abs(compsizeGrayscale)>0).*dcell.original;

end


function array_cor = get_cormatarray(cormatcell)
    niter=length(cormatcell)-1;
    array_cor.array = zeros([size(cormatcell{1,1}), niter+1]);
    for knd=1:(niter+1)
             array_cor.array(:,:,knd)=abs(cormatcell{knd});
    end
    array_cor.original = array_cor.array(:,:,1);
    array_cor.mean=mean(array_cor.array(:,:,2:end), 3);
    array_cor.array = array_cor.array-repmat(array_cor.mean, 1, 1, niter+1);
end

