% We here offer another option for finding landmarks through the watershed
% algorithm.

%load data
datAll = load("../../data/Guoqing2D.mat");
datAll = datAll.dat;
nsubj = 33;
% define the boundry of the data
outbdry = bwboundaries(ones(30));
outbdryIND = sub2ind([30,30], outbdry{1}(:,1), outbdry{1}(:,2));
rad = ones(33,1)*2;
radius = 14; % radius of the searchlight

% subject-specific max number of pixels within each segments
Nmaxs = 10*ones(33,1);
Nmaxs([4,8 19]) = 6;

figure('Position', [0,0, 1300, 1000]);
for i = 1:nsubj
    % perform segmentation
    se = strel('square', 2);
    se2 = strel('square',5);

    dattest = datAll(:,:, i);
    dattest_positive = dattest;
    dattest_positive(dattest <= 0) = 0;

    % circle mask
    mask = 0*dattest == 0;
    [xm, ym]= find(mask);
    distMask = sqrt((xm-mean(xm)).^2+(ym-mean(ym)).^2);
    mask(distMask>radius)=0;

    I =  dattest_positive;
    I(~mask)=0;
    gmag = imgradient(I);
    Ie = imerode(I, se);
    Iobr = imreconstruct(Ie, I);
    Iobrd = imdilate(Iobr,se2);
    Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    fgm = imregionalmax(Iobrcbr);

    imgthreshold = dattest_positive(fgm);
    imgthresholdMean = median(imgthreshold) + 3/(sqrt(2)*erfcinv(3/2)) * median(abs(imgthreshold-median(imgthreshold)));
    imgthd =  fgm & (dattest_positive > imgthresholdMean);

    bwthresh = dattest; 
    bwthresh(~imgthd) = nan;
    bwt = imbinarize(bwthresh, 'global');
    CC = bwconncomp(bwt);

    CCSizes = cellfun(@length, CC.PixelIdxList);
    CCThreshold = find(CCSizes>=4);

    % find the groups with size greater than Nmax.
    localmaxs = [];
    intensMaxs = [];

    for l = 1:length(CCThreshold)
        dat_grp = zeros(size(dattest_positive));
        dat_grp(CC.PixelIdxList{CCThreshold(l)}) = 1;
        bdry = bwboundaries(1-dat_grp, 8);

        bdryLong = cat(1, bdry{:});
        bdryIND =  sub2ind(size(dat_grp), bdryLong(:,1), bdryLong(:,2));
        bdryIndIn = setdiff(bdryIND, outbdryIND);
        dat_grp([CC.PixelIdxList{CCThreshold(l)}; bdryIndIn]) = dattest_positive([CC.PixelIdxList{CCThreshold(l)}; bdryIndIn]);
        [rTemp, cTemp] = find(ismember(dat_grp, max(dat_grp(:))));
        localmaxs = [localmaxs; [rTemp, cTemp]];
        intensMaxs = [intensMaxs; dattest_positive(rTemp, cTemp)];
    end

    [~,iMax] = maxk(intensMaxs, Nmaxs(i));
    nonzeroInd{i,1} = localmaxs(iMax,:);
    Intensities_peaks{i,1} = intensMaxs(iMax);

    subplot(6,6, i);
    imagesc(dattest'); axis xy; colormap jet;  colorbar;
    title(i);
    hold on
    plot(nonzeroInd{i,1}(:,1), nonzeroInd{i,1}(:,2), 'bo');
    hold off
end


