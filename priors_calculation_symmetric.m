% clear all;
datAll = load("../data/Guoqing2D.mat");
datAll = datAll.dat;

%% find local peaks
localpeaks = [];
figure('Position', [0,0, 1500, 1000]);
for i = 1:33
    datThreshed = adaptthresh(datAll(:,:, i), .5);
    localMax = imregionalmax(datThreshed, 8) .* datAll(:,:,i);
    localpeaks(:,:,i) = localMax;
    subplot(6,6, i);
    localMax(localMax == 0) = NaN;
    imagesc(localMax'); axis xy; colormap jet;  colorbar; 
    title(i);
end
figure;
ThreshDat = 1.5e-3;
for i = 1:33
    datThresh = datAll(:,:,i);
    datThresh(datThresh<ThreshDat) = 0;
    localMax = imregionalmax(datThresh, 8) .* datAll(:,:,i);
    subplot(6,6, i);
    localpeaks(:,:,i) = localMax;    
    localMax(localMax == 0) = NaN;
    imagesc(localMax'); axis xy; colormap jet;  colorbar; 
    title(i);
end
nonzeroInd = cell(33,1);
for i = 1:33
   [row, col] = find(localpeaks(:,:,i));
   nonzeroInd{i,1} = [row, col];
end

% If the above steps does not produce any informative landmarks, 
% run find_landmarks_watershed.m file to find the landmarks by watershed
% algtithm.

%% munipulate the reference map
% close all;
% Reference id 20.
pointsAreCollinear = @(xy) rank(xy(2:end,:) - xy(1,:)) == 1;
datAllfig = datAll;
datAllfig(datAll==0) = nan;
[allX, allY] = find(~isnan(datAll(:,:,1)));
datrange = [-max(datAll(:)) max(datAll(:))];

ref = struct;
ref.id = 20;
mean_coord = 15;
mean_ref_coord = [15 15];
ref.coord = [21, 19; 12, 22; 13, 14];

ref.cornerCoord = [10, 13; 25, 13; 25, 24; 10, 24] - mean_ref_coord;
[ref.mesh.X, ref.mesh.Y] = meshgrid((10:24)-mean_ref_coord(1), (13:25) - mean_ref_coord(2));
ref.n = size(ref.coord, 1);
ref.intensities = datAll(:,:, ref.id);
ref.peakIntensities = ref.intensities(sub2ind([30 30],ref.coord(:,1), ref.coord(:,2)));
ref.localIntens = datAll(10:24, 13:25, ref.id)';
ref.coord = [21, 19; 12, 22; 13, 14] - mean_ref_coord;
refsaveMAT = ref.localIntens';
% ref.localIntens(ref.localIntens<0) = 0;
save("../data/reference_3pks.mat", "refsaveMAT");

%% main steps
% figure; imagesc(ref.localIntens); colormap jet;
rs = 0: .1:1;
trs = cell(33, 6);
ninds = 20;
Imin = ones(33,1);
L_min = NaN(33,1);
Loss_RMSE = NaN(33,1);
r = NaN(33,1);

thresh = -1;
ref.localThresheld = ref.localIntens;
ref.localThresheld(ref.localIntens<thresh) = NaN;
% figure; imagesc(ref.localThresheld); colormap jet; colorbar;

srMin = NaN(33,ninds);
ninds1 = zeros(33,1);
coef = NaN(33,1);
coef_inv = NaN(33,1);
l_cauchy = 3;
% save the results to mat. ( prepare for the Rdata)
trsall = struct;
trsall.A = [];
trsall.T = [];
trsall.c = [];
trsall.scale = [];
trsall.scaling = [];
trsall.rotation = [];

alpha = 6;
constraints_b = 0* ones(1,33);
coefmat = NaN(33, 10);
ssdmat = NaN(33, 10);

for i = 1:nsubj
    i
    n_coord = size(nonzeroInd{i, 1}, 1);
    comb = nchoosek(1:n_coord, ref.n);
    n_comb = size(comb, 1);
    
    trA = NaN(n_comb, 1);
    Reg = NaN(n_comb, 1);
    Loss = NaN(n_comb,1);
    Zts = cell(n_comb,1);
    ts = cell(n_comb, 1);
    
    ItrCMin = NaN(n_comb, 1);
    dattemp = datAll(:,:,i)';
    dattemp(dattemp<0) = 0;
    
    % compute the transformation candidates
    for j = 1: n_comb
        nzInd = nonzeroInd{i, 1} - mean_coord;
        %         trCirc = NaN(ref.n,1);
        Tr_reg = NaN(ref.n,1);
        Loss_j = NaN(ref.n,1);
        Zt_j = cell(ref.n, 1);
        t_j = cell(ref.n, 1);
        
        combs = comb(j, :);
        if pointsAreCollinear(nzInd(combs,:)) ~= 1
            k = convhull(nzInd(combs, :));
            if length(unique(k)) == ref.n
                for m = 1: ref.n
                    kcirc = circshift(k(1:ref.n), m-1);
                    ref.peakIntensities=[1,1,1]';
                    intens = ref.peakIntensities;
                    [d, Zt_j{m,1}, t_j{m,1}] = ABC_Procurstean(ref.coord, nzInd(combs(kcirc), :), [ref.peakIntensities, intens]);
                    t = t_j{m,1};
                    
                    % filter out the extreme scaling cases
                    if t.scaling > 0
                        if abs(log(t.scaling(1)/t.scaling(2)))<1
                            Tr_reg(m) = norm(eye(2) - t.A, 'fro')^2;  %S1
                            Loss_j(m) = d;
                        end
                    end
                end
            end
        end
        [Reg(j), ItrCMin(j)] = min(Tr_reg);
        Loss(j) = Loss_j(ItrCMin(j));
        Zts{j} = Zt_j{ItrCMin(j)};
        ts{j} = t_j{ItrCMin(j)};
    end
    
    coef_rs = nan(length(rs),1);
    coef_rs_inv = nan(length(rs),1);
    ssd_rs = nan(length(rs),1);
    inds_rs = nan(length(rs),1);
    
    % loss the constraints at each rs
    for ind = 1: length(rs)
        Ind_reg = find(Reg < rs(ind));
        if ~isempty(Ind_reg)
            inds_loss_feature = find(Loss(Ind_reg)<= min(Loss(Ind_reg))+alpha);
            Ind_trunc = Ind_reg(inds_loss_feature);
            
            coef_r = nan(length(Ind_trunc),1);
            coef_r_inv = nan(length(Ind_trunc),1);
            
            d_r = coef_r;
            for rinds = 1:length(Ind_trunc)
                
                t_r = ts{Ind_trunc(rinds),1};
                
                ZmeshT = t_r.A * [ref.mesh.X(:), ref.mesh.Y(:)]' + t_r.b' + mean_ref_coord';
                intensities = interp2(dattemp, ZmeshT(1,:)', ZmeshT(2,:)', 'spline', NaN);
                
                warpedCorners = t_r.A * (ref.cornerCoord)' + t_r.b' + mean_coord;
                in_idxs = inpolygon(allX, allY, warpedCorners(1,:), warpedCorners(2,:));
                inverse_intensities = dattemp(in_idxs);
                inverse_coord = [allX(in_idxs) allY(in_idxs)];
                Rinverse_coord = inv(t_r.A)*(inverse_coord'-( t_r.b' + mean_ref_coord'))+mean_ref_coord';
                Ref_inv_intensities = interp2(ref.intensities', Rinverse_coord(1,:)', Rinverse_coord(2,:)', 'spline', NaN);
                
                if ~isempty(find(~isnan(intensities), 1))
                    intensities_reshape = reshape(intensities, size(ref.localIntens));
                    intensities_reshape(ref.localIntens<thresh) = NaN;
                    
                    tbl = table(ref.localThresheld(:)/std(ref.localThresheld(:),'omitnan'),...
                        intensities_reshape(:)/std(ref.localThresheld(:),'omitnan'), ...
                        'VariableNames', {'y', 'x'});
                    md  = fitlm(tbl, 'y ~ x-1');
                    coef_r(rinds) = md.Coefficients.Estimate;
                    % inverse L2 loss
                    tbl_inv = table(inverse_intensities/std(ref.localThresheld(:),'omitnan'),...
                        Ref_inv_intensities/std(ref.localThresheld(:),'omitnan'), ...
                        'VariableNames', {'y', 'x'});
                    md_inv = fitlm(tbl_inv, 'y ~ x-1');
                    
                    coef_r(rinds) = md.Coefficients.Estimate;
                    coef_r_inv(rinds) = md_inv.Coefficients.Estimate;
   
                    if coef_r(rinds) > constraints_b(i)
                        d_r(rinds) = .5*(md.RMSE+md_inv.RMSE)*2;
                    end
                end
            end
            
            [ssd_rs(ind), Ind_dmin] = min(d_r);
            coef_rs(ind) = coef_r(Ind_dmin);
            coef_rs_inv(ind) = coef_r_inv(Ind_dmin);
            
            inds_rs(ind) = Ind_trunc(Ind_dmin);
        end
    end
    
    [Loss_RMSE(i),I_rmin] = min(ssd_rs);
    j_min = inds_rs(I_rmin);
    
    coef(i) = coef_rs(I_rmin);
    coef_inv(i) = coef_rs_inv(I_rmin);
    combs_j = comb(j_min, :);
    k_j = convhull(nzInd(combs_j, : ));
    kcirc_j = circshift(k_j(1:ref.n), ItrCMin(j_min)-1);
    trs{i,1} = nzInd(combs_j(kcirc_j), :) + mean_coord;
    trs{i,3} = Zts{j_min,1}+mean_coord;
    trs{i,4} = ts{j_min,1};
    trs{i,5} = trs{i,4}.A * ref.cornerCoord' + trs{i,4}.b' + mean_coord;
    
    % save all results for prior distributions
    trsall.b(i,:) = trs{i,4}.b'; % translation
    trsall.A(:,:,i) = trs{i,4}.A; % A matrix
    trsall.scale(i,1) = coef(i); % intensity correction
    trsall.scale_inv(i,1) = coef_inv(i);% inverse intensity correction
    trsall.scaling(i,:) = trs{i,4}.scaling'; % scaling 
    trsall.rotation(i,:) = trs{i,4}.rotation; % rotation
end

save('../data/symmetric_priors.mat', "trsall");

%% save the maps for inverse registration

for i=1:33
    Ai = trsall.A(:,:,i);
    bi = trsall.b(i,:);
    warpedCorners = Ai * (ref.cornerCoord)' + bi' + mean_coord;
    in_idxs = inpolygon(allX, allY, warpedCorners(1,:), warpedCorners(2,:));
    inverse_intensities = dattemp(in_idxs);
    inverse_coord = [allX(in_idxs) allY(in_idxs)];
    invImgs{i,1} = [inverse_coord, inverse_intensities];
end
save(['../data/symmetric_IMG_priors.mat'], "invImgs");

%% After the compuation from STAN, draw the posterior means 

post =load("../data/posters_30by30_round_3pks_realdata.mat");
post.A(:,:,20) = eye(2);
post.b(20,:) = 0*ones(1,2);
figure('Position', [0,0, 1000, 1100]);
for i = 1:33
    subplot(6,6, i);
    imagesc(datAll(:,:,i)');
    postBox = post.A(:,:,i) * (ref.cornerCoord)' + post.b(i,:)' + mean_coord;
    title("Subject "  +i);
    colormap jet;
    caxis([min(datAll(:)) max(datAll(:))]);
    axis xy;
    hold on;
    plot(postBox(1,[1, 2, 3, 4,1]), postBox(2, [1,2,3,4,1]), '-k');
    hold off;
    caxis(datrange);
end
colp = get(subplot(6,6,18), 'Position');
colorbar('Position', [colp(1) + colp(3) + 0.02  colp(2)  0.01  colp(4)*2.1])
saveas(gcf,'../../figs_paper/realData_posterior.png');

%% Draw the credible regions
post = load("../data/sym_posters_30by30_round_3pks_realdata.mat");
post.A(:,:,20) = eye(2);
post.b(20,:) = 0*ones(1,2);
ref.cornerCoord = [10, 13; 25, 13; 25, 24; 10, 24] - mean_ref_coord;

datAll1 = datAll;
datAll1(datAll==0) = nan;

figure('Position', [0,0, 1000, 1100]);
for subj =1:33
    subplot(6,6, subj);
    CIs = load("../data/CIs/sym_edges_"+subj+".mat");
    CIs = unique(CIs.edges, "rows");
    postBox = post.A(:,:,subj) * (ref.cornerCoord)' + post.b(subj,:)' + mean_coord;
    CI_in_ind = inpolygon(CIs(:,1), CIs(:,2), postBox(1,:), postBox(2,:));
    
    CI_group_out = CIs(~CI_in_ind,:);
    CI_conv_out = convhull(CI_group_out);
    p_out = polyshape(CI_group_out(CI_conv_out,:));
    CI_group_in = CIs(CI_in_ind,:);
    CI_conv_in = convhull(CI_group_in);
    p_in = polyshape(CI_group_in(CI_conv_in,:));
    p_total = polyshape([CI_group_out(CI_conv_out,:) ; CI_group_in(CI_conv_in,:)]);
    
    imagesc(datAll1(:,:,subj)');colormap jet; axis xy;caxis([min(datAll1(:)) max(datAll1(:))]);
    title("Subject "  + subj);
    hold on;
    %         plot(CIs(:,1), CIs(:,2), '.');
    plot(postBox(1,[1, 2, 3, 4,1]), postBox(2, [1,2,3,4,1]), '-k');
    plot(p_total, "FaceAlpha", .5, 'LineStyle','--');
    hold off;
    caxis(datrange);
end

colp = get(subplot(6,6,18), 'Position');
colorbar('Position', [colp(1) + colp(3) + 0.02  colp(2)  0.01  colp(4)*2.1])
saveas(gcf,'figs/sym_realData_posterior.png');
