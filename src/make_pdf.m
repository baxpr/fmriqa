function make_pdf( ...
	magick_path, ...
	out_path, ...
	seg_file, ...
	rcseg_file, ...
	meanfmri_file, ...
	origfmri_file, ...
	rfmri_file, ...
	mvt_file, ...
	tsnr_file, ...
	rp_file, ...
	med_displc, ...
	glob, FD, DVARS, ...
	project, subject, session, scan ...
	)


%% Figure out screen size so the figure will fit
ss = get(0,'screensize');
ssw = ss(3);
ssh = ss(4);
ratio = 8.5/11;
if ssw/ssh >= ratio
        dh = ssh;
        dw = ssh * ratio;
else
        dw = ssw;
        dh = ssw / ratio;
end


%% Page 1: carpetplots

% Compute them
disp('Carpet plot pre')
[cpt_pre,lbl_pre,mask_pre] = carpet(rcseg_file,origfmri_file);
disp('Carpet plot post')
[cpt_post,lbl_post] = carpet(rcseg_file,rfmri_file);
[cpt_slice,slice_img] = carpet_slice(origfmri_file);

% Number of fmri volumes
nvol = size(cpt_pre,2);

% Create figure
f1 = openfig('pdf_figure.fig','new');
set(f1,'Position',[0 0 dw dh]);
figH = guihandles(f1);

set(figH.scan_info, 'String', sprintf( ...
        '%s %s %s %s', ...
        project, subject, session, scan));
set(figH.date,'String',date);

% Carpet plots (note - depends on label ordering in carpet.m)
cptscale = 2;
lblmap = jet(7);
lblmap = lblmap([1 3 5 7 2 4 6],:);

axes(figH.carpet_pre);
imagesc(cpt_pre,[-cptscale cptscale]);
colormap(figH.carpet_pre,gray);
title('Raw images (6mm spatial smoothing)');
set(gca,'XLim',[0.5 nvol+0.5],'XTick',[],'YTick',[]);

axes(figH.label_pre);
imagesc(lbl_pre);
colormap(figH.label_pre,lblmap);
set(gca,'XTick',[],'YTick',[]);

axes(figH.carpet_post);
imagesc(cpt_post,[-cptscale cptscale]);
colormap(figH.carpet_post,gray);
title('Realigned images (6mm spatial smoothing)');
xlabel('Volume');
set(gca,'XLim',[0.5 nvol+0.5],'YTick',[])

axes(figH.label_post);
imagesc(lbl_post);
colormap(figH.label_post,lblmap);
set(gca,'XTick',[],'YTick',[]);

% Label legend
axes(figH.ax_labelbar)
imagesc(1:7)
colormap(figH.ax_labelbar,lblmap)
set(gca,'XTick',[],'YTick',[])
text(0.6,1.9,'L cortex')
text(1.6,1.9,'R cortex')
text(2.6,1.9,'Subcortical')
text(3.6,1.9,'Cerebellum')
text(4.6,1.9,'Supfl. white')
text(5.6,1.9,'Deep white')
text(6.6,1.9,'CSF')

% Label image
axes(figH.slice_seg);
imagesc(imrotate( mask_pre(:,:,round(size(mask_pre,3)/2)), 90 ));
colormap(figH.slice_seg,[0 0 0; lblmap]);
axis image off

% FMRI image
Yfmri = spm_read_vols(spm_vol(meanfmri_file));
axes(figH.slice_fmri);
imagesc(imrotate( Yfmri(:,:,round(size(Yfmri,3)/2)), 90 ));
colormap(figH.slice_fmri,gray);
axis image off

% Colorbar
axes(figH.ax_colorbar);
img = -cptscale:0.1:cptscale;
imagesc(img,[-cptscale cptscale]);
colormap(figH.ax_colorbar,gray);
set(gca,'XTick',[1 (length(img)+1)/2 length(img)], ...
	'XTickLabel',[-cptscale 0 cptscale], ...
	'YTick',[]);
title('Detrended fMRI signal (%)');

% Head motion
axes(figH.motion);
plot((1:nvol)-0.5,med_displc,'r-','LineWidth',1);
set(gca,'XLim',[0.5 nvol+0.5],'XTick',[]);
set(gca,'YLim',[0 3],'YTick',[0 1 2 3]);
title('Median voxel displacement between vols');
ylabel('mm');

% Print
print(f1,'-dpng','-r600',fullfile(out_path,'output1.png'))



%% Page 2, motion location and slice carpetplot

f2 = openfig('pdf_figure_2.fig','new');
set(f2,'Position',[0 0 dw dh]);
figH = guihandles(f2);

set(figH.scan_info, 'String', sprintf( ...
        '%s %s %s %s', ...
        project, subject, session, scan));
set(figH.date,'String',date);

% Global time series
axes(figH.ax_global)
plot((1:nvol),glob,'LineWidth',1);
set(gca,'XLim',[0.5 nvol+0.5]);
ylabel('Global signal (%)');
xlabel('Volume');

% Motion 
Ymvt = spm_read_vols(spm_vol(mvt_file));
ns = size(Ymvt,3);
slices = round(1:ns/9:ns);
for s = 1:9
	ax = ['slice' num2str(s)];
	axes(figH.(ax))
	imagesc(imrotate(Ymvt(:,:,slices(s)),90),[0 2])
	colormap(figH.(ax),parula);
	axis image off
end

% Motion colorbar
axes(figH.mvt_colorbar)
img = 0:0.01:2;
imagesc(img,[0 2])
colormap(figH.mvt_colorbar,parula)
set(gca,'XTick',[1 length(img)], ...
	'XTickLabel',[0 2], ...
	'YTick',[])
title({'Framewise','displacement','95th %ile'})
xlabel('mm')

% TSNR 
Ytsnr = spm_read_vols(spm_vol(tsnr_file));
ns = size(Ytsnr,3);
slices = round(1:ns/9:ns);
for s = 1:9
	ax = ['tsnr' num2str(s)];
	axes(figH.(ax))
	imagesc(imrotate(Ytsnr(:,:,slices(s)),90),[0 150])
	colormap(figH.(ax),jet);
	axis image off
end

% TSNR colorbar
axes(figH.tsnr_colorbar)
img = 0:1:150;
imagesc(img,[0 150])
colormap(figH.tsnr_colorbar,jet)
set(gca,'XTick',1+(0:50:150), ...
	'XTickLabel',(0:50:150), ...
	'YTick',[])
title({'Temporal SNR','(robust)'})

% Slicewise carpet plot
axes(figH.carpet_slice);
imagesc(flip(cpt_slice));
colormap(figH.carpet_slice,gray);
title('Raw image intensity (mean within slice)');
xlabel('Volume')
set(gca,'XLim',[0.5 nvol+0.5],'YTick',[])

% Slice image
axes(figH.slice_fmri);
imagesc(flip(slice_img));
colormap(figH.slice_fmri,gray);
ylabel('Slice')
set(gca,'XTick',[],'YTick',[1 size(cpt_slice,1)], ...
	'YTickLabel',{num2str(size(cpt_slice,1)) '1'});

% Print
print(f2,'-dpng','-r600',fullfile(out_path,'output2.png'))


%% Page 3: Coregistration
coreg_check( ...
	out_path, ...
    meanfmri_file, ...
    seg_file ...
    );


%% Page 4: Missing coverage
coverage_check( ...
	out_path, ...
    seg_file, ...
    meanfmri_file, ...
    rfmri_file ...
    );


%% Page 5: Motion, FD, DVARS
f3 = openfig('pdf_figure_3.fig','new');
set(f3,'Position',[0 0 dw dh]);
figH = guihandles(f3);

set(figH.scan_info, 'String', sprintf( ...
        '%s %s %s %s', ...
        project, subject, session, scan));
set(figH.date,'String',date);

% Motion params
rp = load(rp_file);

% Translation
axes(figH.translation)
plot((1:nvol),rp(:,1:3),'LineWidth',1);
set(gca,'XLim',[0.5 nvol+0.5],'XTickLabel',[]);
ylabel('Translation (mm)');
legend({'X','Y','Z'},'Location','Best')

% Rotation
axes(figH.rotation)
plot((1:nvol),rp(:,4:6)*180/pi,'LineWidth',1);
set(gca,'XLim',[0.5 nvol+0.5]);
ylabel('Rotation (deg)');
xlabel('Volume');
legend({'X','Y','Z'},'Location','Best')

% FD
axes(figH.FD)
plot((1:nvol),FD,'LineWidth',1);
set(gca,'XLim',[0.5 nvol+0.5]);
ylabel('FD (mm)');
xlabel('Volume');

% DVARS
axes(figH.DVARS)
plot((1:nvol),DVARS,'LineWidth',1);
set(gca,'XLim',[0.5 nvol+0.5]);
ylabel('DVARS (%)');
xlabel('Volume');

% Print
print(f3,'-dpng','-r600',fullfile(out_path,'output5.png'))


%% Combine PNGs to PDF
[status,msg] = system([ ...
	'cd ' out_path ' &&' ...
	magick_path '/convert' ...
	' output1.png output2.png output3.png output4.png output5.png' ...
	' fmriqa_v4.pdf' ...
	]);
if status~=0
    warning('Could not cleanly create PDF file from PNG.');
    disp(msg);
end


