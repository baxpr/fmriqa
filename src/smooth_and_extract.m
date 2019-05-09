function data = smooth_and_extract(Y,mask)

osize = size(Y);
if ~all(size(mask)==osize(1:3))
	mask = reshape(mask,osize(1:3));
end
	
sY = nan(prod(osize(1:3)),osize(4));

for t = 1:osize(4)
	tmp = Y(:,:,:,t);
	tmp(~mask) = nan;
	tmp2 = nan(size(tmp));
	spm_smooth(tmp,tmp2,6);
	tmp2(~mask) = nan;
	sY(:,t) = tmp2(:);
end

sY = reshape(sY,[],osize(4))';
data = sY(:,mask);
data = 100 * data ./ repmat(mean(data),osize(4),1);
data = detrend(data)';

