function MRSS = computeRSS(volumename, R1name, R2name)
% function MRSS = computeRSS(volumename, [R1name, R2name])

if (nargin == 1)
	R1name = 'varBefore';
	R2name = 'varAfter';
end
R1 = read_nii_img(R1name);
R2 = read_nii_img(R2name);

Rdiff = (R1-R2)*100./R1;

mask = readvol(volumename);
tmask = zeros(size(mask));
tmask(mask > 2500) = 1;

tmask (isnan(mask)) = 0;
Rdiff(isnan(Rdiff)) =0;
lightbox(tmask);

MRSS = Rdiff(:).*tmask(:);
N = length(find(MRSS));
MRSS = sum(MRSS)/N;

return
