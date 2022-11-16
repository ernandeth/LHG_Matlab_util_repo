function genimagedelta(P)

spm_defaults;
if nargin<1,
    P = spm_get(2,'*.img',['Please select the images to be processed']);
        Q = 'residuals_delta';
        f = '((i1-i2)*100).*(1/i1)';
        flags = {0,0,spm_type('float'),1};
        spm_imcalc_ui(P,Q,f,flags);
        
end
