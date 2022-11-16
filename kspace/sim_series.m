DURATION=312;
TR=1;
hdrname= 'sim0001.hdr';
imgname = 'sim0001.img';

noise_amp=0.1;
resp_amp=0;
card_amp=0;
bold_amp=0;

h =read_hdr(hdrname);
img=read_img_data(h,imgname);

mask=read_img2(h,imgname);
mask = zeros(size(mask));
mask(20:30,20:30, h.zdim/2) = 1;
mask = reshape(mask,size(img));

load phys_phases

resp = squeeze(resp_phases(1,1,:));
resp = resp-mean(resp);
resp = resp/max(resp);
resp = resp*resp_amp;

card = squeeze(card_phases(1,1,:));
card = card-mean(card);
card = card/max(card);
card = card * card_amp;

noise = make1overf(DURATION,1/TR);
noise = noise/max(noise);
noise = noise *noise_amp;

bold =resp_gen2(TR,DURATION,15);
bold = bold / max(bold) ;
bold = bold * bold_amp;

for n=1:DURATION/TR
    
    outimg = img + img * (noise(n) + card(n) + resp(n));
    outimg = outimg + (img .* (bold(n)*mask));
    
    str = sprintf('sim_%04d.hdr',n);
    write_hdr(str,h);
    str = sprintf('sim_%04d.img',n);
    write_img_data(str,outimg,h);
    
end
    
