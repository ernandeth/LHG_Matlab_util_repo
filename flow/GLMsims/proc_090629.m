% this script separates the resting and active parts of the time course,
% then 
TR=4; Ttrans = 1.5; pid=1.5; inv_alpha =0.8; Ttag=2;

rest=[	1:11		28:38		53:63		78:88];
active = [	15:25		40:50		65:75		90:100];

rest_tag = raw(rest(1:2:end),:);
rest_con = raw(rest(2:2:end),:);

act_tag = raw(active(1:2:end),:);
act_con = raw(active(2:2:end),:);

subs = [ ...
	'080828ng';
	'080908ss';
	'080911zk';
	'081123ah';
	'081129ay'];

for s=1:size(subs,1)
	cd(['/net/martin/data/hernan/beta2flow/' subs(s,:) '/visualPC'])
	files = dir('vol*.nii');

	root = files(1).name;

	[raw, h] = read_nii_img(root);
	h = nii2avw_hdr(h);



	h.tdim = size(rest_tag,1);
	write_img('rest_tag.img',abs(rest_tag),h);
	write_img('rest_con.img',abs(rest_con),h);
	write_img('act_tag.img',abs(act_tag),h);
	write_img('act_con.img',abs(act_con),h);

	act_ff = calcPerf('act_con', 'act_tag', 'actPerf.img', TR,Ttag,pid,Ttrans, inv_alpha);
	rest_ff = calcPerf('rest_con', 'rest_tag', 'restPerf.img', TR,Ttag,pid,Ttrans, inv_alpha);

	figure;
	lightbox(reshape(act_ff - rest_ff, 64,64,8), [-100 100],3);
end