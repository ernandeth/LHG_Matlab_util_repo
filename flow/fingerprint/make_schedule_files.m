function make_schedule_files(timing_parms, dirName)
% function make_schedule_files(timing_parms, dirName)

del2           = timing_parms.del2;
del3         = timing_parms.del3;
del1        = timing_parms.del1 ;
isLabel         = timing_parms.labelcontrol ;
order           = timing_parms.order;
t_aq            = timing_parms.t_aq ;
AS_delay        = timing_parms.del3;
doAS            = timing_parms.doArtSup;
RO_type         = timing_parms.RO_type;
RO_time         = timing_parms.RO_time;
label_type      = timing_parms.label_type;

% put the delays in us and multiples of 4 so that the scanner can handle
% it.
del1 = del1*1e6;
del2 = del2*1e6;
del3 = del3*1e6;
del1 = 4*round(del1/4);
del2 = 4*round(del2/4);
del3 = 4*round(del3/4);

BGS0 = zeros(size(del2));

Nframes = length(del3)

figure(1)

subplot(321)
plot(timing_parms.del1); title('Delay 1')
axis tight

subplot(322)
plot(timing_parms.del2); title('Delay 2')
axis tight

subplot(323)
plot(timing_parms.del3); title('Delay 3')
axis tight

subplot(324)
plot(del2); title('Label Duration')
axis tight

subplot(325)
stem(timing_parms.labelcontrol); title('Label/Control')
axis tight

subplot(326)
stem(timing_parms.doArtSup); title('Do Art Sup')
axis tight

%%%%%%%%
% make up for scanner error:
%AS_delay = AS_delay-0.1;
%%%%%%%

!mkdir asl3dflex_mrf/
cd asl3dflex_mrf
mat2txt('tadjusttbl.txt', del1, '%d');
mat2txt('prep1_pldtbl.txt', del2, '%d');
mat2txt('prep2_pldtbl.txt', del3, '%d');

mat2txt('doblksattbl.txt', BGS0, '%d');
mat2txt('prep1_lbltbl.txt', isLabel, '%d');
mat2txt('prep2_lbltbl.txt', doAS, '%d');

mat2txt('label_type.txt', string(label_type), '%s');
mat2txt('RO_type.txt', string(RO_type), '%s');
mat2txt('RO_time.txt', RO_time, '%f');

save order.txt       order -ascii
cd ..

str = ['!rm -r ' dirName];
eval(str)
str = ['!mv asl3dflex_mrf ' dirName]
eval(str)
return


function mat2txt(fname,A,fmt)
% function mat2txt(fname,A)
%
% Part of fmrifrey/mri-devtools software package by David Frey (2023)
%   git@github.com:fmrifrey/mri-devtools.git
%
% Description: function to quickly write a matrix to text file
%
%
% Static input arguments:
%   - fname:
%       - name of text file to write matrix to
%       - string describing name of text-based file (with extension)
%       - no default, argument is required
%   - A:
%       - matrix to write to file
%       - 2D matrix containing float/double data
%       - no default, argument is required
%   - fmt:
%       - string format for writing
%       - char array
%       - default is '%f'

    % Set default fmt
    if nargin<3 || isempty(fmt)
        fmt = '%f';
    end

    % Open file for writing
    fID = fopen(fname,'w');

    % Write matrix as float table
    for row = 1:size(A,1)
        for col = 1:size(A,2)
            if ~isstring(A(row,col)) && iscomplex(any(A(:))) && sign(imag(A(row,col))) >= 0
                fprintf(fID,[fmt,'+',fmt,'i \t'],real(A(row,col)),imag(A(row,col)));
            elseif ~isstring(A(row,col)) && iscomplex(any(A(:))) && sign(imag(A(row,col))) < 0
                fprintf(fID,[fmt,'-',fmt,'i \t'],real(A(row,col)),imag(A(row,col)));
            else
                fprintf(fID,[fmt,' \t'],A(row,col));
            end
        end
        fprintf(fID,'\n');
    end

    % Close the file
    fclose(fID);
    return