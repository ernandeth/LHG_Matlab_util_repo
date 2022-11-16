% (C) Philips Research Europe, Holger Eggers, 2007-05-09
%  File may not be distributed without permission

% Filenames

RawFilename = 'Scan1.raw';
LabFilename = 'Scan1.lab';

% Information to be extracted

ReferenceCardiacPhase = 0;
ReferenceEcho         = 0;

% Initialization

NumberOfProfiles      = 0;
NumberOfSamples       = 0;

clear RawData;

RawFid = fopen(RawFilename);
LabFid = fopen(LabFilename);

fseek (RawFid, 512, 0);

[Label, LabSize] = fread (LabFid,  16, 'uint32');

while LabSize == 16

  DataSize         = Label(1);
  LeadingDummies   = bitand(Label( 2), (2^16-1));
  TrailingDummies  = bitshift (bitand(Label( 2), (2^16-1)*2^16), -16);

  LabelType        = bitshift (bitand(Label( 4), (2^16-1)*2^16), -16);
  ControlType      = bitand(Label( 5), (2^8 -1));

  MeasurementPhase = bitshift (bitand(Label( 5), (2^8 -1)*2^16), -16);
  MeasurementSign  = bitshift (bitand(Label( 5), (2^8 -1)*2^24), -24);

  GainSetting      = bitand(Label( 6), (2^8 -1));
  Mix              = bitshift (bitand(Label( 7), (2^16-1)*2^16), -16);
  Dynamic          = bitand(Label( 8), (2^16-1));
  CardiacPhase     = bitshift (bitand(Label( 8), (2^16-1)*2^16), -16);
  Echo             = bitand(Label( 9), (2^16-1));
  Location         = bitshift (bitand(Label( 9), (2^16-1)*2^16), -16);
  Row              = bitand(Label(10), (2^16-1));
  Extra_Attr       = bitshift (bitand(Label(10), (2^16-1)*2^16), -16);
  Measurement      = bitand(Label(11), (2^16-1));

  E1               = bitshift (bitand(Label(11), (2^16-1)*2^16), -16);
  E2               = bitand(Label(12), (2^16-1));
  E3               = bitshift (bitand(Label(12), (2^16-1)*2^16), -16);

  RfEcho           = bitand(Label(13), (2^16-1));
  GradEcho         = bitshift (bitand(Label(13), (2^16-1)*2^16), -16);

  EncTime          = bitand(Label(14), (2^16-1));
  RandomPhase      = bitshift (bitand(Label(14), (2^16-1)*2^16), -16);

  RRInterval       = bitand(Label(15), (2^16-1));
  RRTopOffset      = bitshift (bitand(Label(15), (2^16-1)*2^16), -16);

  ChannelsActive   = Label(16);


  if LabelType == 32513

    % Standard label

    if ControlType == 0

      % Normal data

      if ((Echo == ReferenceEcho) && (CardiacPhase == ReferenceCardiacPhase))

        [RawProfile, RawSize] = fread (RawFid, (DataSize - LeadingDummies - TrailingDummies) / 2, 'int16');

        fprintf (1, 'Cardiac phase: %3d; Echo: %3d; Profile: %3d\n', CardiacPhase, Echo, E1);

        % Determine number of coils

        NumberOfCoils = 0;

        while ChannelsActive > 0

          if bitand(ChannelsActive, 1)

            NumberOfCoils = NumberOfCoils + 1;

          end

          ChannelsActive = bitshift (ChannelsActive, -1);

        end

        DataSize = (DataSize - LeadingDummies - TrailingDummies) / 4 / NumberOfCoils;

        NumberOfProfiles = max(NumberOfProfiles, E1+1);

        NumberOfSamples = DataSize;

        % Phase correction

        c = exp (- i * pi * (2 * RandomPhase / (2^16-1) + MeasurementPhase / 2));

        for n=1:DataSize

          for m=1:NumberOfCoils

	    if (Measurement > 1)

              RawData(m,E1+1,n) = RawData(m,E1+1,n) + c * (RawProfile(2*(m-1)*DataSize+2*(n-1)+1) + i * RawProfile(2*(m-1)*DataSize+2*(n-1)+2));

            else

              RawData(m,E1+1,n) = c * (RawProfile(2*(m-1)*DataSize+2*(n-1)+1) + i * RawProfile(2*(m-1)*DataSize+2*(n-1)+2));

            end

          end

        end
 
      else

        fseek (RawFid, DataSize - LeadingDummies - TrailingDummies, 0);
  
      end

    else

      fseek (RawFid, DataSize - LeadingDummies - TrailingDummies, 0);
  
    end

  end

  [Label, LabSize] = fread (LabFid,  16, 'uint32');

end

% Write deadface file for testing purposes

Header(1) = 3735943886;
Header(2) = NumberOfSamples;
Header(3) = NumberOfProfiles;
Header(4) = NumberOfCoils;

fid=fopen('Scan1.mdf','wb');
fwrite(fid,Header,'uint32');

for m=1:NumberOfCoils

  for p=1:NumberOfProfiles

    fwrite(fid,real(RawData(m,p,:)),'float32');

  end

  for p=1:NumberOfProfiles

    fwrite(fid,imag(RawData(m,p,:)),'float32');

  end

end

fclose(fid);
