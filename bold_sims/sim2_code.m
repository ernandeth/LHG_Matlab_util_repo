function myDesign = OptimizeDesign;
% outputs a random-ordered list of condition #s that optimizes 3 fMRI considerations
clockStart = clock;					% keeps track of starting time
conditions = [1 2 3 4];				% vector of condition #s (1 for each trial type)
freqConditions = [.35 .35 .15 .15];	% relative freq of conditions - should sum to 1
scanLength = 300;					% length in s of the scan,must be multiple of ISI! 
numCondNoInt = 0;					% 1 = insert a condition of no interest for jitter
ISI = 600;							% inter-stimulus interval
TR = 2;								% TR of your experiment
cbalColinPowerWeights = [1 5 5]; 	% rel. importance of three factors, scale from 0 to 1+
jitter = 0;							% msec value to jitter + or - by
numGenerations = 1000;				% number of generations for genetic algorithm
sizeGenerations = 30;				% number of organisms per generation, must be even number!
lowerLimit = 0.05;					% lower limit of acceptable power frequency, in Hz
maxOrder = 1;						% maximum order of trial sequence counterbalancing to consider
plotFlag = 0;						% flag to subfunctions to plot output; turn on after design is chosen.

% AT THE OUTSET
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
fitnessMatrix = zeros(4,sizeGenerations);
% generate an unsorted list of conditions to use for all random generations
numStim = scanLength / (ISI/1000);
numStimEachCond = numStim * freqConditions;  			% row vector of stim in each cond
unsortedList = [];
for i = 1:size(conditions,2)									% unsorted list of all stims in proper frequencies
   startIndex = size(unsortedList,1)+1;
   unsortedList(startIndex:(startIndex + numStimEachCond(i)-1),1) = conditions(i);	
end

% determine HRF in stim counts, based on ISI
% GAMMA FUNCTION
   n = 4; c = 1;   							% constants for gamma function
   HRF = 1:.001:12;							% 12 seconds, 1 ms resolution
   HRF = HRF.^n .* exp(-HRF .* c);		% gamma function
   HRF = HRF';									% HRF is a column vector
   
   index = 1;
   for i = 1:ISI:size(HRF)					% must resample based on ISI and original sampling rate of HRF.
	  HRF(index) = HRF(i);					% to get HRF in resolution of images
	 index = index + 1;
 	end 
    HRF = HRF(1:index);						% save only the HRF values at points where stimulus presentations occur 
    												% this is the resolution of this program.
   plot(HRF);
   title(['Your gamma-function HRF sampled at ISI = ' num2str(ISI)]);                                      
                                        
                                        
% set up initial set of organisms
disp(['Creating organisms...'])
for z = 1:sizeGenerations % do this for each organism in the FIRST generation
  	% generate random ordering of x conditions
  	stimList = getRandom(unsortedList);
   listMatrix(:,z) = stimList; % a row for each subject
end

% START EVOLUTIONARY ALGORITHM
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
for generation = 1:numGenerations+1
   if generation < numGenerations+1	  
      disp(['Determining fitness of generation ' num2str(generation)]) 
   else 
      disp([num2str(numGenerations) ' generations completed.'])
   end
   % WITHIN A GENERATION
	% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
	for z = 1:sizeGenerations % do this for each organism in the generation
		stimList = listMatrix(:,z);
      
      if (cbalColinPowerWeights(1) > 0)
      	% calculate frequencies each follows each other and deviations from expected random
   		cBal = getCounterBal(stimList, maxOrder,conditions,freqConditions);	 % FUNCTION CALL
      	fitnessMatrix(1,z) = cBal;
      end
      
      if (cbalColinPowerWeights(2) > 0 | cbalColinPowerWeights(3) > 0) 
   		% build predictor set of vectors and convolve with HRF
      	predMatrix = getPredictors(stimList, conditions, HRF);							% FUNCTION CALL
   	end
      
      if (cbalColinPowerWeights(2) > 0)
      	% calculate colinearity of predictors and reduce to a single value
      	% could use cond instead - singular value decomposition.
      	colin = corrcoef(predMatrix);
      	colin = mean(((sum(colin.^2) - 1) / (size(colin,1)-1)) .^0.5); 
      	% average squared correlation minus 1 for self correlation, averaged over non-self corrs 
      	% and converted back to raw (non-squared) correlation
     		fitnessMatrix(2,z) = 1 - colin;
      end  
      
      if (cbalColinPowerWeights(3) > 0)
      	% calculate power spectrum of predictors and collapse to a single value
   		power = getPower(predMatrix, TR, ISI, lowerLimit,plotFlag);
      	fitnessMatrix(3,z) = power;
      end
   end
   % AFTER EACH GENERATION
   % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
   overallFitness = sum(fitnessMatrix,1); % sum the rows in each col
   
   % rank order all organisms and scale according to relative importance of factor
   overallFitnessRank = getFitnessRank(fitnessMatrix, sizeGenerations,cbalColinPowerWeights); 
   
   % breed by rank if not end of last generation
   if generation < numGenerations+1	     
      listMatrix = crossBreed(overallFitnessRank,listMatrix,sizeGenerations,conditions);
   end 
   
  	% determine the most fit organism overall
	mostFit = (overallFitnessRank == max(overallFitnessRank));
	for i = 1:size(mostFit,2)
      if mostFit(i)
         bestFitAcrossGens(generation) = overallFitness(i);
         bestCbal(generation) = fitnessMatrix(1,i);
			bestColin(generation) = fitnessMatrix(2,i);
         bestPower(generation) = fitnessMatrix(3,i);
      	myDesign = listMatrix(:,i);
   	end
   end
   disp(['best of generation ' num2str(generation) ' : Power = ' ...
         num2str(bestPower(generation)) ' Colin. = ' num2str(bestColin(generation)) ...
         ' counterbal. = ' num2str(bestCbal(generation))])
end  % generations loop - finish ev. algorithm

% AFTER ALL GENERATIONS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

save myDesign myDesign			% saves the stim list of chosen design

% plot output
% -------------------- %
plotFlag = 1;
figure
subplot(3,1,1)
plot(bestCbal, 'g')
title('Counterbalancing of best individual in each generation')
subplot(3,1,2)
plot(1-bestColin, 'b')
title('Predictor colinearity of best individual in each generation')
subplot(3,1,3)
plot(bestPower, 'r')
title('Power of best individual in each generation')
% get predictor matrix for chosen one
x = 1:ceil(scanLength / (ISI / 1000));
myPredMatrix = getPredictors(myDesign, conditions,HRF);
figure
for i = 1:size(conditions,2)
   subplot(size(conditions,2),1,i)
   plot(x,myPredMatrix(:,i))
end
subplot(2,1,1);
title('plots of regressors for each condition');

% make power matrix and convert to frequency
[power,myPowerMat] = getPower(myPredMatrix,TR,ISI,lowerLimit,plotFlag); % also plots powerMat
% display summary properties
format short g
format compact
disp('Elapsed Time for this run')
elapsedTime = clock - clockStart
disp('colinearity of predictors')
colin = corrcoef(predMatrix)
disp('fitness and frequency matrix for counterbalancing by order.')
[cBal,freqMatrix] = getCounterBal(stimList, maxOrder,conditions,freqConditions)

for i = 1:size(conditions,2)
   for j = 1:size(conditions,2)
      disp(['condition ' num2str(i) ' followed by condition ' num2str(j) ' = ' num2str(freqMatrix(j,i,1)*100) '%'])
   end
end

return





% ****************************************************************** %

% SUB-FUNCTIONS   
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %   
function stimList = getRandom(stimList)
% gets a random list of conditions in appropriate frequencies
% rand('state',sum(100*clock));								% re-seed the RNG; or not - this makes all orgs the same
randvector = rand(size(stimList,1),1);						% create random vector
stimList(:,2) = randvector;
stimList = sortrows(stimList,2);								% sort the rows by random seed
stimList = stimList(:,1);
return

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ % 
function predMatrix = getPredictors(stimList,conditions, HRF)
for i = 1:size(conditions,2)
   regressor = (stimList == conditions(i));
   regressor = conv(regressor, HRF);
   predMatrix(:,i) = regressor;
end
predMatrix = predMatrix(1:size(stimList),:);      	% eliminate extra values
return

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ % 
function [power,powerMat] = getPower(predMatrix,TR,ISI,lowerLimit,plotFlag)
powerMat = fft(predMatrix);		% powerMat is in units of ISI, not Hz.
powerMat = abs(powerMat);			% make all positive - include reals and imaginaries, all phases
nyquist = 1/(2*TR);										% in Hz
% resample at TR.  can extend to random ISIs...etc..
% compute power and collapse to a single value
newNyquist = ceil(nyquist * size(powerMat,1)* ISI / 1000);	% in predMatrix units (freq. #)
newLowerLimit = ceil(lowerLimit * (size(powerMat,1) * ISI / 1000));
power = mean(powerMat(newLowerLimit:newNyquist,:));	% avg power values in acceptable range for each predictor
power = mean(power);	
% plot output
if plotFlag
   powerMat = fft(predMatrix);
   x =  1:size(powerMat,1);
   x = x / (size(powerMat,1) * ISI / 1000);
   figure
   	for i = 1:size(powerMat,2)
   	subplot(size(powerMat,2),1,i)
      plot(x,powerMat(:,i),'k')
      ymax = abs(max(max(powerMat(2:end,:))));
      axis([0 nyquist 0 ymax]);
      hold on; plot([lowerLimit lowerLimit], [0 max(max(powerMat))],'r');
      plot(x([newLowerLimit:newNyquist]),powerMat([newLowerLimit:newNyquist],i),'b');
   end 
	subplot(2,1,1);
	title('Frequency spectrum of predictors in selected design'); 
end
return
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ % 
function [cBal,freqMatrix] = getCounterBal(stimList, maxOrder,conditions,freqConditions)
cBal = 0;
freqMatrix = zeros(size(conditions,2),size(conditions,2),maxOrder);
for cbalOrder = 1:maxOrder
	for i = 1:size(stimList,1)-cbalOrder	% for every trial besides the last one
		for j = 1:size(conditions,2)			% current element cycles thru conditions
   		for k = 1:size(conditions,2)		% next element cycles thru conditions
      		if stimList(i) == conditions(j) & stimList(i+cbalOrder) == conditions(k) 
            	freqMatrix(k,j,cbalOrder) = freqMatrix(k,j,cbalOrder) + 1;
            	% freqMatrix: condition x condition x order matrix
            	% j = column for each condition; k = index of next trial
         	end
      	end
   	end
	end
   expectedMatrix(:,:,cbalOrder) = freqConditions' * freqConditions;
end
freqMatrix = freqMatrix / size(stimList,1);  % divide count by total stims to get %
diffMatrix = freqMatrix - expectedMatrix;
	MSD = sum(sum(sum(diffMatrix.^2)))/(size(conditions,2)^2*cbalOrder);      % mean squared deviations from expected
		% mean squared deviation per cell of diffMatrix.

cBal = 1 - MSD;
   % absolute max for MSD = 1 * cbalOrder
   % subtracts so that higher numbers are better, and max is 1.
return
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ % 
function overallFitnessRank = getFitnessRank(fitnessMatrix, sizeGenerations,cbalColinPowerWeights) 
   % rank order fitness to decide relative fitness for breeding
   tempFM = fitnessMatrix;				% do this to scale relative importance of values
   for h = 1:3  							% for Cbal, then Colin, then Power
   	for i = sizeGenerations:-1:1
      	myMax = max(tempFM(h,:));
      	for j = 1:sizeGenerations
         	if tempFM(h,j) == myMax & not(tempFM(h,j) == -1) 
               rankMatrix(h,j) = i;
            	tempFM(h,j) = -1;
         	end
      	end
      end
      % scale ranks by relative importance
      rankMatrix(h,:) = rankMatrix(h,:) * cbalColinPowerWeights(h);
   end
   overallFitnessRank = sum(rankMatrix,1); % sum the rows in each col
return
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
% crossbreeding function
function listMatrix = crossBreed(overallFitnessRank,listMatrix,sizeGenerations,conditions)
% SORT OUT WHO STAYS IN THE GENERATION AND CROSSBREED
% disp('   breeding generation...')
numOrgs = size(listMatrix,2);
numStim = size(listMatrix,1);
% keep the best one to ensure hill climbing
bestOne = listMatrix(:,(overallFitnessRank == max(overallFitnessRank)));
bestOne  = bestOne(:,1);					% in case of multiple matches to max

% random adjust to prevent median probs
for z = 1:numOrgs					
   overallFitnessRank(z) = overallFitnessRank(z) + (rand / 100000);
end 

% drop the bad ones (bottom half!)
listMatrix(:,(overallFitnessRank < median(overallFitnessRank))) = [];

% fill out list with clones of good ones
listMatrix = [listMatrix listMatrix];				

if size(listMatrix,2) > sizeGenerations				% Error
   disp('ERROR: LIST MATRIX HAS BECOME TOO LONG!!!')
end   
if mod(size(listMatrix,2),2)								% if odd # of survivors, drop the last one
   listMatrix(:,size(listMatrix,2)) = [];
   disp('ERROR: generation size must be even!!!')
end

% 1% chance of point mutation across whole population
for z = 1:ceil(numStim*numOrgs / 100)
   whichStim = 0; whichOrg = 0; condNum = 0;
   while whichStim < 1 | whichOrg < 1 | condNum < 1 
   	whichStim = round(rand*numStim);
   	whichOrg = round(rand*numOrgs);
      condNum = round(rand*size(conditions,2));
   end   
   listMatrix(whichStim,whichOrg) = condNum;
end

% pair off and swap upper halves
for z = 1:2:numOrgs - 1										% pair off
   crossPoint = 0;											% ensure no 0 crossover point
   while crossPoint < 1
      crossPoint = floor(rand * size(listMatrix,1));	% set crossover point
   end
   % disp(['switching cols ' num2str(z) ' and ' num2str(z+1)])
   a = listMatrix(1:crossPoint,z);
   b = listMatrix(1:crossPoint,z+1); % swap chromosome pieces
   listMatrix(1:crossPoint,z) = b;
   listMatrix(1:crossPoint,z+1) = a;
end

% insert best one randomly somewhere in the list - it survives across gens
whichToMutate = 0;											% ensure no 0 mutate column
while whichToMutate < 1
   whichToMutate = floor(rand * numOrgs);
end
listMatrix(:,whichToMutate) = bestOne; 

return
