timeres = 10;

for i = 10:10:90
    sf = [];
    sf{1} = zeros(1000,1);
    sf{1}(round(rand(1,i)*1000),1) = 1;

    [ih,idata,idx]=ideal_Deconv([],sf,8,2,timeres,1);
    title([num2str(i) ' percent'])
end