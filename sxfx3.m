load 'dataphage'
%%
sizeslur03KIT03 = size(dataslur03KIT03);
sizeT4T6 = size(dataT4T6);
sizeT6slur03 = size(dataT6slur03);
%%
%T4-KIT03
maxslur03KIT03 = zeros( 1 , sizeslur03KIT03(1) );
locationslur03KIT03 = zeros( 1 , sizeslur03KIT03(1) );
for i = 1 : sizeslur03KIT03(1)
    maxslur03KIT03(i) = max(dataslur03KIT03(i,:));
    temp1 = dataslur03KIT03(i,:);
    temp2 = find(temp1==max(temp1));
    locationslur03KIT03(i) = temp2(1);
end
%T4-T6
maxT4T6 = zeros( 1 , sizeT4T6(1) );
locationT4T6 = zeros( 1 , sizeT4T6(1) );
for i = 1 : sizeT4T6(1)
    maxT4T6(i) = max(dataT4T6(i,:));
    temp1 = dataT4T6(i,:);
    temp2 = find(temp1==max(temp1));
    locationT4T6(i) = temp2(1);
end
%T6-slur03
maxT6slur03 = zeros( 1 , sizeT6slur03(1) );
locationT6slur03 = zeros( 1 , sizeT6slur03(1) );
for i = 1 : sizeT6slur03(1)
    maxT6slur03(i) = max(dataT6slur03(i,:));
    temp1 = dataT6slur03(i,:);
    temp2 = find(temp1==max(temp1));
    locationT6slur03(i) = temp2(1);
end