clc,clear
%%
%‘ÿ»Î ˝æ›
tic
load phage
% KIT03 = getgenbank('AP018932');%Escherichia phage KIT03
% T2 = getgenbank('AP018813');%Enterobacteria phage T2
% T6 = getgenbank('AP018814');%Enterobacteria phage T6
% T4 = getgenbank('NC_000866');%T4 phage
% slur03 = getgenbank('NC_042129');%slur03
%%
%œ‘ æø…∑≠“Îª˘“Ú ˝ƒø
numT4=numel(T4.CDS);
numT2=numel(T2.CDS);
numT6=numel(T6.CDS);
numKIT03=numel(KIT03.CDS);
numslur03=numel(slur03.CDS);
%%
%∂‘ø’∑≠“ÎÃÓ≥‰T4
missingpn=find(cellfun(@isempty,{T4.CDS.translation}));
allcds=featureparse(T4,'Feature','CDS','Sequence',true);
missingeqs=cellfun(@nt2aa,{allcds(missingpn).Sequence},'uniform',false);
[T4.CDS(missingpn).translation]=deal(missingeqs{:});
%∂‘ø’∑≠“ÎÃÓ≥‰T2
missingpn=find(cellfun(@isempty,{T2.CDS.translation}));
allcds=featureparse(T2,'Feature','CDS','Sequence',true);
missingeqs=cellfun(@nt2aa,{allcds(missingpn).Sequence},'uniform',false);
[T2.CDS(missingpn).translation]=deal(missingeqs{:});
%∂‘ø’∑≠“ÎÃÓ≥‰T6
missingpn=find(cellfun(@isempty,{T6.CDS.translation}));
allcds=featureparse(T6,'Feature','CDS','Sequence',true);
missingeqs=cellfun(@nt2aa,{allcds(missingpn).Sequence},'uniform',false);
[T6.CDS(missingpn).translation]=deal(missingeqs{:});
%∂‘ø’∑≠“ÎÃÓ≥‰KIT03
missingpn=find(cellfun(@isempty,{KIT03.CDS.translation}));
allcds=featureparse(KIT03,'Feature','CDS','Sequence',true);
missingeqs=cellfun(@nt2aa,{allcds(missingpn).Sequence},'uniform',false);
[KIT03.CDS(missingpn).translation]=deal(missingeqs{:});
%∂‘ø’∑≠“ÎÃÓ≥‰slur03
missingpn=find(cellfun(@isempty,{slur03.CDS.translation}));
allcds=featureparse(slur03,'Feature','CDS','Sequence',true);
missingeqs=cellfun(@nt2aa,{allcds(missingpn).Sequence},'uniform',false);
[slur03.CDS(missingpn).translation]=deal(missingeqs{:});
%%

for i = 1 : numT4
    for j = 1 : numT6
        
        [sc30,globAlign]=nwalign(T4.CDS(i).translation ,T6.CDS(j).translation);
        
        a=size(globAlign);
        num=0;
        numgap = 0;
        for k = 1 : a(2)
            if globAlign(2,k) == '|'
                num = num + 1;
            end
        end
        for k = 1 : a(2)
            if globAlign(3,k) == '-'
                numgap = numgap + 1;
            end
        end
        dataT4T6(i,j)=num/a(2);
    end
end
for i = 1 : numT6
    for j = 1 : numslur03
        
        [sc30,globAlign]=nwalign(T6.CDS(i).translation ,slur03.CDS(j).translation);
        
        a=size(globAlign);
        num=0;
        numgap = 0;
        for k = 1 : a(2)
            if globAlign(2,k) == '|'
                num = num + 1;
            end
        end
        for k = 1 : a(2)
            if globAlign(3,k) == '-'
                numgap = numgap + 1;
            end
        end
        dataT6slur03(i,j)=num/a(2);
    end
end
for i = 1 : numslur03
    for j = 1 : numKIT03
        
        [sc30,globAlign]=nwalign(slur03.CDS(i).translation ,KIT03.CDS(j).translation);
        
        a=size(globAlign);
        num=0;
        numgap = 0;
        for k = 1 : a(2)
            if globAlign(2,k) == '|'
                num = num + 1;
            end
        end
        for k = 1 : a(2)
            if globAlign(3,k) == '-'
                numgap = numgap + 1;
            end
        end
        dataslur03KIT03(i,j)=num/a(2); 
    end
end
toc