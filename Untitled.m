%载入数据
% T4=getgenbank('NC_000866');%T4 phage
% T3=getgenbank('NC_042129');%slur03
% T=getgenbank('AP018932');%Escherichia phage KIT03
% T2=getgenbank('AP018813');%Enterobacteria phage T2
% T6=getgenbank('AP018814');%Enterobacteria phage T6
% T7=getgenbank('NC_003298');%T7 phage
% TK=getgenbank('NC_027995');%Escherichia phage vB_EcoM_ECO1230-10
%%
%显示可翻译基因数目
a=numel(T4.CDS);
% b=numel(T2.CDS);
c=numel(T3.CDS);
% d=numel(T.CDS);
% e=numel(T6.CDS);
% f=numel(T7.CDS);
% g=numel(TK.CDS);
%%
%对空翻译填充
missingpn=find(cellfun(@isempty,{T4.CDS.translation}));
allcds=featureparse(T4,'Feature','CDS','Sequence',true);
missingeqs=cellfun(@nt2aa,{allcds(missingpn).Sequence},'uniform',false);
[T4.CDS(missingpn).translation]=deal(missingeqs{:});
%对空翻译填充
missingpn=find(cellfun(@isempty,{T4.CDS.translation}));
allcds=featureparse(T4,'Feature','CDS','Sequence',true);
missingeqs=cellfun(@nt2aa,{allcds(missingpn).Sequence},'uniform',false);
[T3.CDS(missingpn).translation]=deal(missingeqs{:});
%%
[sc30,globAlign]=nwalign(T4.CDS(1).translation ,T3.CDS(1).translation);
showalignment(globAlign);