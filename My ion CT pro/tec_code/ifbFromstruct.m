function ifb=ifbfromStruct(IECstruct,tgdAll)

% computes the receiver ifb from iec structure

warning off MATLAB:nearlySingularMatrix;

%tgdAll=get_tgd_ionex(ionexfile);

felder=fieldnames(IECstruct);
ifb=[];
for n=1:length(felder)
      name=cell2mat(felder(n));
      prnNo=str2num(name(4:length(name)));
      tgd=tgdAll(find(tgdAll(:,1)==prnNo),2);
      eval(['P=IECstruct.' name ';']);
      rowsP=size(P,1);
      if mod(rowsP,2)==1, P(rowsP,:)=[];end		
      for nn=1:2:size(P,1)-1
      ifbtmp=get_ifb(P(nn:nn+1,13),tgd,P(nn:nn+1,12));
      ifb=[ifb; ifbtmp];
      end
end
ifb=mean(ifb(find(~isnan(ifb))));

