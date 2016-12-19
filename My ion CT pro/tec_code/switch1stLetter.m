function siteNew=switch1stLetter(site);

% the first letter of a sitename cannot be a number
% this changes the number to a character

siteNew=site;
if ~isstr(site)
	disp('Input is not a string, ending switch1stLetter !')
	return
elseif length(site)~=4
        disp('Input is a string, but to long for a sitename, ending switch1stLetter !')
	return
end
firstL=site(1);
NumberStrings=['1' '2' '3' '4' '5' '6' '7' '8' '9' '0'];
LetterStrings=['l' 'z' 'e' 'h' 's' 'g' 'y' 'q' 'p' 'o'];

location=find(NumberStrings==firstL);

if ~isempty(location)
	siteNew(1)=LetterStrings(location);
end
