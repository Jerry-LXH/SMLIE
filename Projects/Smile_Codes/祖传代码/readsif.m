%reading Andor EMCCD sif files by Yunxiang Zhang
%email : yxzhang@gmail.com
%compatibility : Andor sif version 4.21~4.30
function [info] = readsif(file)
f=fopen(file,'r'); % read the file
if f < 0
   errordlg('Could not open the file.'); % Unable or damaged 
end
tline=fgetl(f);
if ~isequal(tline,'Andor Technology Multi-Channel File')
   fclose(f);
   errordlg('Not an Andor SIF image file.');
end
skipBytes(f,8);
o=fscanf(f,'%f',6);
info.date=datestr(o(5)/86400 + 719529);
info.temperature=o(6);
skipBytes(f,10);
o=fscanf(f,'%f',5);
info.exposureTime=o(2);
info.cycleTime=o(3);
info.accumulateCycles=o(5);
info.accumulateCycleTime=o(4);
skipBytes(f,2);
o=fscanf(f,'%f',5);
info.stackCycleTime=o(1);
info.pixelReadoutTime=o(2);
info.gainDAC=o(5);
%disp(info.gainDAC);%print em gain
o=fscanf(f,'%f',5);
skipBytes(f,4);
o=fscanf(f,'%f',15);
info.VertShift=o(13);
info.PreAmpGain=o(15);
info.SerialNumber=fscanf(f,'%d',1);
%disp({info.VertShift,info.PreAmpGain,info.SerialNumber});%print serial
o=fscanf(f,'%f',9);
o=fscanf(f,'%d',5);
subver=o(2);
info.SoftwareVersion = [o(1) o(2) o(3) o(4) o(5)];
%disp(info.SoftwareVersion);
%o=fscanf(f,'%d',18);skipBytes(f,1);
fgetl(f);
info.detectorType=deblank(fgetl(f));
%disp(info.detectorType);
info.detectorSize=fscanf(f,'%d',[1 2]);
%disp(info.detectorSize);

info.fileName=readString(f);
csize=fscanf(f,'%d',[1 2]);
chunk=fread(f,[1 csize(2)+1]);
skipLines(f,28);
info.frameAxis=readString(f);
% fprintf('frame axis : %s\n',info.frameAxis);
info.dataType=readString(f);
% fprintf('data type : %s\n',info.dataType);
info.imageAxis=readString(f);
% fprintf('image axis : %s\n',info.frameAxis);
o=fscanf(f,'%d',16);
info.imageArea=[o(2) o(5) o(7);o(4) o(3) o(6)];
info.frameArea=[o(11) o(14);o(13) o(12)]; % pixel frame 256x 256
info.frameBins=[o(16) o(15)];
s=(1 + diff(info.frameArea))./info.frameBins;
z=1 + diff(info.imageArea(5:6));
% if prod(s) ~= o(9) || o(9)*z ~= o(8)
%    fclose(f);
%    error('Inconsistent image header.');
% end
o=readString(f);
% fprintf('%s\n',o);
skipLines(f,z);
if subver > 25
    skipLines(f,z);
end
tline = fgetl(f);
% fprintf('last line ');
% fprintf('%s\n',tline);
% fprintf('%d %d %d\n', s, z);
info.imageData=reshape(fread(f,prod(s)*z,'*single'),[s z]);


function skipBytes(f,N)
[n]=fread(f,N,'*uint8');

function skipLines(f,N)
for n=1:N
   if isequal(fgetl(f),-1)
      fclose(f);
      error('Inconsistent image header.');
   end
end

function o=readString(f)
n=fscanf(f,'%f',1);  
if isempty(n) || n < 0 || isequal(fgetl(f),-1)
   fclose(f);
   error('Inconsistent string.');
end
o=fread(f,[1 n],'uint8=>char');