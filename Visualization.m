

function [D,QA]=Visualization(D,P,is_disp,is_show,save_fig,is_video)
%%%% Changed by Yael in 24.02.20
%%%% this function generates and save the images that will be displayed on
%%%% the viewer

QA=[];
new_cols=P.new_cols;
new_rows=P.new_rows;
enface_col = 293;
I_all=D.Im.I_align_all;
y_pos=D.Pos.y_reg_sort;
y_max=max(abs(y_pos));
[num_lines,~,~]=size(I_all);
Idisp_all=nan(num_lines,new_rows,new_cols+enface_col,3);
ImageFrame_all=cell(num_lines,1);

if save_fig
    SavePath=[P.SavePathMain '\Visualization'];
    if~isfolder(SavePath);mkdir(SavePath);end
end
if is_video
    SavePathVideo=[P.SavePathMain '\Video'];
    if~isfolder(SavePathVideo);mkdir(SavePathVideo);end
end

for n=1:num_lines
    if is_disp==1;disp(['visualization ' num2str(n) '\' num2str(num_lines)]);end
    I=squeeze(I_all(n,:,:));
    I = imresize(I,[new_rows,new_cols],'lanczos2'); 
    Iline=LocationDisplay(y_pos(n),new_rows,n,num_lines,y_max);
    Iline=double(Iline);
    Iline=Iline/max(Iline(:));
    I_Iline_adj=AdjustDisplay(Iline);
    Irgb=cat(3,I,I,I);
    Idisp=[I_Iline_adj Irgb];
    if is_show==1
        figure(1);clf;imshow(Idisp);
        if save_fig==1
           fig1=gcf;
           saveas(fig1.Number,[SavePath '\frame_' num2str(n) '.jpg']); 
        end
    end
    Idisp_all(n,:,:,:)=Idisp;
    Idisp(Idisp<0)=0;Idisp(Idisp>1)=1;
    if is_video==1;ImageFrame_all{n}=im2frame(Idisp);end 
end
   if is_video==1
       VideoDisp(ImageFrame_all,SavePathVideo,P.scan_str);
   end
   D.Im.Idisp_all=Idisp_all;
end

 
function PosImage=LocationDisplay(y_location,rows,n_line,num_lines,y_max)
color = [255 0 0];
x_vec=[-1514 1510];
a=figure(1);clf;
a.Color = 'white';
b=axes;
b.Color='Black';
set(gcf, 'Position',[300,300, rows, rows]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);

% hold on;
xlim([x_vec(1) x_vec(2)]);
y_max=max(1540,y_max);
ylim([-y_max y_max]);
hold on
x = 0;
y = 0;
r_vec=[500 1500 3000]; % ETDRS
pi_vec=[0  pi/2 pi (pi/2+pi) 2*pi]+pi/4;
green = [0 180 0]/255;
for l=1:3
    r=r_vec(l);
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(xunit, yunit,'color',green,'LineWidth' ,1.5);
    %plot line
end
for p=1:length(pi_vec)
    x=[r_vec(1)* cos(pi_vec(p)), r_vec(end)* cos(pi_vec(p))];
    y = [r_vec(1) * sin(pi_vec(p)), r_vec(end) * sin(pi_vec(p))];
    hold on;plot(x,y,'color',green,'LineWidth' ,1.5);
end
line(x_vec,[y_location,y_location],'Color',color/255,'LineWidth' ,3);
title([num2str(n_line) '/' num2str(num_lines)],'fontsize',22);
F = getframe(gcf);
[PosImage, ~] = frame2im(F);
% PosImage_pad =  padarray(PosImage,[50 50],255,'both');
close 1;
% A=sum(PosImage,3);
% PosImage(A==765,:)=[];
end

 
function Ifinal=AdjustDisplay(Iline)
%%%params
gray_lev=180;
% [r_I,c_I,~] = size(Iline);
new_size = 338;

II = imresize(Iline,[new_size,new_size]);
pad_size = (size(Iline,1)-size(II,1))/2;
II = padarray(II,[pad_size 0],1);
I2 = II(:,:,1);
I3 = II(:,:,2);
I4 = II(:,:,3);

[r,c] = find(I2 & I3 & I4); % white frame pixels
Ifinal =II;
px_white_r = 77:100;
px_white_c = 148:204;
for i=1:length(r)
    if ~(ismember(c(i),px_white_c) && ismember(r(i),px_white_r))
        Ifinal(r(i),c(i),:)=gray_lev/255;
    end
end
cut_num = 16;
left_cut = 14;
Ifinal = Ifinal(:,left_cut+cut_num:end-cut_num,:);
end

function VideoDisp(ImageFrame_all,SavePath,scan_str)

N=length(ImageFrame_all);
file_path= [SavePath '\' scan_str '.avi'];

v = VideoWriter(file_path,'Uncompressed AVI');
v.FrameRate=1;

open(v);

for n=1:N
    frameI=ImageFrame_all{n};
    writeVideo(v,frameI);
end

close(v);
        
% old_cd=pwd;cd(vP.path);
% implay(vP.file);cd(old_cd);

end
