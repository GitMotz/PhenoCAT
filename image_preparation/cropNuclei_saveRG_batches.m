function cropNuclei_saveRG_batches(pathToData,metaData,pathToResults,channels,clims,perimeter,n_batches)

% cropNuclei_saveRG saves croppped images of good nuclei per
% screen plate in 2 colors (as specified by channels, first position = G, sec pos = R)
% illumination correction
% B/C adjustments (as specified by clims)

% splits output of one plate into n_batches, creates separate folders

[~,plate_name] = fileparts(pathToData);

% filter meta data for 1 plate
plate_expr = '^.*_CP(\d\d\d)-.*';
curr_plate = cellfun(@str2double, regexp(plate_name,plate_expr,'tokens'));
ix_plate = metaData(:,4) == curr_plate;
metaData_plate = metaData(ix_plate,[1 2 6 7]);
clear metaData
clear ix_plate

[site_index,~, site_GI] = unique(metaData_plate(:,[1 2 3]),'rows');
cellIDs_plate = metaData_plate(:,4);
clear metaData_plate

n_images = size(site_index,1);

% create n_batches containing equal number of random images
rand_num = randi(n_images,n_images,1);
if mod(n_images, n_batches) == 0
    im_per_batch = n_images/n_batches;
    my_edges = (1:n_batches)*im_per_batch;
    batch_ix = discretize(rand_num, [0 my_edges]);
    
else % if not dividable, add remaining images to last batch
    remainder = mod(n_images, n_batches);
    im_per_batch = (n_images-remainder)/n_batches;
    my_edges = (1:n_batches)*im_per_batch;
    my_edges(end) = my_edges(end)+remainder;
    batch_ix = discretize(rand_num, [0 my_edges]);
end


% create variables for all images
channel_expr = '^.*C0(\d)';
channel_red = cellfun(@str2double, regexp(channels{1,1},channel_expr,'tokens'));
channel_green = cellfun(@str2double, regexp(channels{1,2},channel_expr,'tokens'));

path_segmentation = fullfile(pathToData,'SEGMENTATION');
path_images = fullfile(pathToData,'TIFF');

% load data for illum correction
path_BATCH = fullfile(pathToData,'BATCH');
green_Illum_name = fullfile(path_BATCH,sprintf('Measurements_batch_illcor_channel%03d_zstack000.mat',channel_green));
red_Illum_name = fullfile(path_BATCH,sprintf('Measurements_batch_illcor_channel%03d_zstack000.mat',channel_red));

green_Illum = load(green_Illum_name);
red_Illum = load(red_Illum_name);

green_Mean = green_Illum.stat_values.mean;
green_Std = green_Illum.stat_values.std;
clear green_Illum
red_Mean = red_Illum.stat_values.mean;
red_Std = red_Illum.stat_values.std;
clear red_Illum

% loop over batches
for i = 1:max(batch_ix);
    
    % create folder
    batch_folder_name = sprintf('batch_%d',i);
    mkdir(pathToResults,batch_folder_name)
    
    % filter metadata
    curr_batch_ix = batch_ix == i;
    curr_row_GI = find(batch_ix == i);
    curr_site_index = site_index(curr_batch_ix,:);
    
    
    
    % loop over sites from meta data (= images)
    for n = 1:size(curr_site_index,1);
        
        current_image_meta = curr_site_index(n,:);
        
        % convert site index to image file name
        row_conversion = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O'};
        current_row = row_conversion(1,current_image_meta(1,1));
        current_col = sprintf('%02d',current_image_meta(1,2)); % field width for digits = 2
        if rem(current_image_meta(1,3),36) == 0
            current_site = '36';
        elseif rem(current_image_meta(1,3),36) > 0
            current_site = sprintf('%02d',rem(current_image_meta(1,3),36));
        end
        % build file name
        curr_file_name = sprintf('CN-CP%d_%s%s_T0001F0%sL01A01Z01C01',curr_plate,current_row{1,1},current_col,current_site);
        fprintf('Working on image %s\n',curr_file_name)
        
        % load nuclear segmentation
        curr_nuclearSegmentation = double(imread(fullfile(path_segmentation,sprintf('%s_SegmentedNuclei.png',curr_file_name))));
        
        % get cell IDs
        ix_cellIDs = site_GI == curr_row_GI(n);
        curr_cellID_list = cellIDs_plate(ix_cellIDs,1);
        
        % load and process images
        curr_green_name = sprintf('CN-CP%d_%s%s_T0001F0%sL01%s.png',curr_plate,current_row{1,1},current_col,current_site,channels{1,2});
        curr_red_name = sprintf('CN-CP%d_%s%s_T0001F0%sL01%s.png',curr_plate,current_row{1,1},current_col,current_site,channels{1,1});
        
        % check if images exist (isfile works in 2017b plus)
        %if isfile(fullfile(path_images,curr_green_name)) && isfile(fullfile(path_images,curr_red_name))
        if exist(fullfile(path_images,curr_green_name),'file')==2 && exist(fullfile(path_images,curr_red_name),'file')==2
            
            try
                greenImage = imread(fullfile(path_images,curr_green_name)); % original
                greenImage = IllumCorrect(greenImage,green_Mean,green_Std,1); % illum corrected
                greenImage = greenImage-clims(2,1); % background substract
                greenImage(greenImage > clims(2,2)) = clims(2,2); % max B/C
                
                redImage = imread(fullfile(path_images,curr_red_name)); % original
                redImage = IllumCorrect(redImage,red_Mean,red_Std,1); % illum corrected
                redImage = redImage-clims(1,1); % background substract
                redImage(redImage > clims(1,2)) = clims(1,2); % max B/C
                
                % loop over cell IDs per image
                for m = 1:size(curr_cellID_list,1);
                    curr_cellID = curr_cellID_list(m,1);
                    
                    % get nuclear segmentation
                    ix_curr_cellID = curr_nuclearSegmentation == curr_cellID;
                    curr_props = regionprops(ix_curr_cellID,'BoundingBox');
                    curr_LeftCorner = curr_props.BoundingBox([1 2]);
                    curr_dimensions = curr_props.BoundingBox([3 4]);
                    curr_rows = round([curr_LeftCorner(1):(curr_LeftCorner(1)+curr_dimensions(1))]);
                    curr_cols = round([curr_LeftCorner(2):(curr_LeftCorner(2)+curr_dimensions(2))]);
                    % crop ROI
                    ix_curr_ROI = ix_curr_cellID(curr_cols,curr_rows);
                    
                    % crop and rescale images
                    curr_green_image = greenImage(curr_cols,curr_rows);
                    curr_green_image = curr_green_image./clims(2,2);
                    curr_green_image = curr_green_image.*ix_curr_ROI;
                    
                    curr_red_image = redImage(curr_cols,curr_rows);
                    curr_red_image = curr_red_image./clims(1,2);
                    curr_red_image = curr_red_image.*ix_curr_ROI;
                    
                    % generate background
                    curr_RGB = zeros(perimeter,perimeter,3);
                    
                    % place image in center or crop image if too big
                    if size(curr_green_image,1) < perimeter && size(curr_green_image,2) < perimeter
                        row_startpoint = round(perimeter/2-size(curr_green_image,1)/2);
                        final_rows = row_startpoint:(row_startpoint + size(curr_green_image,1)-1);
                        col_startpoint = round(perimeter/2-size(curr_green_image,2)/2);
                        final_cols = col_startpoint:(col_startpoint + size(curr_green_image,2)-1);
                        
                        curr_RGB(final_rows,final_cols,1) = curr_red_image;
                        curr_RGB(final_rows,final_cols,2) = curr_green_image;
                        
                    elseif  size(curr_green_image,1) >= perimeter || size(curr_green_image,2) >= perimeter
                        % crop both images in x and y
                        rowsToCrop = round((size(curr_green_image,1)-perimeter)/2);
                        if rowsToCrop >=0
                            curr_green_image = curr_green_image((rowsToCrop+1:size(curr_green_image,1)-rowsToCrop),:);
                            curr_red_image = curr_red_image((rowsToCrop+1:size(curr_red_image,1)-rowsToCrop),:);
                        end
                        
                        colsToCrop = round((size(curr_green_image,2)-perimeter)/2);
                        if colsToCrop >=0
                            curr_green_image = curr_green_image(:,(colsToCrop+1:size(curr_green_image,2)-colsToCrop));
                            curr_red_image = curr_red_image(:,(colsToCrop+1:size(curr_red_image,2)-colsToCrop));
                        end
                        
                        % place in center
                        row_startpoint = round(perimeter/2-size(curr_green_image,1)/2);
                        if row_startpoint == 0 % in case image is exactly 256
                            row_startpoint = row_startpoint+1;
                        end
                        final_rows = row_startpoint:(row_startpoint + size(curr_green_image,1)-1);
                        
                        col_startpoint = round(perimeter/2-size(curr_green_image,2)/2);
                        if col_startpoint == 0 % in case image is exactly 256
                            col_startpoint = col_startpoint+1;
                        end
                        final_cols = col_startpoint:(col_startpoint + size(curr_green_image,2)-1);
                        
                        curr_RGB(final_rows,final_cols,1) = curr_red_image;
                        curr_RGB(final_rows,final_cols,2) = curr_green_image;
                        
                    end
                    
                    imwrite(curr_RGB,sprintf('%s/%s/%s_%03d.png',pathToResults,batch_folder_name,curr_file_name,curr_cellID),'png','BitDepth',16)
                    
                end
            
            catch
                % if error appears do nothing (writes no image)    
            end
        else % if images dont exist write error message in log file
            fid = fopen(fullfile(pathToResults,'logfile.txt'),'a+'); % a+ = opens or writes new file and adds error message to end
            fprintf(fid, 'ERROR: Image %s was not written.\n',curr_file_name);
            fclose(fid);
        end
    end
end

end


