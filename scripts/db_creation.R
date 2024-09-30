# set seed
set.seed(240930)

#libraries
library(stringr)
library(googledrive)

#data import
key <- read.delim('input/scan_model_mouse_rx.csv', sep = ";")
url <- "https://drive.google.com/drive/folders/1zh1UBNfhjiwsUA-zdnEagEuR2XUDTG4v"
drive_deauth()
folder_id <- drive_get(as_id(url)) #get subfolders id
files <- drive_ls(path = folder_id) #find files in folder


#loop dirs and download files inside them
for (i in seq_along(files$name)) {
    i_dir = drive_ls(files[i, ], type = "csv") #list files
    #download files
    for (file_i in seq_along(i_dir$name)) {
        drive_download(
        as_id(i_dir$id[file_i]),
        path = str_c("./input/", i_dir$name[file_i])
        )}
    }

# list the measurements files
measurements <- list.files(path = "./input/", 
                        pattern = "Measurements*", 
                        full.names = TRUE)

db <- data.frame() # initialize the db

for (file in measurements) {
    csv_data <- read.csv(file) # Read the CSV file
    csv_data$Scan <- str_sub(basename(file), start= 14L, end = -5L) # add case
    db <- rbind(db, csv_data) # Append the data to the combined_data data frame
}


# Merge the 'db' and 'key' data frames by 'Scan' column
merged_data <- merge(db, key, by = "Scan")

write.csv(merged_data, file = "./input/db.csv") # write the  final db
file.remove(measurements) # delete all the intermediate files
