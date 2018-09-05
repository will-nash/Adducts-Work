library("RSQLite")
library("dplyr")


top10_intense <- function(Code, Ion_mode){

#Import the raw peak matrix
raw_peaks <- (read.csv(paste0("Z:/users/wxn264/Adducts_Data_Sets/", Code, "/", Ion_mode, ".csv")))

#Subset just the intensities and calculate the mean
average_intensity <- rowMeans(raw_peaks[,3:ncol(raw_peaks)], na.rm = T)

#Attach mean intensity to the peak matrix table
peaks_avg_intensity <- cbind(raw_peaks, average_intensity)

#Order the table by the mean intensity and only include the top 10% highest intensity peaks
ordered_intensity <- (peaks_avg_intensity[order(-average_intensity),])[1:(round(nrow(raw_peaks)/10)),]

#Take the ID for the top 10% highest intensity peaks
top_10 <- rownames(ordered_intensity)

#######################################################################################
### Importing of the significant peak differences ####################################

#Connects to the database with the significant feature pairs
db = dbConnect(SQLite(), "Z:/users/wxn264/Adducts_Data_Sets/results_0.8_0.05_20.sqlite")

#Generates a character vector for the right results file from the database
corr_results <- dbListTables(db)
corr_results <- corr_results[grep(Code, corr_results)]
corr_results <- corr_results[grep(Ion_mode, corr_results)]
corr_results <- corr_results[grep("counts", corr_results, invert = T)]

#Extracts the correct results file from the database and creates a dataframe
test <- dbReadTable(db, corr_results)

#Creates a blank dataframe for storage of the results
test5 <- data.frame()

#Loop that takes all rows of the results table where the feature ID for feature1 matches up to any of the feature ID's
# off of the top 10% list

for(T in top_10){

part <- subset(test, (test$Feature_1 == T))
test5 <- rbind(test5, part)

}

#Creates a blank data frame for the storage of results
test6 <- data.frame()

#Loop that takes all rows of the results table where Feature1 was filtered and filters by feature 2 for all rows which it matches up
# to any of the feature ID's off the top 10% list 

for(T in top_10){
  
  part <- subset(test5, (test5$Feature_2 == T))
  test6 <- rbind(test6, part)
  
}

#Saves the pairs which meet the top10 % criteria as a csv
write.csv(test6, paste0("Z:/users/wxn264/Adducts_Data_Sets/Heterogeneous_Dimers/Top_10%_pairs_", Code, "_", Ion_mode, ".csv"))

}

####################################################################################################################################
####################################################################################################################################
#Extracts the saved table of high intensity pairs for the dataset and determines which are due to homogenous or heterogeneous dimers

dimers_table <- function(Code, ion_mode){

#Imports the saved lists of pairs comprised of features in the top 10% highest intensity for that dataset
top10_pairs <- read.csv(paste0("Z:/users/wxn264/Adducts_Data_Sets/Heterogeneous_Dimers/Top_10%_pairs_", Code, "_", ion_mode, ".csv"))

#Creates and sorts a vector of the unique RT bins associated with each pair in the dataset
RT_bins_list <- sort(unique(top10_pairs$Bin))

output <- data.frame()

    for(RT in RT_bins_list){

      #Subsets the dataset by the value of the RT Bin
      RT_Bin <- subset(top10_pairs, top10_pairs$Bin == RT)
     
         #Adds or subtracts the mass of a proton from all features in the bin dependent on whether the dataset is in pos or neg mode
         if (ion_mode == "POS"){
         
           neutral_masses <- c((RT_Bin$mz_feature1 - 1.0078),(RT_Bin$mz_feature2 - 1.0078))
         
           } else { neutral_masses <- c((RT_Bin$mz_feature1 + 1.0078),(RT_Bin$mz_feature2 + 1.0078))} 

              #Calculate the minimum and maximum values for the mz difference between the two features assuming 10 ppm error
              for(mz in RT_Bin$mz_distance.1){
  
               min_dist <- mz - (mz/100000)
               max_dist <- mz + (mz/100000)
  
                   #checks if each monomer neutral  mass falls between the min and max error values calculated for each feature
                   # and saves the relelvant info in a table if a match is found
                   for(monomer in neutral_masses){
    
                         if (monomer <= max_dist & monomer >= min_dist){
      
                         print("YIPPEKIYAAAAYYY MOTHERDUCKERS")
                         result <- cbind(monomer, RT_Bin[RT_Bin$mz_distance.1 == mz,])
                         output <- rbind(output, result)
      
                         }    else print("FUUUUUUUUUUUCCCCK")
                   }       
              }
  
        }

#removes duplication of results in the table that come from repetition of features/neutral masses overlapping RT bin values
output <- distinct(output)

#Saves the output to a csv. file
write.csv(output, paste0("Z:/users/wxn264/Adducts_Data_Sets/Heterogeneous_Dimers/Dimers", Code, "_", ion_mode, ".csv"))

}

dimers_table("PCB12", "POS")
dimers_table("PCB2", "POS")
dimers_table("ST236", "POS")
dimers_table("PCB3", "NEG")
dimers_table("PCB5", "NEG")
dimers_table("MTBLS372", "NEG")

#####################################################################################################################
######################################################################################################################
############ Factoring in BEAMS information  ########################################################################









