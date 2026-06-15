# extracting information from WoRMS

# defining functions

## extracts unique identifiers for all taxon used in WoRMS called AphiaIDs -------------------------------
get_wormsid_noerror <- function(sp_name){
  wm_id <- try(taxize::get_wormsid(sp_name,
                                   accepted = TRUE, 
                                   searchtype= "scientific", 
                                   marine_only = FALSE,
                                   ask = FALSE,
                                   row = 1,
                                   message = FALSE),
               silent = TRUE)
  
  #remove end word to attempt matching genus when previous no match
  
  if(class(wm_id) == "try-error"){
    sp_name <- str_remove(sp_name, "(\\s+\\w+)") 
    wm_id <- try(taxize::get_wormsid(sp_name,
                                     accepted = TRUE, 
                                     searchtype= "scientific", 
                                     marine_only = FALSE,
                                     ask = FALSE,
                                     row = 1,
                                     message = FALSE),
                 silent = TRUE)
  } 
  
  #convert non-matched species to NA AphiaID
  
  if(class(wm_id) == "try-error"){
    wm_id <- NA
  } 
  tibble(sciname = sp_name, aphiaID = as.double(unlist(wm_id)))
}

## accesses and formats any taxonomic meta data based on Aphia IDs -------------------------------
get_wormsmeta <- function(aphia_input){
  
  #split into smaller chunks for wm_record() to work
  taxadf_split <- split(aphia_input, (seq(nrow(aphia_input))-1) %/% 50) #split into smaller groups for worms to work
  temp_df <- data.frame()
  
  #run wm_record() through split list
  for (i in taxadf_split) {
    taxa_temp <- wm_record(id = i$aphiaID)
    temp_df = rbind(temp_df, taxa_temp)
  }
  tibble(temp_df)
}

# get list of taxa
eDNA_taxa_unique_list <- unique(master_long_df$taxa_name_clean) #Get unique taxa for WoRMs input as list
eDNA_taxa_unique_df <- data.frame(taxonomy = tolower(unique(master_long_df$taxa_name_clean)), 
                                  stringsAsFactors = FALSE) ##Get unique taxa for WoRMs input as df

# get aphia ids
aphiaID_worms_output <- eDNA_taxa_unique_list %>%
  purrr::map(get_wormsid_noerror, .progress = TRUE) %>%
  enframe() %>%
  unnest(cols = everything()) %>%
  select(-name)

# format
aphiaID_worms_output$taxa_name_clean = eDNA_taxa_unique_df$taxa_name_clean #add previous taxonomy name for joining later (order remains the same)
aphiaID_worms_output$aphiaID = as.integer(aphiaID_worms_output$aphiaID) #convert to integer for joining

# add back to main dataset
join_alphiaID <-aphiaID_worms_output
colnames(join_alphiaID)[1] <- "taxa_name_clean"
master_long_df <- full_join(master_long_df, join_alphiaID, by = join_by(taxa_name_clean))

# get taxonomic metadata from worms
worms_meta_output <- get_wormsmeta(aphiaID_worms_output)
str(worms_meta_output)

# formatting
worms_meta_output_tidy <- worms_meta_output %>%
  dplyr::rename(aphiaID = AphiaID) %>%
  filter(!is.na(aphiaID)) %>% #remove NAs
  distinct() #remove duplicates

# add into main dataset
master_long_df <- right_join(master_long_df, worms_meta_output_tidy, by = "aphiaID") #worms meta
which(is.na(master_long_df$aphiaID), arr.ind=TRUE) #check NAs in IDs
print(any(duplicated(master_long_df))) #check for duplicates