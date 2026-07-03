# script to tidy taxonomic names and extract worms data if it exists

# get highest available level of taxonomy
master_long_df <- master_long_df %>%
  mutate(across(
    c(species_db, genus_db, family_db, order_db, class_db, phylum_db, kingdom_db),
    as.character
  )) %>%
  mutate(across(
    c(species_db, genus_db, family_db, order_db, class_db, phylum_db, kingdom_db),
    ~na_if(., "")
  )) %>%
  mutate(
    taxa_name = coalesce(
      species_db,
      genus_db,
      family_db,
      order_db,
      class_db,
      phylum_db,
      kingdom_db
    ),
    taxa_name = replace_na(taxa_name, "unclassified")
  )

# make new cleaned column
master_long_df <- master_long_df %>%
  mutate(
    taxa_name_clean = janitor::make_clean_names(taxa_name, allow_dupes = TRUE),
    taxa_name_clean = sub("_\\d+$", "", taxa_name_clean),
    species_db = janitor::make_clean_names(species_db, allow_dupes = TRUE) # clean up species_db names too
  )

# function to get metadata from worms
get_worms_record_safe <- function(sp_name) {
  
  rec <- try(
    wm_records_name(sp_name, fuzzy = FALSE),
    silent = TRUE
  )
  
  if (inherits(rec, "try-error") || is.null(rec) || nrow(rec) == 0) {
    return(tibble(
      input_name = sp_name,
      AphiaID_worms = NA_real_,
      species_worms = NA_character_,
      valid_name_worms = NA_character_,
      status_worms = NA_character_,
      rank_worms = NA_character_,
      citation_worms = NA_character_,
      isExtinct_worms = NA,
      kingdom_worms = NA_character_,
      phylum_worms = NA_character_,
      class_worms = NA_character_,
      order_worms = NA_character_,
      family_worms = NA_character_,
      genus_worms = NA_character_,
      marine = NA,
      brackish = NA,
      freshwater = NA,
      terrestrial = NA
    ))
  }
  
  rec <- rec[1, ]
  
  use_id <- ifelse(
    is.na(rec$valid_AphiaID),
    rec$AphiaID,
    rec$valid_AphiaID
  )
  
  meta <- try(
    wm_record(use_id),
    silent = TRUE
  )
  
  if (inherits(meta, "try-error") || is.null(meta)) {
    return(tibble(
      input_name = sp_name,
      AphiaID_worms = use_id,
      species_worms = NA_character_,
      valid_name_worms = NA_character_,
      status_worms = NA_character_,
      rank_worms = NA_character_,
      citation_worms = NA_character_,
      isExtinct_worms = NA,
      kingdom_worms = NA_character_,
      phylum_worms = NA_character_,
      class_worms = NA_character_,
      order_worms = NA_character_,
      family_worms = NA_character_,
      genus_worms = NA_character_,
      marine = NA,
      brackish = NA,
      freshwater = NA,
      terrestrial = NA
    ))
  }
  
  tibble(
    input_name        = sp_name,
    AphiaID_worms     = use_id,
    species_worms     = dplyr::case_when(meta$rank == "Species" ~ meta$scientificname, TRUE ~ NA_character_),
    valid_name_worms  = meta$valid_name,
    status_worms      = meta$status,
    rank_worms        = meta$rank,
    citation_worms    = meta$citation,
    isExtinct_worms   = as.logical(meta$isExtinct),
    kingdom_worms     = meta$kingdom,
    phylum_worms      = meta$phylum,
    class_worms       = meta$class,
    order_worms       = meta$order,
    family_worms      = meta$family,
    genus_worms       = meta$genus,
    marine            = as.logical(meta$isMarine),
    brackish          = as.logical(meta$isBrackish),
    freshwater        = as.logical(meta$isFreshwater),
    terrestrial       = as.logical(meta$isTerrestrial)
  )
}

# call function and safe new data
unique_species <- unique(master_long_df$taxa_name_clean)
worms_lookup <- map_df(unique_species, get_worms_record_safe)
master_long_worms_df <- master_long_df %>%
  left_join(worms_lookup, by = c("taxa_name_clean" = "input_name"))

# Create final taxonomy using available worms information, or db information if worms is not available
master_long_worms_df <- master_long_worms_df %>%
  mutate(
    across(
      .cols = where(is.character) & ends_with("_worms"),
      .fns = ~ na_if(.x, "")
    )
  ) %>%
  mutate(
    species = coalesce(species_worms, species_db),
    genus   = coalesce(genus_worms, genus_db),
    family  = coalesce(family_worms, family_db),
    order   = coalesce(order_worms, order_db),
    class   = coalesce(class_worms, class_db),
    phylum  = coalesce(phylum_worms, phylum_db),
    kingdom = coalesce(kingdom_worms, kingdom_db)
  )

# fix Actinopteri and petromyzonti issue - still occurs when worms is missing
master_long_worms_df$class[
  master_long_worms_df$class == "Actinopteri"
] <- "Teleostei"

master_long_worms_df$class[
  master_long_worms_df$class == "Petromyzonti"
] <- "Petromyzontida"

# get new tidied highest available taxonomy
master_long_worms_df <- master_long_worms_df %>%
  mutate(across(
    c(species, genus, family, order, class, phylum, kingdom),
    as.character
  )) %>%
  mutate(across(
    c(species, genus, family, order, class, phylum, kingdom),
    ~na_if(., "")
  )) %>%
  mutate(
    taxa_name_final = coalesce(
      species,
      genus,
      family,
      order,
      class,
      phylum,
      kingdom
    ),
    taxa_name_final = replace_na(taxa_name_final, "unclassified")
  )

# final clean so that taxa_names_final and species are consistent across merge
# make new cleaned column
master_long_worms_df <- master_long_worms_df %>%
  mutate(
    taxa_name_final = janitor::make_clean_names(taxa_name_final, allow_dupes = TRUE),
    species = janitor::make_clean_names(species, allow_dupes = TRUE)
  )

# write new
write.csv(master_long_worms_df, file = "Data/Processed/master_long_worms_df.csv")
