library(dplyr)
library(tidyr)

# Create the dataframe
df <- data.frame(
  common_col_for_merging = c("ADENYYK-5", "ADENYYK-6", "CDREYMN-5", "PAZLSGTREE-6"),
  E2_A1_R1 = c(NA, NA, NA, NA),
  E2_A1_R2 = c(NA, NA, NA, NA),
  E2_A1_R3 = c(47, 47, 47, NA),
  E2_A2_R1 = c(NA, NA, NA, NA),
  E2_A2_R2 = c(63, NA, 63, NA),
  E2_A2_R3 = c(74, NA, 74, NA),
  E2_A3_R1 = c(NA, NA, NA, 11),
  E2_A3_R2 = c(14, 11, 14, 11),
  E2_A3_R3 = c(45, 45, 45, 45),
  E2_A4_R1 = c(77, 77, 77, 77),
  E2_A4_R2 = c(17, NA, NA, NA),
  E2_A4_R3 = c(14, 25, 14, 25),
  E2_A5_R1 = c(23, NA, 23, NA),
  E2_A5_R2 = c(NA, 14, NA, 25),
  E2_A5_R3 = c(NA, NA, 12, 47)
)

# Reshape the dataframe to long format
df_long <- df %>%
  pivot_longer(cols = starts_with("E2_A"),
               names_to = "col",
               values_to = "value")

# Count NA values based on prefixes and replace values
df_modified <- df_long %>%
  group_by(col = sub("^(E2_A[^_]*).*", "\\1", col)) %>%
  mutate(na_count = sum(is.na(value))) %>%
  # mutate(value = case_when(
  #   na_count == 3 ~ 0,
  #   na_count == 2 ~ 1,
  #   na_count == 1 ~ 2,
  #   TRUE ~ value
  # )) %>%
  ungroup() %>%
  select(-na_count) %>%
  pivot_wider(names_from = col, values_from = value)

# Combine the modified columns with the remaining columns
df_output <- bind_cols(df["common_col_for_merging"], df_modified)

# Print the modified dataframe
df_output




# Read the data into a data frame
df <- read.table(text = "common_col_for_merging	E2_A1_R1	E2_A1_R2	E2_A1_R3	E2_A2_R1	E2_A2_R2	E2_A2_R3	E2_A3_R1	E2_A3_R2	E2_A3_R3	E2_A4_R1	E2_A4_R2	E2_A4_R3	E2_A5_R1	E2_A5_R2	E2_A5_R3
ADENYYK-5	NA	NA	47	NA	63	74	NA	14	45	77	17	14	23	NA	NA
ADENYYK-6	NA	NA	47	NA	NA	NA	NA	11	45	77	NA	NA	25	NA	14
CDREYMN-5	NA	NA	47	NA	63	74	NA	14	45	77	NA	14	23	NA	12
PAZLSGTREE-6	NA	NA	NA	NA	NA	NA	11	11	45	77	NA	NA	25	NA	47", header = TRUE, stringsAsFactors = FALSE)


# Reshape the data from wide to long format
df_long <- df %>% 
  select(common_col_for_merging,all_of(experiment_name)) %>%
  pivot_longer(cols = -common_col_for_merging, names_to = "column", values_to = "value") %>%
  separate(column, into = c("col_group", "col_index", "row_index"), sep = "_")
  #mutate(col_index = str_remove(col_index, "A"))

# Count the number of NA values in each row for each column group
na_counts <- df_long %>%
  group_by(common_col_for_merging, col_index) %>%
  summarise(na_count = sum(is.na(value))) %>% 
  #pivot_wider(names_from = col_index, values_from = na_count)

# Replace the NA values based on the count
df_result <- df_long %>%
  left_join(na_counts, by = c("common_col_for_merging", "col_index")) %>%
  mutate(value = case_when(
    na_count == 3 ~ 0,
    
    na_count == 2 & col_index == "A1" ~ impute_values_2NA[1],
    na_count == 1 & col_index == "A1" ~ impute_values_1NA[1],
    
    na_count == 2 & col_index == "A2" ~ impute_values_2NA[2],
    na_count == 1 & col_index == "A2" ~ impute_values_1NA[2],
    
    na_count == 2 & col_index == "A3" ~ impute_values_2NA[3],
    na_count == 1 & col_index == "A3" ~ impute_values_1NA[3],
    
    na_count == 2 & col_index == "A4" ~ impute_values_2NA[4],
    na_count == 1 & col_index == "A4" ~ impute_values_1NA[4],
    
    na_count == 2 & col_index == "A5" ~ impute_values_2NA[5],
    na_count == 1 & col_index == "A5" ~ impute_values_1NA[5],
    TRUE ~ value
  )) %>%
  select(-na_count) %>% unite(new_col,col_index, row_index,sep = "_") %>%
  pivot_wider(names_from = new_col, values_from = value) %>%
  select(-col_group) 

# Print the resulting data frame
print(df_result)





























  