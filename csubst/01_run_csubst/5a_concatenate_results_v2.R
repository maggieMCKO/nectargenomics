library(tidyverse)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
#args 1 is where the dir with csubst output
#args 2 is interested run
#args 3 is number of cores

project_path = args[1]
# project_path = "/home/mpg08/mko/Nectar/analysis/csubst/01_run_csubst/out_v3/" 

interested_run = args[2]
# interested_run = "run3"

core = args[3]
# core = 1

concat_func = function(path){
  # tmp_path = target_files_run1[1]
  tmp_path = path 
  tmp_dir = dirname(tmp_path)
  tmp = unlist(strsplit(gsub(project_path, "", tmp_dir), "/"))
  tmp_transcript = tmp[2]
  tmp_run = tmp[3]
  tb = read_tsv(tmp_path, show_col_types = FALSE) 
  tb2 = bind_cols(transcript = tmp_transcript, run = tmp_run) %>% bind_cols(tb)
}


# 1. csubst_cb_stats.tsv ====
target_files = list.files(project_path, "csubst_cb_stats.tsv", full.names = TRUE, recursive = TRUE); # target_files
target_files_run1 = target_files[grepl(interested_run, target_files)]; # target_files_run1

out = mclapply(target_files_run1,  concat_func, mc.cores= core, mc.preschedule = FALSE) %>%
  bind_rows()

path = paste0(project_path, "/../run1_full_csubst_cb_stats_", Sys.Date(), ".tsv")
write_tsv(out, path, na = "")
rm(out)

# 2. csubst_b.tsv ====
target_files = list.files(project_path, "csubst_b.tsv", full.names = TRUE, recursive = TRUE); # target_files
target_files_run1 = target_files[grepl(interested_run, target_files)]; # target_files_run1

out = mclapply(target_files_run1,  concat_func, mc.cores= core, mc.preschedule = FALSE) %>%
  bind_rows()

path = paste0(project_path, "/../run1_full_csubst_b_", Sys.Date(), ".tsv")
write_tsv(out, path, na = "")

# 3. csubst_cb_2.tsv ====
target_files = list.files(project_path, "csubst_cb_2.tsv", full.names = TRUE, recursive = TRUE); # target_files
target_files_run1 = target_files[grepl(interested_run, target_files)]; # target_files_run1

out = mclapply(target_files_run1,  concat_func, mc.cores= core, mc.preschedule = FALSE) %>%
  bind_rows()

path = paste0(project_path, "/../run1_full_csubst_cb_2_", Sys.Date(), ".tsv")
write_tsv(out, path, na = "")
rm(out)


# 4. csubst_cb_3.tsv ====
target_files = list.files(project_path, "csubst_cb_3.tsv", full.names = TRUE, recursive = TRUE); # target_files
target_files_run1 = target_files[grepl(interested_run, target_files)]; # target_files_run1

out = mclapply(target_files_run1,  concat_func, mc.cores= core, mc.preschedule = FALSE) %>%
  bind_rows()

path = paste0(project_path, "/../run1_full_csubst_cb_3_", Sys.Date(), ".tsv")
write_tsv(out, path, na = "")
rm(out)


# 5. csubst_cb_4.tsv ====
target_files = list.files(project_path, "csubst_cb_4.tsv", full.names = TRUE, recursive = TRUE); # target_files
target_files_run1 = target_files[grepl(interested_run, target_files)]; # target_files_run1

out = mclapply(target_files_run1,  concat_func, mc.cores= core, mc.preschedule = FALSE) %>%
  bind_rows()

path = paste0(project_path, "/../run1_full_csubst_cb_4_", Sys.Date(), ".tsv")
write_tsv(out, path, na = "")
rm(out)

# 6. csubst_cb_5.tsv ====
target_files = list.files(project_path, "csubst_cb_5.tsv", full.names = TRUE, recursive = TRUE); # target_files
target_files_run1 = target_files[grepl(interested_run, target_files)]; # target_files_run1

out = mclapply(target_files_run1,  concat_func, mc.cores= core, mc.preschedule = FALSE) %>%
  bind_rows()

path = paste0(project_path, "/../output/", interested_run, "_full_csubst_cb_5_", Sys.Date(), ".tsv")
write_tsv(out, path, na = "")

# 7. csubst_cb_6.tsv ====
target_files = list.files(project_path, "csubst_cb_6.tsv", full.names = TRUE, recursive = TRUE); # target_files
target_files_run1 = target_files[grepl(interested_run, target_files)]; # target_files_run1

out = mclapply(target_files_run1,  concat_func, mc.cores= core, mc.preschedule = FALSE) %>%
  bind_rows()

path = paste0(project_path, "/../output/", interested_run, "_full_csubst_cb_6_", Sys.Date(), ".tsv")
write_tsv(out, path, na = "")

# 8. csubst_cb_7.tsv ====
target_files = list.files(project_path, "csubst_cb_7.tsv", full.names = TRUE, recursive = TRUE); # target_files
target_files_run1 = target_files[grepl(interested_run, target_files)]; # target_files_run1

out = mclapply(target_files_run1,  concat_func, mc.cores= core, mc.preschedule = FALSE) %>%
  bind_rows()

path = paste0(project_path, "/../output/", interested_run, "_full_csubst_cb_7_", Sys.Date(), ".tsv")
write_tsv(out, path, na = "")

# 9. csubst_cb_8.tsv ====
target_files = list.files(project_path, "csubst_cb_8.tsv", full.names = TRUE, recursive = TRUE); # target_files
target_files_run1 = target_files[grepl(interested_run, target_files)]; # target_files_run1

out = mclapply(target_files_run1,  concat_func, mc.cores= core, mc.preschedule = FALSE) %>%
  bind_rows()

path = paste0(project_path, "/../output/", interested_run, "_full_csubst_cb_8_", Sys.Date(), ".tsv")
write_tsv(out, path, na = "")

# 10. csubst_cb_9.tsv ====
target_files = list.files(project_path, "csubst_cb_9.tsv", full.names = TRUE, recursive = TRUE); # target_files
target_files_run1 = target_files[grepl(interested_run, target_files)]; # target_files_run1

out = mclapply(target_files_run1,  concat_func, mc.cores= core, mc.preschedule = FALSE) %>%
  bind_rows()

path = paste0(project_path, "/../output/", interested_run, "_full_csubst_cb_9_", Sys.Date(), ".tsv")
write_tsv(out, path, na = "")

