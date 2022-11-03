#!/usr/bin/env Rscript

.libPaths(c("/root/miniconda3/envs/single-cell-test-Rbase/lib/R/library"))

args = commandArgs(trailingOnly=TRUE)

get_flg_arg <- function(flag, arg_list, var_type, mandetory){
  
  found_flag_num <- 0
  
  found_result <- ""
  
  define_type <- function(found_string, var_type, flag){
    return_var <- ""
    if (var_type == "Boolean"){
      tryCatch(
        {
          return_var <- as.logical(found_string)
        },
        error=function(e){
          stop(paste0("Failed to convert: ",flag, ": ", found_string, " to Boolean"))
        }
      )
    }else if(var_type == "String"){
      tryCatch(
        {
          return_var <- found_string
        },
        error=function(e){
          stop(paste0("Failed to convert: ",flag, ": ", found_string, " to String"))
        }
      )
    }else if(var_type == "Numeric"){
      tryCatch(
        {
          return_var <- as.numeric(found_string)
          if(is.na(return_var)){
            stop()
          }
        },
        error=function(e){
          stop(paste0("Failed to convert: ",flag, ": ", found_string, " to Numeric"))
        }
      )
      
    }else{
      stop(paste0("Undefined variable type: ", var_type))
    }
    
    return(return_var)
  }
  
  for (args in arg_list){
    
    parsed_input <- strsplit(grep(paste0("\\b",flag,"\\b*"), 
                                  args, value = TRUE), 
                             split = '=')
    
    if(length(parsed_input)!=0){
      found_string <- parsed_input[[1]][2]
      found_flag <- parsed_input[[1]][1]
      
      if(found_flag == flag){
        found_flag_num <- found_flag_num + 1
        found_result <- define_type(found_string, var_type, flag)
      }
    }
  }
  
  if(found_flag_num ==0){
    if (mandetory){
      stop(paste0("Flag: ", flag, " not found."))
    }else{
      return("")
    }
  }else if(found_flag_num >1){
    stop(paste0("Flag: ", flag, " found ", found_flag, " times."))
  }else{
    return(found_result)
  }
  
  
  
}

arg_job_name <- get_flg_arg("--job-name", args, 
                              var_type = "String", mandetory = "True")

arg_output_dir <- get_flg_arg("--output", args, 
                              var_type = "String", mandetory = "True")

#dir.create(file.path(getwd(), arg_output_dir), showWarnings = FALSE)

#arg_output_dir <- file.path(getwd(), arg_output_dir)

arg_Conserve_memory <- get_flg_arg("--Conserve_memory", args, 
                               var_type = "Boolean", mandetory = "True")

arg_Input_file <- get_flg_arg("--Input_file", args, 
                               var_type = "String", mandetory = "True")

arg_Number_of_components <- get_flg_arg("--Number_of_components", args, 
                              var_type = "Numeric", mandetory = "True")

arg_Exclude_sample <- get_flg_arg("--Exclude_sample", args, 
                                  var_type = "Numeric", mandetory = "True")

arg_tempdir_sample <- get_flg_arg("--tmp_dir", args, 
                                  var_type = "String", mandetory = "True")

tmpdir_loc=paste0(arg_tempdir_sample)
print(tmpdir_loc)
# Sys.setenv(TEMP=tmpdir_loc)
# Sys.setenv(TMPDIR=tmpdir_loc)
# Sys.setenv(TMP=tmpdir_loc)

print(tempdir())

arg_work_dir <- get_flg_arg("--work_dir", args, 
                            var_type = "String", mandetory = "True")

setwd(paste0(arg_work_dir))

source("/mnt/projects/CCBR-Pipelines/pipelines/Combine_Renorm_NIDAP/combine_and_renormalized.R")

sink(paste0(arg_output_dir,"/","combine_renorm_console_log.txt"))

RObjectdata <- readRDS(paste0(arg_Input_file))

combine_renorm_result <- combine_and_renormalized(
                                RObjectdata, 
                                Conserve_memory = arg_Conserve_memory, 
                                Number_of_components = arg_Number_of_components,
                                Exclude_sample = arg_Exclude_sample, 
                                image_out_dir = paste0(arg_output_dir)
                                )

saveRDS(combine_renorm_result,paste0(arg_output_dir,"/seurat_object.rds"))

tryCatch(
  {
    file.remove(paste0(arg_work_dir,"/Rplots.pdf"))
  },
  error=function(e){
    print(paste0("Failed to remove file: ",Rplots.pdf))
  }
)

print(sessionInfo())

sink()
