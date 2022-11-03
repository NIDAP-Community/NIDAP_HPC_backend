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

arg_Input_file <- get_flg_arg("--Input_file", args, 
                               var_type = "String", mandetory = "True")

arg_tempdir_sample <- get_flg_arg("--tmp_dir", args, 
                                  var_type = "String", mandetory = "True")

# Filename

arg_Organism <- get_flg_arg("--Organism", args, 
                              var_type = "String", mandetory = "True")

arg_Legend_Dot_Size <- get_flg_arg("--Legend_Dot_Size", args, 
                                  var_type = "Numeric", mandetory = "True")

arg_Perform_fine_tuning <- get_flg_arg("--Perform_fine_tuning", args, 
                                  var_type = "Boolean", mandetory = "True")



tmpdir_loc=paste0(arg_tempdir_sample)
print(tmpdir_loc)
# Sys.setenv(TEMP=tmpdir_loc)
# Sys.setenv(TMPDIR=tmpdir_loc)
# Sys.setenv(TMP=tmpdir_loc)

print(tempdir())

arg_work_dir <- get_flg_arg("--work_dir", args, 
                            var_type = "String", mandetory = "True")

setwd(paste0(arg_work_dir))

source("/mnt/projects/CCBR-Pipelines/pipelines/Single_R/Single_R_function.R")

sink(paste0(arg_output_dir,"/","Single_R_console_log.txt"))

RObjectdata <- readRDS(paste0(arg_Input_file))

pipeline_function_result <- Single_R_function(
                                RObjectdata, 
                                Organism = arg_Organism, 
                                Legend_Dot_Size = arg_Legend_Dot_Size,
                                Perform_fine_tuning = arg_Perform_fine_tuning,
                                image_out_dir = arg_output_dir
                                )

saveRDS(pipeline_function_result,paste0(arg_output_dir,"/Single_R_output.rds"))

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
