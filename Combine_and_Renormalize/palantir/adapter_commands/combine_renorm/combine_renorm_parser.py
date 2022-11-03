
#import pathlib
import sys
import os
import logging

file_dir = os.path.dirname(os.path.abspath(__file__)) 
sys.path.append(os.path.join("..", file_dir))


default_config_file_name = "combine_renorm_config_default.cfg"
local_config_file_name = "combine_renorm_config_local.cfg"

from utils.general_utils import get_merged_config, common_prefix_list

def add_argument(attribute, value=None, new_line=True):
  """
  Add command line arguments to the pipeline command 
  @param attribute <str>
    The attribute to add "--" included 

  @param value <str>
    The value of the argument

  @param new_line <bool>
    Whether to add the character \ and a new line to the argument

  @return argument <str>
    The added argument
  """
 
  argument = "    {0}=".format(attribute)
  if value != None:
    argument += value
  if new_line:
    argument += " \\\n"
  return argument

def remove_argument(attribute, value=None, new_line=True):
  """
  Remove command line arguments to the pipeline command 
  @param attribute <str>
    The attribute to remove "--" included 

  @param value <str>
    The value of the argument

  @param new_line <bool>
    Whether to add the character \ and a new line to the argument

  @return argument <str>
    The added argument
  """

  argument = "    {0} ".format(attribute)
  if value != None:
    argument += value
  if new_line:
    argument += " \\\n"
  return argument


def setup_input_files(dataset_rid, transection_rid, file_name, input_dir, NIDAP_token, NDIAP_downloader_script):
  """
  Browse and Download data from NIDAP from selected dataset and transection

  @param NDIAP_downloader_script <str>
    The loacation of the NIDAP downloader script in bash

  @param dataset_rid <str>
    The rid of dataset
    
  @param transection_rid <str>
    The rid of transection
  
  @param file_name <str>
    The name of selected file

  @param job_dir <str>
    The path to the job work directory 

  @param NIDAP_token <str> 
    The DME authentication token

  @return download_to <list(str)>
    The path of the downloaded file
  """

  #Download input files from NIDAP
  import subprocess
  try:
    download_to = input_dir+"/"+str(file_name)
    script_args = [NIDAP_token, dataset_rid, transection_rid, str(file_name), download_to ]
    print(script_args)
    #cmd = [NDIAP_downloader_scriptscript_args]
    proc = subprocess.Popen(['/mnt/projects/CCBR-Pipelines/palantir/adapter_commands/utils/Download_from_NIDAP.sh',str(NIDAP_token),str(dataset_rid),str(transection_rid),str(file_name),str(download_to)], 
                            shell=False)
    proc.communicate()
  except subprocess.CalledProcessError as e:
    logging.error(e)
    
  return download_to


def process_combine_renorm(args):
  """
  Process the command line arguments to setup the RNA seek job to run on the HPC cluster

  @param args <argparse.Parser() object>
    The Namespace for the command options

  @return combine_renorm_cmd  <str>
    The full command for the RNA seek pipeline
  """
  import configparser

  config = configparser.ConfigParser()

  default_config_path =  os.path.join(file_dir, default_config_file_name)
  local_config_path =  os.path.join(file_dir, local_config_file_name)
  combine_renorm_config = get_merged_config(default_config_path, local_config_path)

  #Parse the user email to get the user name
  user_email = args.user_email
  user_name = user_email.split("@")[0]

  #Get the current datetime
  from datetime import datetime
  now = datetime.now()
  datetime_string = now.strftime("%Y-%m-%d-%H-%M-%S")

  #Setup the work directory
  top_work_dir = combine_renorm_config["general"]["workdir"]
  jobname = user_name + "--" + datetime_string
  job_dir = os.path.join(top_work_dir, jobname)
  os.mkdir(job_dir)
  logging.info("Creating job directory:" + job_dir)

  #Script line starter
  combine_renorm_command = combine_renorm_config["general"]["wrapper_script"]
  

  #Set up mode (--mode slurm)
  # mode = combine_renorm_config["general"]["mode"]
  # combine_renorm_command += add_argument("--mode", mode)

  #Set up the input directory
  input_dir = os.path.join(job_dir, "input")
  os.mkdir(input_dir)

  #Download and add input files
  NIDAP_token = args.NIDAP_token
  
  NIDAP_file_name = args.file_name
  
  NIDAP_input_dataset_rid = args.input_dataset_rid
  
  NIDAP_input_transaction_rid = args.transection_rid
  
  nidap_downloader_script=combine_renorm_config["general"]["nidap_downloader"]
  
  Input_file_dir = setup_input_files(NIDAP_input_dataset_rid, 
                                NIDAP_input_transaction_rid, 
                                NIDAP_file_name,
                                input_dir,
                                NIDAP_token,
                                nidap_downloader_script)

  #Set up the job name a1
  combine_renorm_command += add_argument("--job-name", jobname)
  
  #Set up the work dir a2
  combine_renorm_command += add_argument("--work_dir", top_work_dir)
  
  #set up environment image to use a3
  sif_cache = combine_renorm_config["general"]["sif_cache"]
  sif_local_line = "    " + sif_cache + "/" + args.environment + " \\\n"
  combine_renorm_command += sif_local_line
  
  #Set up the output directory (--output OUTPUT) a4
  output_dir  = os.path.join(job_dir, "output")
  os.mkdir(output_dir)
  combine_renorm_command += add_argument("--output", output_dir)

  #Set up the output_file dir a5
  ouput_dir_line = "    " + output_dir + " \\\n"
  combine_renorm_command += ouput_dir_line
  
  #Set up the NIDAP_token a6
  NIDAP_token_line = "    " + NIDAP_token + " \\\n"
  combine_renorm_command += NIDAP_token_line
  
  #Set up the NIDAP output dataset a7
  output_dataset_rid=args.output_dataset_rid
  NIDAP_output_data_rid_line = "    " + output_dataset_rid + " \\\n"
  combine_renorm_command += NIDAP_output_data_rid_line
  
  #Set up the tmp directory (--output tmp) a8
  tmp_dir = combine_renorm_config["general"]["tmp_dir"]
  tmp_dir_line = "    " + tmp_dir + " \\\n"
  combine_renorm_command += tmp_dir_line
  
  #Set up the tmp directory (--output tmp)
  combine_renorm_command += add_argument("--tmp_dir", tmp_dir)
  
  #Set up the input directory (--Input_file tmp)
  combine_renorm_command += add_argument("--Input_file", Input_file_dir)
  
  #Set up combine_renorm input arguements
  
  #set up sif-cache (--sif-cache SIF_CACHE)
  combine_renorm_command += add_argument("--Conserve_memory", args.Conserve_memory) 

  #set up sif-cache (--sif-cache SIF_CACHE)
  combine_renorm_command += add_argument("--Number_of_components", args.Number_of_components) 
  
  #set up sif-cache (--sif-cache SIF_CACHE)
  combine_renorm_command += add_argument("--Exclude_sample", args.Exclude_sample) 

  job_script_name = os.path.join(job_dir, "job-script.sh")
  with open(job_script_name, "w") as job_script:
    job_script.write(combine_renorm_command)

  logging.info("Writing the RNA seek script to:" + job_script_name )
  print(job_script_name)

  with open(args.output_location, "w") as output_loc:
    output_loc.write(output_dir)

  #Call the pipeline
  import subprocess
  try:
    #cmd = "bash {0}".format(job_script_name)
    #subprocess.run(cmd, shell=True)
    cmd = ["bash",job_script_name]
    proc = subprocess.Popen(cmd, shell=False)
    proc.communicate()
  except subprocess.CalledProcessError as e:
    logging.error(e)


def combine_renorm_subcommand(subparsers):
  """
  Add the command line options to setup the RNA-seek (combine_renorm) pipeline
  Only options that needs to process files are parse. Other pipeline options are passed down

  @param subparsers <argparse.SubParser() object>
    The subparser to add the combine_renorm options
  """

  combine_renorm = subparsers.add_parser("combine_renorm", help="Launch the RNA-seek pipeline")
  
  combine_renorm.add_argument('NIDAP_token', help="The NIDAP authentication token")
  combine_renorm.add_argument('output_location', help="A path of the file to communicate back where the output of the job should be")
  combine_renorm.add_argument('user_email', help="The email of the user who launched the job")

  combine_renorm.add_argument('--normal',help='The genome to search')
  combine_renorm.add_argument('--environment', help='The samples path on DME')
  combine_renorm.add_argument('--input_dataset_rid', help='The samples path on DME')
  combine_renorm.add_argument('--transection_rid', help='The samples path on DME')
  combine_renorm.add_argument('--file_name', help='The samples path on DME')
  combine_renorm.add_argument('--job-name', help='The samples path on DME')
  combine_renorm.add_argument('--Input_file', help='The samples path on DME')
  combine_renorm.add_argument('--output_dataset_rid', help='The samples path on DME')
  
   # Combine and renormalized function
   
  combine_renorm.add_argument('--Conserve_memory', help='The samples path on DME')
  combine_renorm.add_argument('--Number_of_components', help='The samples path on DME')
  combine_renorm.add_argument('--Exclude_sample', help='The samples path on DME')
