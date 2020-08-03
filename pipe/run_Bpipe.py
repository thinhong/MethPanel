from configobj import ConfigObj
import sys, os, re, argparse, subprocess

BASEDIR_full = os.path.dirname(os.path.realpath(__file__))
BASEDIR = BASEDIR_full.replace("/pipe", "")
error = ["ERROR", "do not exist", "No such file or directory", "Exception", "corrupted"]

# read in parameters
parser = argparse.ArgumentParser()
parser.add_argument("sample_configp", help="Path to sample config file", type=str)
parser.add_argument("system_configp", help="Path to system config file", type=str)
args = parser.parse_args()

# check parameters sample config 
print 'Sample config file is at ', args.sample_configp
if not os.path.exists(args.sample_configp): print "The sample config file at", args.sample_configp, " does not exist!" 
# check parameters system config
print 'System config file is at ', args.system_configp
if not os.path.exists(args.system_configp): print "The sample config file at", args.system_configp, " does not exist!" 
# check BASEDIR
print 'The BASEDIR is at ', BASEDIR

# read sample parameter file
sample_config = ConfigObj(args.sample_configp)
try:
	project = sample_config.keys()[0]
except:
	print "Please provide the project name in []!"
	sys.exit(0)

samples = [i for i in sample_config[project]]
if len(samples) < 1:
	print "Please provide the sample name in [[]]"
	sys.exit(0)

# read system parameter file
system_config = {}
for i in open(args.system_configp):
	if "=" in i: 
		system_config[i.strip().split("=")[0]] = i.strip().split("=")[1][1:-1]

# create directories
try:
	os.system("mkdir -p " + system_config["projectp"] + "/" + project)
	os.system("mkdir -p " + system_config["projectp"] + "/" + project + "/bigTable")
except:
	print "There is NO projectp path in the system config file (", system_configp, "). Please provide the path to place the project!"

# run the pipeline
PROJECT = open(system_config["projectp"] + "/" + project + "/bigTable/project.csv", "w")
PROJECT.write("Sample,Group\n")
print "\n=================================================== SAMPLE INFORMATION ===================================================="
for sample in samples:
	print "+++++ Sample is", sample
	samplep = system_config["projectp"] + "/" + project + "/raw/" + sample + "/"
	os.system("mkdir -p " + samplep)
	PROJECT.write(sample + "," + sample + "\n")
	for lane in sample_config[project][sample]:
		print "--- Lane is", lane

print "================================================ END OF SAMPLE INFORMATION =================================================\n\n"
PROJECT.close()

# add BASEDIR to system config file
system_config_new_path = system_config["projectp"] + "/" + project + "/config/" + os.path.basename(args.system_configp) + ".post"
OUT = open(system_config_new_path, "w")
for line in open(args.system_configp): OUT.write(line + "\n")
OUT.write("\n// amplicon coordinates after processed\n")
OUT.write("\nAMP='" + system_config["projectp"] + "/" + project + "/ref/" + system_config["panel_name"] + "_sorted_mapped.bed'\n")
OUT.write("\n// path to CpG sites of all amplicons\n")
OUT.write("\nFRAME='" + system_config["projectp"] + "/" + project + "/ref/" + system_config["build_name"] + "_" + system_config["panel_name"] + ".CpG.mapped.bed'\n")
OUT.write("\nBASEDIR=" + '"' + BASEDIR + '"' + "\n")
OUT.write("\n")
OUT.write("\n")
OUT.close()

# copy sample config file to output path
cmd = "cp " + args.sample_configp + " " + system_config["projectp"] + "/" + project + "/config/" + os.path.basename(args.sample_configp) + ".post"
os.system(cmd)

# copy Bpipe config file to output path
cmd = "cp " + BASEDIR_full + "/bpipe.config" + " " + system_config["projectp"] + "/" + project + "/config/"
os.system(cmd)


print "======================================================= MethPanel Pipeline ===================================================\n"
projectp = system_config["projectp"] + "/" + project 

cmd = "cd " + projectp + "; grep open /proc/self/limits; ulimit -n 4096; grep open /proc/self/limits; " + "bpipe" + " run " + "-p SYSTEM_CONFIG=" + system_config_new_path + " " + BASEDIR_full + "/main.MethPanel.bpipe "  + " raw/*/*_R*.fastq.gz"
print cmd + "\n"
os.system(cmd)

