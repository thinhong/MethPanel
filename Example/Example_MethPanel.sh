### 1. Create sample config file
project="Example"
path="/path/to/${project}/raw/"
sample_config="/path/to/${project}/config/sample.${project}.pre.config"

echo -e "[ $project ]\n" > $sample_config
for i in `ls $path`; do
echo -e "\t[[ $i ]]" >> $sample_config;
echo >> $sample_config;
done

head $sample_config
# [ Example ]

# 	[[ 001-P1 ]]

# 	[[ 002-P1 ]]

# 	[[ 003-P1 ]]

# 	[[ 004-P1 ]]

### 2. Create system config file
system_config="/path/to/${project}/config/system.${project}.pre.config"
cat $system_config


### 3. Run MethPanel
# Load python version 3.7.4, java version jdk-13.33, bpipe version 0.9.9.9
module load python/3.7.4
module load java/jdk-13.33
module load bpipe/0.9.9.9

project="Example"
sample_config="/path/to/${project}/config/sample.${project}.pre.config"
system_config="/path/to/${project}/config/system.${project}.pre.config"
python "/path/to/pipe/run_Bpipe.py" $sample_config $system_config
