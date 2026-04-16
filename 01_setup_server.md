

## 1. create group called **master_b2s**

```bach
sudo groupadd master_b2s
```

## 2. create all users (x user) from txt file.

```bach
sudo adduser master_b2s  ## for one
```

**for multiple**
```bash
nano create_course_users.sh
```

```bash
#!/bin/bash

# Define the input file
INPUT_FILE="students.txt"

# 1. Check if run as root
if [ "$EUID" -ne 0 ]; then
  echo "Error: Please run this script with sudo."
  exit 1
fi

# 2. Check if file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: $INPUT_FILE not found!"
    exit 1
fi

echo "Starting user creation process..."
echo "---------------------------------"

# 3. Read the file, skipping the header row
tail -n +2 "$INPUT_FILE" | while IFS=$'\t' read -r username password fullname group || [ -n "$username" ]
do
    # Remove hidden Windows carriage returns (\r) just in case
    username=$(echo "$username" | tr -d '\r')
    password=$(echo "$password" | tr -d '\r')
    fullname=$(echo "$fullname" | tr -d '\r')
    group=$(echo "$group" | tr -d '\r')

    # Skip empty lines
    if [ -z "$username" ]; then
        continue
    fi

    # 4. Check if the course group exists, create it if it doesn't
    if ! getent group "$group" >/dev/null; then
        echo "📁 Group '$group' does not exist. Creating it now..."
        groupadd "$group"
    fi

    # 5. Check if the user already exists
    if id "$username" &>/dev/null; then
        echo "⚠️  User '$username' already exists. Skipping."
    else
        # 6. Create user: -m (home dir), -s (shell), -c (fullname), -G (add to shared course group)
        useradd -m -s /bin/bash -c "$fullname" -G "$group" "$username"
        
        # 7. Set the password safely (quotes ensure special characters like " and ; are protected)
        echo "$username:$password" | chpasswd
        
        # 8. Force the student to change this password on their first login
        chage -d 0 "$username"
        
        echo "Created user: $username | Added to group: $group"
    fi
done

echo "---------------------------------"
echo "User creation complete!"

```


```bash
chmod +x create_course_users.sh
sudo ./create_course_users.sh
```




metadata file look like

username password fullname group

## 3. add users to group

```bash
# Create the group
sudo groupadd master_b2s
sudo adduser yk01

# Add your existing students to the group (repeat for each user)
sudo usermod -aG master_b2s yk01
sudo usermod -aG master_b2s user01

```

## 4. Create the Directory Structure

```bash
# Change ownership to the shared group
#sudo chgrp -R students /courses/software /courses/miniconda3
mkdir -p /courses/{software,miniconda3,HH,master_b2s}
# Give them read and execute (traversal) rights
sudo chown root:root /courses
sudo chmod 711 /courses
sudo chmod -R 755 /courses/software
sudo chmod -R 755 /courses/miniconda3

```

```bash
# Create the base directory
sudo mkdir -p /courses/master_b2s/{00_scripts,references,raw_data}

# 1. Software directory
sudo mkdir -p /courses/master_b2s/00_scripts/{s01_qc,s02_vc,s03_scrna,s04_microbiom}

# 2. References directory
sudo mkdir -p /courses/master_b2s/references/{GRCh,silva,clinvar}

# 3. Raw Data directory
sudo mkdir -p /courses/master_b2s/raw_data/{s01_qc,s02_vc,s03_scrna,s04_microbiom}



```

## 5. Assign Ownership to the Group

```bash
sudo chown -R root:master_b2s /courses/master_b2s
```

## 6. Set the Collaborative Permissions

```bash

# Give full permissions (7) to Owner and Group, and Read/Execute (5) to Others
sudo chmod -R 755 /courses/master_b2s

# Apply the SGID bit to all directories inside the course folder
sudo find /courses/master_b2s -type d -exec chmod g+s {} +

```

## 7. move R library to courses

```bash
sudo mkdir -p /courses/software/R_libs
sudo chmod 755 /courses/software/R_libs
```

Tell R to use this folder globally:

- install R and Rstudio 
```bash
apt update
apt install build-essential -y

```
```bash
echo 'export R_LIBS_SITE="/courses/software/R_libs"' | sudo tee /etc/profile.d/course_r_libs.sh
```

To install new package do

```r
install.packages("ggplot2", lib="/courses/software/R_libs")

install.packages("Seurat", lib="/courses/software/R_libs")
```

install.packages(c("tidyverse", "Matrix", "RCurl", "scales", "cowplot", "patchwork", "devtools", "BiocManager", "hdf5r", "Seurat"), lib="/courses/software/R_libs")

\ user01 Wr;R0"/ user 01 master_b2s
user02	A!b@c#D$	user 02	master_b2s
where **Wr;R0"/** is password and **user 01** is full name



## 8. **Activate miniconda**
```bash
conda init
source ~/.bashrc

```

**Student**

```bash
cd /courses/miniconda3/bin/
./conda init
source ~/.bashrc
```




## 9. Create the environment with the tools we discussed

```bash
conda create -p /courses/miniconda3/envs/b2s_20260504 

conda activate /courses/miniconda3/envs/b2s_20260504

conda install -c bioconda -c conda-forge \
    fastqc multiqc fastp \
    bwa samtools bcftools \
    freebayes vcftools \
    gatk4 openjdk picard snpeff -y

sudo chmod -R 755 /courses/miniconda3

```



## 10. Software necessary

- FASTQC
- MultiQC
- GATK
- IGV
- Qiime2
- trimmomatic
- 


conda env create \
  --name qiime2-amplicon-2026.1 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml

conda env create \
  --name qiime2-moshpit-2026.1 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/moshpit/released/qiime2-moshpit-ubuntu-latest-conda.yml


conda env create \
  --name qiime2-pathogenome-2026.1 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/pathogenome/released/qiime2-pathogenome-ubuntu-latest-conda.yml