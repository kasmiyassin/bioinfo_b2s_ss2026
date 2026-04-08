# 1. Move to your software directory
cd /courses/master_b2s/software

# 2. Download the latest Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# 3. Install it into a specific path (-b for batch/silent, -p for path)
bash Miniconda3-latest-Linux-x86_64.sh -b -p /courses/master_b2s/software/miniconda3

# 4. Remove the installer file to save space
rm Miniconda3-latest-Linux-x86_64.sh

# Give the group read and execute permissions to the entire miniconda folder
sudo chmod -R 755 /courses/master_b2s/software/miniconda3

# 5. Initialize conda for bash and activate the base environment
source /courses/master_b2s/software/miniconda3/bin/activate
conda init bash
source ~/.bashrc

# 6. Update conda to the latest version
conda update conda

# 7. Create a new conda environment for QIIME 2 
mkdir -p /courses/master_b2s/software/qiime2
cd /courses/master_b2s/software/qiime2

## 7.1. QIIME 2 amplicon environment (2026.1)
### 7.1.1. Install the latest QIIME 2 amplicon environment (2026.1) using the provided YAML file
conda env create \
  --name qiime2-amplicon-2026.1 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml

### 7.1.2. Activate the newly created QIIME 2 environment and check the installation
conda activate qiime2-amplicon-2026.1
qiime info

### 7.1.3. Deactivate the newly created QIIME 2 environment
conda deactivate

## 7.2. QIIME 2 metagenomics environment (2026.1) 
### 7.2.1. Install the latest QIIME 2 metagenomics environment (2026.1) using the provided YAML file
conda env create \
  --name qiime2-moshpit-2026.1 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/moshpit/released/qiime2-moshpit-ubuntu-latest-conda.yml
### 7.2.2. Activate the newly created QIIME 2 environment and check the installation
conda activate qiime2-moshpit-2026.1
qiime info

### 7.2.3. Deactivate the newly created QIIME 2 environment
conda deactivate

## 7.3. QIIME 2 pathogenome environment (2026.1)
### 7.3.1. Install the latest QIIME 2 pathogenome environment (2026.1)
conda env create \
  --name qiime2-pathogenome-2026.1 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/pathogenome/released/qiime2-pathogenome-ubuntu-latest-conda.yml

### 7.3.2. Activate the newly created QIIME 2 environment and check the installation
conda activate qiime2-pathogenome-2026.1
qiime info

### 7.3.3. Deactivate the newly created QIIME 2 environment
conda deactivate




cat <<EOF | sudo tee /usr/local/bin/qiime-amplicon
#!/bin/bash
source /courses/master_b2s/software/miniconda3/bin/activate qiime2-amplicon-2026.1
exec qiime "\$@"
EOF

cat <<EOF | sudo tee /usr/local/bin/qiime-moshpit
#!/bin/bash
# Note: You named the env moshpit in step 7.2.1
source /courses/master_b2s/software/miniconda3/bin/activate qiime2-moshpit-2026.1
exec qiime "\$@"
EOF

cat <<EOF | sudo tee /usr/local/bin/qiime-pathogenome
#!/bin/bash
source /courses/master_b2s/software/miniconda3/bin/activate qiime2-pathogenome-2026.1
exec qiime "\$@"
EOF

sudo chmod +x /usr/local/bin/qiime-amplicon /usr/local/bin/qiime-moshpit /usr/local/bin/qiime-pathogenome


conda update -n base -c conda-forge conda -y


sudo chmod -R 775 /courses/master_b2s/software/miniconda3

# Students run this:
/courses/master_b2s/software/miniconda3/bin/conda init bash
conda env list
conda activate qiime2-amplicon-2026.1
qiime info