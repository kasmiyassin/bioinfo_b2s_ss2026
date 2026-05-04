## Before start wirte the following command

```bash
tmux new -A -s session1
```

# Exc 01

#### 01. `mkdir` – make directory create Directory

in chatgpt write **"Comment créer un dossier (directory) dans Ubuntu 24.04 via la ligne de commande ?"**


```bash
mkdir  session01/mygit    # To create directory  without -p
mkdir -p session01/mygit  # To create directory 

```
The `-p` flag creates all intermediate directories if they don't exist yet.  
Without `-p`, the command would fail if `session01/` does not already exist.


#### 02. `cd` – Change Directory

Use `cd` to move between directories.

```bash
cd                        # Go to your home directory (~)
cd session01/mygit        # Go into session01/mygit (relative path)
cd ..                     # Go one level up
cd /home                  # Go to /root (absolute path)
cd                     # Go  level
```

> 💡 **Tip:** Press `Tab` to auto-complete directory and file names.



#### 03. `ls` – List directory contents

```bash
ls                        # List files and folders in the current directory
ls /courses/master_b2s/raw_data/s01_qc/*.fastq.gz                   # List only files ending in .fq (FASTQ files)
ls /courses/master_b2s/ -l                     # Long format: shows permissions, size, date
ls -lh                    # Long format with human-readable file sizes
ls -a                     # Show hidden files (starting with .)
```

**Example output of `ls`:**

```
aboutme.md   data/   scripts/   README.md
```

**Example output of `ls *.fq`:**

```
sample1.fq   sample2.fq   sample3.fq
```

- What files are in your current directory?
- Are there any `.fq` files?



#### `nano` – Simple terminal text editor

```bash
cd session01/mygit
nano aboutme.md
```

This opens (or creates) the file `aboutme.md` in the `nano` editor.

| Shortcut | Action |
|----------|--------|
| `Ctrl + O` | Save the file |
| `Ctrl + X` | Exit nano |
| `Ctrl + K` | Cut a line |
| `Ctrl + U` | Paste a line |
| `Ctrl + W` | search |

> 📝 **Exercise:** Write a few lines about yourself in `aboutme.md` – name, background, why you are here.



### Viewing the directory tree

#### `tree` – Display directory structure

```bash
cd
tree
tree -d
```

Shows the current directory and all subdirectories as a visual tree.

**Example output:**

```
.
└── session01
    └── mygit
        └── aboutme.md

2 directories, 1 file
```

> If `tree` is not installed: `sudo apt install tree`


### Copying and moving files

#### `cp` – Copy files or directories

```bash
cp -r session01/mygit/aboutme.md session01/aboutme.md
```

| Flag | Meaning |
|------|---------|
| `-r` | Recursive – required when copying directories |

This copies `aboutme.md` from inside `mygit/` to the parent `session01/` folder.

#### `mv` – Move (or rename) files and directories

```bash
mv session01/mygit/  ~/mygit_old
```

This moves (renames) the `mygit` directory to `/root/mygit_old`.  
`mv` works for both files and directories – no `-r` flag needed.

> ⚠️ **Warning:** `mv` does not ask for confirmation. Double-check your paths!

---

## Check the HPC architucture

### 1. View Your Hard Drives and Partitions (The Disk Structure)

```bash
lsblk
```

What it tells you: It outputs a visual "tree" of all the physical hard drives attached to your server (usually named sda, nvme0n1, etc.) and the partitions inside them. Look at the SIZE column to find your biggest drives, and the MOUNTPOINTS column to see where they are attached (like /, /boot, or /data).

#### see space avaialble or free
To see how much space is actually used vs free on those mounted drives, use:

```bash
df -h
```

#### 2. Check Your Processors (CPU)
Variant calling with BWA and GATK relies heavily on multi-threading. 

```bash
lscpu

```

What to look for:

Look for the line that says CPU(s):. This is the total number of processing threads your server can handle simultaneously.


#### 3. Check Your Memory (RAM)
GATK is built on Java, which is notoriously hungry for RAM (you often have to set Java heap size using -Xmx8G). You need to know your total RAM to prevent the server from crashing due to "Out of Memory" errors.


```bash

free -h
```

What to look for: Look at the Mem: row and the total column. The -h makes it human-readable (e.g., 128Gi means 128 Gigabytes of RAM).

#### 4. Monitor the Server Live (The "Dashboard")
To see a live overview of your CPU, RAM, and what  are running in real-time, Ubuntu has a built-in tool called top, but I highly recommend installing its much better, colorful cousin, htop.

Install it first:

```bash
# sudo apt update
# sudo apt install htop
```
Then run it:

```bash
htop
```
What it does: You will see a live bar chart of every CPU core, your RAM usage, and a list of every active process. When you  start running  GATK commands, you will see their usernames pop up here along with exactly how much CPU and RAM their specific task is consuming. Press q to exit.


---

## Part 2 – Git & GitHub

### 1. Check your Git installation

```bash
git --version
```

Expected output (version may differ):

```
git version 2.43.0
```

If Git is not installed:

```bash
# sudo apt update && sudo apt install git -y
```

---

### 2. Clone a repository

First, navigate to (or create) your working directory:

```bash
cd | mkdir -p ~/session01/git
cd ~/session01/git
```

> The `|` here is used as shorthand to chain commands. Alternatively, run them on separate lines.

Now clone the course repository, checking out a specific branch:

```bash
git clone --branch session01_git https://github.com/kasmiyassin/bioinfo_b2s_ss2026.git
```

| Part | Meaning |
|------|---------|
| `git clone` | Download a repository to your local machine |
| `--branch session01_git` | Check out the `session01_git` branch |
| `https://github.com/...` | The URL of the remote repository |

After cloning, a new folder `bioinfo_b2s_ss2026/` will appear. Enter it:

```bash
cd bioinfo_b2s_ss2026
ls
```

```bash
mkdir -p ~/session01/git
cd ~/session01/git
git --version
git clone --branch session01_git https://github.com/kasmiyassin/bioinfo_b2s_ss2026.git
cd bioinfo_b2s_ss2026
ls
tree
```
