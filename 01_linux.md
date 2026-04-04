# Session 01 – Introduction to Linux & Git

> **Course:** Bioinformatics B2S SS2026  
> **Topics:** Linux (Ubuntu) command line basics · Git & GitHub  
> **Prerequisites:** None – bring a working terminal!

---

## Table of Contents

1. [Part 1 – Linux Basics](#part-1--linux-basics)
   - [Navigating the file system](#1-navigating-the-file-system)
   - [Listing files](#2-listing-files)
   - [Creating directories and files](#3-creating-directories-and-files)
   - [Viewing the directory tree](#4-viewing-the-directory-tree)
   - [Copying and moving files](#5-copying-and-moving-files)
2. [Part 2 – Git & GitHub](#part-2--git--github)
   - [Check your Git installation](#1-check-your-git-installation)
   - [Clone a repository](#2-clone-a-repository)
3. [Practical Exercises](#practical-exercises)
4. [Quick Reference Cheat Sheet](#quick-reference-cheat-sheet)

---

## Part 1 – Linux Basics

### 1. Navigating the file system

#### `mkdir` – create Directory

in chatgpt write **"comment cree un dossier (directory) dans ubunutu 24.04 par command"**


```bash
mkdir -p session01/mygit  # To create directory 
cd session01/mygit        # Go into session01/mygit (relative path)
cd /root                  # Go to /root (absolute path)
cd ..                     # Go one level up
cd                        # Go one level 0

```


Use `cd` to move between directories.

```bash
cd                        # Go to your home directory (~)
cd session01/mygit        # Go into session01/mygit (relative path)
cd /root                  # Go to /root (absolute path)
cd ..                     # Go one level up
```

#### `cd` – Change Directory

Use `cd` to move between directories.

```bash
cd                        # Go to your home directory (~)
cd session01/mygit        # Go into session01/mygit (relative path)
cd /root                  # Go to /root (absolute path)
cd ..                     # Go one level up
```

> 💡 **Tip:** Press `Tab` to auto-complete directory and file names.

---

### 2. Listing files

#### `ls` – List directory contents

```bash
ls                        # List files and folders in the current directory
ls *.fq                   # List only files ending in .fq (FASTQ files)
ls -l                     # Long format: shows permissions, size, date
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

---

### 3. Creating directories and files

#### `mkdir -p` – Make Directory (with parents)

```bash
mkdir -p session01/mygit
```

The `-p` flag creates all intermediate directories if they don't exist yet.  
Without `-p`, the command would fail if `session01/` does not already exist.

#### `nano` – Simple terminal text editor

```bash
nano aboutme.md
```

This opens (or creates) the file `aboutme.md` in the `nano` editor.

| Shortcut | Action |
|----------|--------|
| `Ctrl + O` | Save the file |
| `Ctrl + X` | Exit nano |
| `Ctrl + K` | Cut a line |
| `Ctrl + U` | Paste a line |

> 📝 **Exercise:** Write a few lines about yourself in `aboutme.md` – name, background, why you are here.

---

### 4. Viewing the directory tree

#### `tree` – Display directory structure

```bash
tree
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

---

### 5. Copying and moving files

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
mv session01/mygit/  /root/mygit_old
```

This moves (renames) the `mygit` directory to `/root/mygit_old`.  
`mv` works for both files and directories – no `-r` flag needed.

> ⚠️ **Warning:** `mv` does not ask for confirmation. Double-check your paths!

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
sudo apt update && sudo apt install git -y
```

---

### 2. Clone a repository

First, navigate to (or create) your working directory:

```bash
cd | mkdir -p /root/session01/git
cd /root/session01/git
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

---

## Practical Exercises

Work through these steps in order. Use the commands above as reference.

### Exercise 1 – Explore your environment

```bash
ls
ls *.fq
```

- What files are in your current directory?
- Are there any `.fq` files?

---

### Exercise 2 – Create your workspace

```bash
cd
mkdir -p session01/mygit
cd session01/mygit
nano aboutme.md
```

Write at least 3 lines about yourself. Save with `Ctrl+O`, exit with `Ctrl+X`.

---

### Exercise 3 – Check the structure

```bash
cd
tree
```

Confirm the directory tree looks like this:

```
.
└── session01
    └── mygit
        └── aboutme.md
```

---

### Exercise 4 – Copy and reorganise

```bash
cp -r session01/mygit/aboutme.md session01/aboutme.md
mv session01/mygit/ /root/mygit_old
tree
```

What changed in the directory structure?

---

### Exercise 5 – Clone the course repository

```bash
mkdir -p /root/session01/git
cd /root/session01/git
git --version
git clone --branch session01_git https://github.com/kasmiyassin/bioinfo_b2s_ss2026.git
cd bioinfo_b2s_ss2026
ls
tree
```

---

## Quick Reference Cheat Sheet

| Command | What it does |
|---------|-------------|
| `ls` | List files in current directory |
| `ls *.fq` | List files matching a pattern |
| `cd <dir>` | Change into a directory |
| `cd` | Go to home directory |
| `mkdir -p <path>` | Create directory (and parents) |
| `nano <file>` | Open/create a file in the nano editor |
| `tree` | Show directory structure as a tree |
| `cp -r <src> <dst>` | Copy file or directory |
| `mv <src> <dst>` | Move or rename a file or directory |
| `git --version` | Check Git is installed |
| `git clone --branch <branch> <url>` | Clone a specific branch of a repo |

---

*Happy coding! 🐧*