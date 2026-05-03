
```bash
# 1. Configuration
PREFIX="strato"
PASS="P1ssw0rd"
GROUP="master_b2s"

# 2. Create the group if it doesn't already exist
sudo groupadd -f $GROUP

# 3. The Loop
for i in $(seq -f "%02g" 1 25); do
    USERNAME="${PREFIX}${i}"
    
    # Check if user already exists to avoid "already exists" errors
    if id "$USERNAME" &>/dev/null; then
        echo "User $USERNAME already exists, skipping..."
    else
        # -m: create home dir
        # -g: set primary group
        # -s: set default shell to bash (so you get a nice prompt)
        sudo useradd -m -g "$GROUP" -s /bin/bash "$USERNAME"
        
        # Set the password
        echo "$USERNAME:$PASS" | sudo chpasswd
        
        echo "User $USERNAME created successfully."
    fi
done
```


remove

```bash
for i in $(seq -f "%02g" 1 25); do
    sudo userdel -r "strato$i"
done

```