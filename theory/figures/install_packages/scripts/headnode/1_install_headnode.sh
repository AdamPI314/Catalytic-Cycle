#!/bin/sh

stty sane
echo /tmp/install_headnode.log.$$ > /tmp/install_headnode.log.$$ 2>&1

sudo apt-get -y update >> /tmp/install_headnode.log.$$ 2>&1
sudo apt-get -y upgrade >> /tmp/install_headnode.log.$$ 2>&1
sudo apt-get -y update >> /tmp/install_headnode.log.$$ 2>&1
sudo apt-get -y autoremove >> /tmp/install_headnode.log.$$ 2>&1

# Install MariaDB for Slurm accounting
sudo apt-get -y install mariadb-server >> /tmp/install_headnode.log.$$ 2>&1
sudo systemctl enable mysql >> /tmp/install_headnode.log.$$ 2>&1
sudo systemctl start mysql >> /tmp/install_headnode.log.$$ 2>&1

# Assume a folder ~/Downloads/install_packages exists, there folder contains
# all the install scripts and source pakcages and configuration files 
# Assume ~/Downloads/install_packages/slurm-17.02.9.tar.bz2 exists

cd ~/Downloads/install_packages/ >> /tmp/install_headnode.log.$$ 2>&1
tar xvjf slurm-17.02.9.tar.bz2 >> /tmp/install_headnode.log.$$ 2>&1
cd slurm-17.02.9 >> /tmp/install_headnode.log.$$ 2>&1
./configure --prefix=/usr --sysconfdir=/etc/slurm --enable-pam --with-pam_dir=/lib/x86_64-linux-gnu/security/ >> /tmp/install_headnode.log.$$ 2>&1
make >> /tmp/install_headnode.log.$$ 2>&1
make contrib >> /tmp/install_headnode.log.$$ 2>&1
sudo make install >> /tmp/install_headnode.log.$$ 2>&1

# NFS setup, head node
sudo apt-get -y install nfs-server >> /tmp/install_headnode.log.$$ 2>&1

sudo apt-get -y upgrade >> /tmp/install_headnode.log.$$ 2>&1
sudo apt-get -y update >> /tmp/install_headnode.log.$$ 2>&1
sudo apt-get -y autoremove >> /tmp/install_headnode.log.$$ 2>&1
