#!/bin/sh

stty sane
echo /tmp/slurm_conf_compute_node.log.$$ > /tmp/slurm_conf_compute_node.log.$$ 2>&1

slurm_user="slurm"
if ! id -u $slurm_user > /dev/null 2>&1; then
    sudo useradd slurm 
fi

dname1="/etc/slurm"
if ! [ -d $dname1 ]; then
	sudo mkdir $dname1 >> /tmp/slurm_conf_compute_node.log.$$ 2>&1
fi

dname2="/var/spool/slurm"
if ! [ -d $dname2 ]; then
	sudo mkdir $dname2 >> /tmp/slurm_conf_compute_node.log.$$ 2>&1
fi

dname3="/var/spool/slurm/d"
if ! [ -d $dname3 ]; then
	sudo mkdir $dname3 >> /tmp/slurm_conf_compute_node.log.$$ 2>&1
fi

sudo chown slurm $dname1 $dname2 $dname3 >> /tmp/slurm_conf_compute_node.log.$$ 2>&1

# Edit /storage/ubuntu-slurm/slurm.conf and replace AccountingStoragePass=slurmdbpass with
# AccountingStoragePass=/var/run/munge/munge.socket.2
sudo cp ~/Downloads/install_packages/conf/slurm.conf /etc/slurm/ >> /tmp/slurm_conf_compute_node.log.$$ 2>&1

# cgroups
cp ~/Downloads/install_packages/conf/cgroup.conf /etc/slurm/ >> /tmp/slurm_conf_compute_node.log.$$ 2>&1

# Assume a folder ~/Downloads/install_packages exists, there folder contains
# all the install scripts and source pakcages and configuration files 
# Assume ~/Downloads/install_packages/conf/slurmbd.service exists
sudo cp ~/Downloads/install_packages/conf/slurmd.service /etc/systemd/system/ >> /tmp/slurm_conf_compute_node.log.$$ 2>&1

