#!/bin/sh

stty sane
echo /tmp/slurm_conf.log.$$ > /tmp/slurm_conf.log.$$ 2>&1

slurm_user="slurm"
if ! id -u $slurm_user > /dev/null 2>&1; then
    sudo useradd slurm 
fi

dname1="/etc/slurm"
if ! [ -d $dname1 ]; then
	sudo mkdir $dname1 >> /tmp/slurm_conf.log.$$ 2>&1
fi

dname2dot5="/var/spool/slurm"
if ! [ -d $dname2dot5 ]; then
	sudo mkdir $dname2dot5 >> /tmp/slurm_conf.log.$$ 2>&1
fi

dname2="/var/spool/slurm/ctld"
if ! [ -d $dname2 ]; then
	sudo mkdir $dname2 >> /tmp/slurm_conf.log.$$ 2>&1
fi

dname3="/var/spool/slurm/d"
if ! [ -d $dname3 ]; then
	sudo mkdir $dname3 >> /tmp/slurm_conf.log.$$ 2>&1
fi

dname4="/var/log/slurm"
if ! [ -d $dname4 ]; then
	sudo mkdir $dname4 >> /tmp/slurm_conf.log.$$ 2>&1
fi

sudo chown slurm $dname1 $dname2 $dname3 $dname4 >> /tmp/slurm_conf.log.$$ 2>&1

# Assume a folder ~/Downloads/install_packages exists, there folder contains
# all the install scripts and source pakcages and configuration files 
# Assume ~/Downloads/install_packages/conf/slurmbd.service exists
sudo cp ~/Downloads/install_packages/conf/slurmdbd.service /etc/systemd/system/ >> /tmp/slurm_conf.log.$$ 2>&1
sudo cp ~/Downloads/install_packages/conf/slurmctld.service /etc/systemd/system/ >> /tmp/slurm_conf.log.$$ 2>&1

# Edit /storage/ubuntu-slurm/slurm.conf and replace AccountingStoragePass=slurmdbpass with
# AccountingStoragePass=/var/run/munge/munge.socket.2
sudo cp ~/Downloads/install_packages/conf/slurm.conf /etc/slurm/ >> /tmp/slurm_conf.log.$$ 2>&1

# Edit /storage/ubuntu-slurm/slurmdbd.conf and replace StoragePass=slrumdbpass with the DB password you used
# in the above SQL section.
sudo cp ~/Downloads/install_packages/conf/slurmdbd.conf /etc/slurm/ >> /tmp/slurm_conf.log.$$ 2>&1
