#!/bin/sh

# terminal setting
stty sane
echo /tmp/install_common.log.$$ > /tmp/install_common.log.$$ 2>&1

sudo apt-get -y update >> /tmp/install_common.log.$$ 2>&1
sudo apt-get -y upgrade >> /tmp/install_common.log.$$ 2>&1
sudo apt-get -y update >> /tmp/install_common.log.$$ 2>&1
sudo apt-get -y autoremove >> /tmp/install_common.log.$$ 2>&1

# tools sets
sudo apt-get -y install vim >> /tmp/install_common.log.$$ 2>&1
sudo apt-get -y install build-essential >> /tmp/install_common.log.$$ 2>&1
sudo apt-get -y install gfortran >> /tmp/install_common.log.$$ 2>&1
sudo apt-get -y install openssh-server >> /tmp/install_common.log.$$ 2>&1
# sudo apt-get -y install libmpich-dev >> /tmp/install_common.log.$$ 2>&1
# sudo apt-get -y install mpich >> /tmp/install_common.log.$$ 2>&1
sudo apt-get -y install openmpi-bin >> /tmp/install_common.log.$$ 2>&1
sudo apt-get -y install libboost-all-dev >> /tmp/install_common.log.$$ 2>&1

# Slurm
sudo apt-get -y install autoconf >> /tmp/install_common.log.$$ 2>&1
sudo apt-get -y install libffi-dev >> /tmp/install_common.log.$$ 2>&1
sudo apt-get -y update >> /tmp/install_common.log.$$ 2>&1
sudo apt-get -y install git gcc make ruby ruby-dev libpam0g-dev libmariadb-client-lgpl-dev libmysqlclient-dev >> /tmp/install_common.log.$$ 2>&1
sudo gem install fpm >> /tmp/install_common.log.$$ 2>&1

# Munge
sudo apt-get -y install libmunge-dev libmunge2 munge >> /tmp/install_common.log.$$ 2>&1
sudo systemctl enable munge >> /tmp/install_common.log.$$ 2>&1
sudo systemctl start munge >> /tmp/install_common.log.$$ 2>&1

sudo apt-get -y upgrade >> /tmp/install_common.log.$$ 2>&1
sudo apt-get -y update >> /tmp/install_common.log.$$ 2>&1
sudo apt-get -y autoremove >> /tmp/install_common.log.$$ 2>&1