#!/bin/bash

# distribute ssh key from headnode to other nodes

# Basic info
date > /tmp/munge_key_distribute.log.$$ 2>&1
whoami >> /tmp/munge_key_distribute.log.$$ 2>&1
echo $@ >> /tmp/munge_key_distribute.log.$$ 2>&1

# Parameters
MASTER_NAME="headnode"
#MASTER_IP="192.168.1.2"
MASTER_IP="128.138.143.154"

NUM_OF_VM=1

# array of worker node names, worker node ip
WORKER_NAME_A=("node1")

WORKER_NAME_FOR_SLURM="node1"

WORKER_IP_A=("128.138.143.223")

ADMIN_USERNAME="sohr"
ADMIN_PASSWORD="sohr666888314"

# start munged
sudo systemctl enable munge >> /tmp/munge_key_distribute.log.$$ 2>&1 
sudo systemctl restart munge >> /tmp/munge_key_distribute.log.$$ 2>&1
# also push munge key and slurm.conf to them
echo "Prepare the local copy of munge key" >> /tmp/munge_key_distribute.log.$$ 2>&1 

mungekey=/tmp/munge.key.$$
sudo cp -f /etc/munge/munge.key $mungekey
sudo chown $ADMIN_USERNAME $mungekey

echo "Start looping all workers" >> /tmp/munge_key_distribute.log.$$ 2>&1 

i=0
while [ $i -lt $NUM_OF_VM ]
do
   echo "SCP to ${WORKER_NAME_A[$i]}"  >> /tmp/munge_key_distribute.log.$$ 2>&1 
   sudo -u $ADMIN_USERNAME scp $mungekey $ADMIN_USERNAME@${WORKER_IP_A[$i]}:/tmp/munge.key >> /tmp/munge_key_distribute.log.$$ 2>&1 

   echo "Remote execute on ${WORKER_NAME_A[$i]}" >> /tmp/munge_key_distribute.log.$$ 2>&1 
   sudo -u $ADMIN_USERNAME ssh $ADMIN_USERNAME@${WORKER_IP_A[$i]} >> /tmp/munge_key_distribute.log.$$ 2>&1 <<'ENDSSH1'
      sudo cp -f /tmp/munge.key /etc/munge/munge.key
      sudo chown munge /etc/munge/munge.key
      sudo chgrp munge /etc/munge/munge.key
      sudo rm -f /tmp/munge.key
      sudo systemctl enable munge
      sudo systemctl restart munge
ENDSSH1
   i=`expr $i + 1`
done
rm -f $mungekey