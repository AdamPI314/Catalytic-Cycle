#!/bin/sh

stty sane
echo /tmp/create_user.log.$$
echo /tmp/create_user.log.$$ > /tmp/create_user.log.$$ 2>&1

USER_NAME=$1
USER_PASSWORD=$2
GROUP=$3

sudo adduser ${USER_NAME}
sudo usermod -a -G ${GROUP} ${USER_NAME}

# add to compute-cluster : account

# in case cluster not exists
sudo sacctmgr add cluster compute-cluster
sudo sacctmgr add account compute-account description="Compute accounts" Organization=UCB
sudo sacctmgr create user ${USER_NAME} account=compute-account adminlevel=None

echo "user created successfull"
