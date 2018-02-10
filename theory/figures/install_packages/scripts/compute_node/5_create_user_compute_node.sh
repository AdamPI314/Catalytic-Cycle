#!/bin/sh

stty sane
echo /tmp/create_user_compute_node.log.$$
echo /tmp/create_user_compute_node.log.$$ > /tmp/create_user_compute_node.log.$$ 2>&1

USER_NAME=$1
USER_PASSWD=$2
GROUP_NAME=$3

# create user if not exists
id -u ${USER_NAME} &>/dev/null || sudo adduser ${USER_NAME} --gecos "${USER_NAME},RoomNumber,WorkPhone,HomePhone" --disabled-password
echo "${USER_NAME}:${USER_PASSWD}" | sudo chpasswd

# create group if not exists
getent group ${GROUP_NAME} || sudo groupadd ${GROUP_NAME}
sudo usermod -a -G ${GROUP_NAME} ${USER_NAME}

# add to compute-cluster : account

# in case cluster not exists
# sudo sacctmgr add cluster compute-cluster
# sudo sacctmgr add account compute-account description="Compute accounts" Organization=UCB
# sudo sacctmgr create user ${USER_NAME} account=compute-account adminlevel=None

echo "user created successfull"
