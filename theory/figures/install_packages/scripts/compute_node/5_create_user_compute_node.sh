#!/bin/sh

stty sane
echo /tmp/create_user_compute_node.log.$$
echo /tmp/create_user_compute_node.log.$$ > /tmp/create_user_compute_node.log.$$ 2>&1

USER_NAME=$1
USER_PASSWD=$2
GROUP_NAME=$3

# create user if not exists
id -u ${USER_NAME} >/dev/null 2>&1 || sudo adduser ${USER_NAME} --gecos "${USER_NAME},RoomNumber,WorkPhone,HomePhone" --disabled-password
echo "${USER_NAME}:${USER_PASSWD}" | sudo chpasswd

# create group if not exists
getent group ${GROUP_NAME} || sudo groupadd ${GROUP_NAME}
sudo usermod -a -G ${GROUP_NAME} ${USER_NAME}

# # add to sohr-cluster : account
# # cluster name: sohr
# # account name: sohr-account

# # add to sohr-cluster : account
# # in case cluster not exists
# # -i, commits change immediately

# sudo sacctmgr -i add cluster sohr-cluster
# sudo sacctmgr -i add account sohr-account description="Compute accounts" Organization=UCB
# sudo sacctmgr -i create user ${USER_NAME} account=sohr-account adminlevel=None

echo "user created successfull"
