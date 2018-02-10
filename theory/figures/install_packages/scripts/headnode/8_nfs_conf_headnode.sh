#!/bin/sh

stty sane
echo /tmp/nfs_conf.log.$$ > /tmp/nfs_conf.log.$$ 2>&1
#sudo blkid to see UUID of disk/folder
#sudo vim /etc/exports, add two lines

str1="/ssd2t *(rw,sync)"
pattern1="\/ssd2t \*(rw,sync)"
if ! grep -q "${pattern1}" /etc/exports; then
    sudo echo $str1 >> /etc/exports | xargs >> /tmp/nfs_conf.log.$$ 2>&1
fi

str2="/hdd4t *(rw,sync)"
pattern2="\/hdd4t \*(rw,sync)"
if ! grep -q "${pattern2}" /etc/exports; then
    sudo echo $str2 >> /etc/exports | xargs >> /tmp/nfs_conf.log.$$ 2>&1
fi

# share /home directory
str3="/home *(rw,sync)"
pattern3="\/home \*(rw,sync)"
if ! grep -q "${pattern3}" /etc/exports; then
    sudo echo $str3 >> /etc/exports | xargs >> /tmp/nfs_conf.log.$$ 2>&1
fi

sudo service nfs-kernel-server restart >> /tmp/nfs_conf.log.$$ 2>&1



