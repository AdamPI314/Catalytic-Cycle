#!/bin/sh

stty sane
echo /tmp/start_head.log.$$ > /tmp/set_up_cgroups.log.$$ 2>&1

sudo sed -i 's/^GRUB_CMDLINE_LINUX=.*$/GRUB_CMDLINE_LINUX="cgroup_enable=memory swapaccount=1"/' /etc/default/grub
sudo update-grub

sudo reboot
